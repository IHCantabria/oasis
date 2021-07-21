
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include "Lines.hpp"
#include "../SEM_math/quadrule.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"


int Line::GetId()
{
	return id;
}

Line::Line(int incId, double incG, double incRhoW, double incFondo)
{
	id = incId;
	g = incG;
	rhoW = incRhoW;
	fondo = incFondo;
}


void Line::ReadPropertiesASCII(FILE* pFilePointer) 
{	
	// Declare variables
	char buffer_line [1000];
	double dtemp;

	//Ignoro las tres primeras lineas, donde pone "New line"
	for(int ii=0; ii<3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFilePointer);
	}

	//Leo todo
	fscanf(pFilePointer, "%d %[^\n]\n", &lineType, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &flag_tension, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &nNodos, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &p, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &L, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &rho0, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &d, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &flag_stiffness, buffer_line);
	if(flag_stiffness>1){
		strain_data = arma::zeros(flag_stiffness,1);
		stress_data = arma::zeros(flag_stiffness,1);
		for (int ii=0; ii<this->flag_stiffness; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1) 
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the strain data for line: " << this->GetId() <<"\n";
				throw ValueError(ss.str());
			}
			strain_data(ii,0) = dtemp;
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
		for (int ii=0; ii<this->flag_stiffness; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1) 
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the stress data for line: " << this->GetId() <<"\n";
				throw ValueError(ss.str());
			}
			stress_data(ii,0) = dtemp;
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
	} else {
		fscanf(pFilePointer, "%lf %[^\n]\n", &EA, buffer_line);
	}
	fscanf(pFilePointer, "%lf %[^\n]\n", &beta, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &CB, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Cmn, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Cdn, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Cdt, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &GK, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &GC, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Gmu, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Gvc, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Dz, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &BCP_N, buffer_line);
	BCP_N -= 1;
	indexBcps[1] = BCP_N;
	fscanf(pFilePointer, "%d %[^\n]\n", &BCP_1, buffer_line);
	BCP_1 -= 1;
	indexBcps[0] = BCP_1;

	if (lineType==1){
		floor_flag = 1;
	}
	fscanf(pFilePointer, "%d %[^\n]\n", &smoothstep, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &frictionModel, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &vth, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &us, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &ud, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &deltamax, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &vstatic, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &vdynamic, buffer_line);

	A = arma::datum::pi*d*d*0.25;
	dL = L/(nNodos-1);
	dL0 = dL;
	N = p*(nNodos-1)+1;
	Kn = Cmn * A * rhoW;

	dampCoef = 2.0 * sqrt(rho0 * GK * d ) ;
	VR= 0.01*(d*d*GK)/(dampCoef);


	pos = arma::zeros(3*N);
	vel = arma::zeros(N,3);
	acc = arma::zeros(N,3);
	F = arma::zeros(N,3);
	s = arma::zeros(N);
	xc = arma::zeros(N);
	zc = arma::zeros(N);
	dxcds = arma::zeros(N);
	dzcds = arma::zeros(N);
	Te = arma::zeros(N);
	roots = arma::zeros(p+1);
	weights = arma::zeros(p+1);
	FF = arma::zeros(N,3);
	ff = arma::zeros(N,3);
	t = arma::zeros(N,3);
	e_z = arma::zeros(1,3);
	e_z(0,2) = 1.0;
	posFriccion = arma::zeros(N,2);
	isSlip = arma::zeros(N,1);
	ffslid=arma::zeros(1,2);
	ffstick=arma::zeros(1,2);

	if (frictionModel==1){
		//voy a resolver un solo sistema que dara los coeficientes por si hay friccion de velocidad al principio del problemna.
		//es el step fuerte
		//Los coeficientes van a ser siempre esos asi que solo hace falta hacerlo una vez y asi se obtiene el polinomio
		arma::mat tmpA = arma::zeros(4,4);
		arma::mat tmpB = arma::zeros(4,1);
		double x1= 0;
		double x2= 1e-04;
		tmpA = {{  pow(x1,3), pow(x1,2), 1*x1, 1},
	       		{3*pow(x1,2),    2*(x1),    1, 0},
				{    pow(x2,3),   pow(x2,2),   1*x2, 1},
				{  3*pow(x2,2),        2*x2,    1, 0}};
		tmpB(0,0) = 1e-05; tmpB(1,0) = 0.0; tmpB(2,0) = x2; tmpB(3,0) = 1;
		std::cout << "tmpA="<< std::endl << tmpA<< std::endl;
		std::cout << "tmpB="<< std::endl << tmpB<< std::endl;
		a_1 = arma::solve(tmpA,tmpB);
		std::cout << "a1="<< std::endl << a_1<< std::endl;

	}


	double * roots_temp   = new double[p+1];
	double * weights_temp = new double[p+1];

	lobatto_set(p+1,roots_temp,weights_temp);

	for(int ii=0; ii<p+1; ii++)
	{
		roots(ii) = roots_temp[ii];
		weights(ii) = weights_temp[ii];
	}

	int kk;
	for(int ii=0; ii<N; ii++){
		kk = ii % p;
		s(ii,0) = dL * ( (ii-kk)/p + ( roots(kk) + 1.0 ) * 0.5 );
	}
}


void Line::SEM_getBaseFunctions(void)
{

	pos.reshape(3,N);
	pos = pos.t();

	this->SEM_coefficients();

	arma::mat D_local = SEM_get_D_local();
	D = arma::zeros(N,N);
	for(int ii=0;ii<nNodos-1;ii=ii+1){
		D.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) = D.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) + D_local;
	}
	for(int ii=1;ii<nNodos-1;ii=ii+1){
		D.row(ii*p) = D.row(ii*p) * 0.5;
	}

	int nIntegrate = p + 2;
	double * roots_temp   = new double[nIntegrate+1];
	double * weights_temp = new double[nIntegrate+1];
	lobatto_set(nIntegrate+1,roots_temp,weights_temp);

	arma::mat MassMatrix_local = arma::zeros(p+1,p+1);
	arma::mat StiffMatrix_local = arma::zeros(p+1,p+1);
	arma::mat MSMatrix_local = arma::zeros(p+1,p+1);

	for(int i=0;i<p+1;i=i+1){
		for(int j=0;j<p+1;j=j+1){
			for(int k=0;k<nIntegrate+1;k=k+1){
				MassMatrix_local(i,j) = MassMatrix_local(i,j) + weights_temp[k]*this->SEM_poly(roots_temp[k],i)*this->SEM_poly(roots_temp[k],j);
				StiffMatrix_local(i,j) = StiffMatrix_local(i,j) + weights_temp[k]*this->SEM_poly_first_derivative(roots_temp[k],i)*this->SEM_poly_first_derivative(roots_temp[k],j);
				MSMatrix_local(i,j) = MSMatrix_local(i,j) + weights_temp[k]*this->SEM_poly_first_derivative(roots_temp[k],i)*this->SEM_poly(roots_temp[k],j);
			}		
		}
	}

	MassMatrix = arma::zeros(N,N);
	StiffMatrix = arma::zeros(N,N);
	MSMatrix = arma::zeros(N,N);

	for(int ii=0;ii<nNodos-1;ii=ii+1){
		MassMatrix.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) = MassMatrix.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) + MassMatrix_local;
		StiffMatrix.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) = StiffMatrix.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) + StiffMatrix_local;
		MSMatrix.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) = MSMatrix.submat(ii*p,ii*p,(ii+1)*p,(ii+1)*p) + MSMatrix_local;
	}

	MM = 0.5*(rho0 + Kn) * MassMatrix;
	inv_MM = arma::solve(MM,arma::eye(N,N));
	MM.row(0) = arma::zeros(1,N);
	MM(0,0) = 1.0;
	inv_MM_1 = arma::solve(MM,arma::eye(N,N));
	MM.row(N-1) = arma::zeros(1,N);
	MM(N-1,N-1) = 1.0;
	inv_MM_1N = arma::solve(MM,arma::eye(N,N));

	MM = 0.5*(rho0 + Kn) * MassMatrix;
	MM.row(N-1) = arma::zeros(1,N);
	MM(N-1,N-1) = 1.0;
	inv_MM_N = arma::solve(MM,arma::eye(N,N));

	MM = 0.5 * (rho0 + Kn) * MassMatrix;

	D_sp = arma::sp_mat(D);
	MassMatrix_sp = arma::sp_mat(MassMatrix);
	StiffMatrix_sp = arma::sp_mat(StiffMatrix);
	MSMatrix_sp = arma::sp_mat(MSMatrix);
	MM_sp = arma::sp_mat(MM);
}


void Line::SEM_coefficients(void)
{

	C = arma::zeros(p+1,p+1);
	arma::mat M = arma::zeros(p+1,p+1);
	for(int ii=0;ii<p+1;ii=ii+1){
		M.insert_cols(ii,arma::pow(roots,ii));
	}
	arma::mat I = arma::eye(p+1,p+1);
	C = arma::solve(M,I);
	C = C.t();
}


double Line::SEM_poly(double x, int i)
{
	double y = 0.0;
	arma::mat cc = C.row(i);
	for(int ii=0;ii<p+1;ii=ii+1){
		y = y + cc(ii)*pow(x,ii);
	}
	return y;
}


double Line::SEM_poly_first_derivative(double x, int i)
{
	double y = 0.0;
	arma::mat cc = C.row(i);
	for(int ii=1;ii<p+1;ii=ii+1){
		y = y + ii*cc(ii)*pow(x,ii-1);
	}
	return y;
}


arma::mat Line::SEM_get_D_local(void)
{
	arma::mat D = arma::zeros(p+1,p+1);
	for(int ii=0;ii<p+1;ii=ii+1){
		for(int jj=0;jj<p+1;jj=jj+1){
			D(ii,jj) = this->SEM_poly_first_derivative(roots(ii),jj);
		}
	}
	return D;
}


void Line::SEM_computeF(void)
{
	//drds = (D_sp * pos) * (2.0/dL0);
	//drdsdt = (D_sp * vel) * (2.0/dL0);

	drds = (D * pos) * (2.0/dL0);
	drdsdt = (D * vel) * (2.0/dL0);

	norm_drds = sqrt(pow(drds.col(0),2) + pow(drds.col(1),2) + pow(drds.col(2),2));
	arma::mat strain = 0.5*((drds.col(0) % drds.col(0) + drds.col(1) % drds.col(1) + drds.col(2) % drds.col(2)) - dL/dL0);

	dedt = drds.col(0) % drdsdt.col(0) + drds.col(1) % drdsdt.col(1) + drds.col(2) % drdsdt.col(2);

	if(flag_stiffness>1){
		arma::mat temp_strain = norm_drds - dL/dL0;
		T = interp1(strain_data,stress_data,temp_strain);
		arma::mat temp_EA = T/temp_strain;
		T = T + temp_EA % (beta * dedt);
	} else {
		T = EA * (norm_drds - dL/dL0 + beta * dedt);
		//T = EA * (strain + beta * dedt);
	}

	if ((flag_tension == 2) && (flag_stiffness<=1)){
		T = 0.5*(T + arma::abs(T));
		double T0 = 0.1; // size of the smoothing region in newtons
		double a2 = 2.0/T0, a3 = -1.0/(T0*T0);
		arma::mat temp_T = a3*arma::pow(T,3)+a2*arma::pow(T,3);
		arma::umat ind = arma::find(((T>0.0)&&(T<T0)));
		T.elem(ind) = temp_T.elem(ind);
	}
	posFriccion= pos.cols(0,1);
	for(int k=0;k<N;k=k+1){
		t.row(k) = drds.row(k) / norm_drds(k);
		FF.row(k) = T(k) * t.row(k);

		fg = (rho0 - rhoW * A) * g / norm_drds(k);
		ff.row(k) = -fg * e_z;

		v = vel.row(k);
		vt = (v * t.row(k).t())* t.row(k);
		vn = v - vt;
		ff.row(k) = ff.row(k) - 0.5 * Cdt * d * rhoW * arma::norm(vt,2) * vt;
		ff.row(k) = ff.row(k) - 0.5 * Cdn * d * rhoW * arma::norm(vn,2) * vn;

		if (floor_flag == 1){
	
			if (smoothstep == 1){
				paramNormal= step(fondo - pos(k,2), -d/2, 0, 0, 1);
				paramVel= step(v(2), -VR, 1, 0, 0);
				parammuelle1= step(fondo - pos(k,2), 0, 0, d/2, 1);
		    	parammuelle2= GK*d*(fondo-pos(k,2)) ;
				ff(k,2)= ff(k,2) + fg*paramNormal + parammuelle1*parammuelle2 - GC*paramNormal*dampCoef*paramVel*v(2);
			} else {

				// comparativa con palm
				// if ((fondo-pos(k,2))>= 0) {
				// 	ff(k,2)= ff(k,2) + fg +  GK*d*(fondo-pos(k,2)) - GC*dampCoef*(std::min(v(2),0.0));
				// }

				//comparativa con el anterior
				fs = abs(fg) * exp(- GK * d * (pos(k,2) - fondo)/abs(fg));
			    GC = 2.0 * sqrt(rho0 * GK * d ) / (abs(fg) * d);
				fd= fs* GC*d *pow(std::min(v(2), 0.0), 2);
				ff(k,2)=ff(k,2) + fs + fd;

			}
			if (frictionModel==2) {
				arma::mat vxy= arma::zeros(1,2);
				vxy(0,0)=v(0);
				vxy(0,1)=v(1);
				double normvxy= arma::norm(vxy, 2);
				fn=fg*step(fondo - pos(k,2), -d/2, 0, 0, 1); //se ha suavizado
				if (isSlip(k,0)==1){
					if(normvxy > vth){
						isSlip(k,0)=0;
						posFriccion(k,0)=pos(k,0);
						posFriccion(k,1)=pos(k,1);
					}
				}
				if(isSlip(k,0)==0){
					if (normvxy <= vth){
						isSlip(k,0)=1;
						//no guarda posicion, ya que solo la guarda si pasa de split a slip
					}
				}
				arma:: mat posActual= arma::zeros(1,2);
				arma:: mat deformacion = arma::zeros(1,2);
				posActual(0,0)=pos(k,0);
				posActual(0,1)=pos(k,1);
				deformacion= posActual - posFriccion.row(k);
				double delta= arma::norm(deformacion, 2);
				double betaStick, usstep, udstep;

				// parametros distintos si es stick o slip
				if(isSlip(k,0) == 0){
					betaStick=1;
					usstep=0;
					udstep=ud;

				}else{
					betaStick=step(normvxy, -vth, -1.0, vth, 1.0);
					usstep=step(delta, -deltamax, -us, deltamax, us);
					udstep=step(normvxy, -vth, -ud, vth, ud);
				}
				//hay que calcular las fuerzas
				ffslid= -vxy* udstep* fn/std::max(normvxy,1e-7);
				ff.submat(k,0,k,1) = ff.submat(k,0,k,1) + ffslid;
				if((isSlip(k,0)==1)) {
					ffstick= -deformacion*(1-betaStick)*usstep*fn/std::max(delta,1e-7);;
					ff.submat(k,0, k ,1) = ff.submat(k, 0, k, 1) + ffstick;
				}
			}
			if (frictionModel==1){
				arma::mat vxy= arma::zeros(1,2);
				vxy(0,0)=v(0);
				vxy(0,1)=v(1);
				double normvxy= arma::norm(vxy, 2);
				fn=abs(fg*step(fondo - pos(k,2), -d/2, 0, 0, 1)); //se ha suavizado
				double coef;
				double normaSmooth;
				arma:: mat ffriccionVel = arma::zeros(1,2);

				if (normvxy < vstatic) {
					coef = step(normvxy, -vstatic, -us, vstatic, us);

				} else if ((vstatic<=normvxy) && (normvxy<vdynamic)) {
					coef = step(normvxy, vstatic, us, vdynamic, ud);

				} else {
					coef = ud;
				}
				if(normvxy<=1e-04) {
					normaSmooth = arma::as_scalar(a_1(0,0)*pow(normvxy,3)+a_1(1,0)*pow(normvxy,2)+a_1(2,0)*normvxy+a_1(3,0));
				} else{
					normaSmooth = normvxy;
				}
				ffriccionVel= -vxy* (coef /normaSmooth)*fn;
				ff.submat(k,0, k ,1) = ff.submat(k, 0, k, 1) + ffriccionVel;	
			}
		}
	}


	//F = 0.5 * dL * (MassMatrix_sp * ff) - (MSMatrix_sp * FF);
	F = 0.5 * dL * (MassMatrix * ff) - (MSMatrix * FF);
	// F.row(0) = F.row(0) + FF.row(0);
	// F.row(N-1) = F.row(N-1) - FF.row(N-1);

	if (F.has_nan()){
		std::cout << std::endl << "ERROR: NaN Detected on line with id = " << id << std::endl;
		std::cout << std::endl << "  F = " << F << std::endl;
		std::cout << std::endl << "  ff = " << ff << std::endl;
		std::cout << std::endl << "  FF = " << FF << std::endl;
		std::cout << std::endl << "  pos = " << pos << std::endl;
		std::cout << std::endl << "  vel = " << vel << std::endl;
		
		throw std::exception();
	}

	ten_1 = FF.row(0).t();
	ten_N = FF.row(N-1).t();
	F_1 = F.row(0).t();
	F_N = F.row(N-1).t();
	// F_1 = d*ff.row(0).t() + ten_1;
	// F_N = d*ff.row(N-1).t() - ten_N;


	pLineBcps[0]->forceBcp.rows(0,2) = pLineBcps[0]->forceBcp.rows(0,2) + F_1;
	pLineBcps[1]->forceBcp.rows(0,2) = pLineBcps[1]->forceBcp.rows(0,2) + F_N;
	pLineBcps[0]->temp = pLineBcps[0]->forceBcp; 
	pLineBcps[1]->temp = pLineBcps[1]->forceBcp;
}


void Line::print_out (void)
{

		std::cout << "Para la linea " << this->id << " , se ha leido:" << std::endl << std::endl;
		std::cout << "nNodos   " << this->nNodos << std::endl;
		std::cout << "p        " << this->p << std::endl;
		std::cout << "L        " << this->L << std::endl;
		std::cout << "rho0     " << this->rho0 << std::endl;
		std::cout << "d        " << this->d << std::endl;
		std::cout << "EA       " << this->EA << std::endl;
		std::cout << "beta     " << this->beta << std::endl;
		std::cout << "CB       " << this->CB << std::endl;
		std::cout << "Cmn      " << this->Cmn << std::endl;
		std::cout << "Cdn      " << this->Cdn << std::endl;
		std::cout << "Cdt      " << this->Cdt << std::endl;
		std::cout << "GK       " << this->GK << std::endl;
		std::cout << "GC       " << this->GC << std::endl;
		std::cout << "Gmu      " << this->Gmu << std::endl;
		std::cout << "Gvc      " << this->Gvc << std::endl;
		std::cout << "Dz       " << this->Dz << std::endl;
		std::cout << "pos_N  " << this->pos_N(0,0) << " " << this->pos_N(1,0) << " " << this->pos_N(2,0) << std::endl;
		std::cout << "pos_1  " << this->pos_1(0,0) << " " << this->pos_1(1,0) << " " << this->pos_1(2,0) << std::endl;
		std::cout << "smoothstep      " << this->smoothstep << std::endl;
		std::cout << "frictionModel      " << this->frictionModel << std::endl;
		std::cout << "vth      " << this->vth << std::endl;
		std::cout << "us      " << this->us << std::endl;
		std::cout << "ud      " << this->ud << std::endl;
		std::cout << "deltamax      " << this->deltamax << std::endl;
		std::cout << "vstatic      " << this->vstatic << std::endl;
		std::cout << "vdynamic      " << this->vdynamic << std::endl << std::endl;
}

void Line::initLine (void) 
{
	// Load BCP data
	pos_1 = pLineBcps[0]->posG_BCP;
	pos_N = pLineBcps[1]->posG_BCP;

	// Init Line
	if (pos_N(2,0)<fondo || pos_1(2,0)<fondo) throw 0;

	if (lineType == 1){
		int flag = 0;
		int flagTense = 0;
		xF=sqrt(pow((pos_N(0,0)-pos_1(0,0)),2)+pow((pos_N(1,0)-pos_1(1,0)),2));
		zF=(pos_N(2,0)-pos_1(2,0));
		if(xF==0){
			cosa=0; sina=0;		
		}else{
			cosa=(pos_N(0,0)-pos_1(0,0))/xF;
			sina=(pos_N(1,0)-pos_1(1,0))/xF;
		}
		double Li= sqrt(pow(xF,2)+pow(zF,2));
		if ( Li<L) {
			if ( pos_1(2,0)==fondo && pos_N(2,0)==fondo ) throw 2;
			if (xF<1e-5) throw 3;
			floor_flag = -1;
			this->qs_GetTen();
			this->qs_Solution();
			if(isnan(xc(0,0))) flagTense = 1;
			for(int ii=0;ii<N;ii=ii+1) {
				pos(3*ii,0) = pos_1(0,0)+cosa*xc(ii,0);
				pos(3*ii+1,0) = pos_1(1,0)+sina*xc(ii,0);
				pos(3*ii+2,0) = pos_1(2,0)+zc(ii,0);
				if (pos(3*ii+2,0)<=fondo && ii>=1) flag = 1; 
			}
			if ((flag == 1)&&pos_1(2,0)>fondo) throw 6;
			if ((flag == 1)&&((pos_1(2,0)==fondo)||(pos_N(2,0)==fondo))) {
				floor_flag = 1;
				this->qs_GetTen();
				this->qs_Solution();
				if(isnan(xc(0,0))) throw 5;
				for(int ii=0;ii<N;ii=ii+1) {
					pos(3*ii,0) = pos_1(0,0)+cosa*xc(ii,0);
					pos(3*ii+1,0) = pos_1(1,0)+sina*xc(ii,0);
					pos(3*ii+2,0) = pos_1(2,0)+zc(ii,0);
				}
			}
			arma::mat vec1(pos.rows(3,5)-pos.rows(0,2));
			double nvec1 = norm(vec1);
			ten_1 = (Te(0)/nvec1) * vec1;
			arma::mat vecN(pos.rows(3*N-6,3*N-4) - pos.rows(3*N-3,3*N-1));
			double nvecN = norm(vecN);
			ten_N = (Te(N-1)/nvecN) * vecN;
		}
		if ((Li>=L)||(flagTense == 1)) {
			double cost=xF/Li;
			double sint=zF/Li;
			for(int ii=0;ii<N;ii=ii+1) {
				xc(ii,0)=cost*(Li/L)*s(ii);
				zc(ii,0)=sint*(Li/L)*s(ii);
			}
			for(int ii=0;ii<N;ii=ii+1) {
				pos(3*ii,0) = pos_1(0,0)+cosa*xc(ii,0);
				pos(3*ii+1,0) = pos_1(1,0)+sina*xc(ii,0);
				pos(3*ii+2,0) = pos_1(2,0)+zc(ii,0);
			}
			arma::mat vec1(pos.rows(3,5)-pos.rows(0,2));
			double nvec1 = norm(vec1);
			ten_1 = (EA*(nvec1 - 1.0)/(0.5*(roots(1)+1)*dL*nvec1)) * vec1;
			arma::mat vecN(pos.rows(3*N-6,3*N-4) - pos.rows(3*N-3,3*N-1));
			double nvecN = norm(vecN);
			ten_N = (EA*(nvecN - 1.0)/(0.5*(roots(1)+1)*dL*nvecN)) * vecN;
			if ( Li==L && pos_1(2,0)==fondo && pos_N(2,0)==fondo ) {
				std::cout << "     WARNING: Mooring line " << id << " is laying on the floor " << std::endl;
			}else if (xF<1e-5){
				std::cout << "     WARNING: Mooring line " << id << " is vertical. " << std::endl;
			}else{
				std::cout << "     WARNING: Mooring line " << id << " tension is high. " << std::endl;
			}
		}
		floor_flag = 1;
		std::cout << "    Tension at anchor for line: " << id << " ; is: (" << ten_1(0) << " , " << ten_1(1) << " , " << ten_1(2) << " ) N" << std::endl;
		std::cout << "    Tension at fairlead for line: " << id << " ; is: (" << ten_N(0) << " , " << ten_N(1) << " , " << ten_N(2) << " ) N" << std::endl << std::endl;
	} else if (lineType == 2){
		floor_flag = -1;
		int flagTense = 0;
		xF=sqrt(pow((pos_N(0,0)-pos_1(0,0)),2)+pow((pos_N(1,0)-pos_1(1,0)),2));
		zF=(pos_N(2,0)-pos_1(2,0));
		if(xF==0){ // Vertical case
			cosa=0; sina=0;		
		}else{ // other
			cosa=(pos_N(0,0)-pos_1(0,0))/xF;
			sina=(pos_N(1,0)-pos_1(1,0))/xF;
		}
		double Li= sqrt(pow(xF,2)+pow(zF,2));
		if (Li<L) {
			this->qs_GetTen();
			this->qs_Solution();
			if(isnan(xc(0,0))) flagTense = 1;
			for(int ii=0;ii<N;ii=ii+1) {
				pos(3*ii,0) = pos_1(0,0)+cosa*xc(ii,0);
				pos(3*ii+1,0) = pos_1(1,0)+sina*xc(ii,0);
				pos(3*ii+2,0) = pos_1(2,0)+zc(ii,0);
				if (pos(3*ii+2,0)<=fondo) throw 1;			
			}
			arma::mat vec1(pos.rows(3,5)-pos.rows(0,2));
			double nvec1 = norm(vec1);
			ten_1 = (Te(0)/nvec1) * vec1;
			arma::mat vecN(pos.rows(3*N-6,3*N-4) - pos.rows(3*N-3,3*N-1));
			double nvecN = norm(vecN);
			ten_N = (Te(N-1)/nvecN) * vecN;
		}
		if ((Li>=L)||(flagTense == 1)) {
			double cost=xF/Li;
			double sint=zF/Li;
			for(int ii=0;ii<N;ii=ii+1) {
				xc(ii,0)=cost*(Li/L)*s(ii);
				zc(ii,0)=sint*(Li/L)*s(ii);
			}
			for(int ii=0;ii<N;ii=ii+1) {
				pos(3*ii,0) = pos_1(0,0)+cosa*xc(ii,0);
				pos(3*ii+1,0) = pos_1(1,0)+sina*xc(ii,0);
				pos(3*ii+2,0) = pos_1(2,0)+zc(ii,0);
			}
			arma::mat vec1(pos.rows(3,5)-pos.rows(0,2));
			double nvec1 = norm(vec1);
			ten_1 = (EA*(nvec1 - 1.0)/(0.5*(roots(1)+1)*dL*nvec1)) * vec1;
			arma::mat vecN(pos.rows(3*N-6,3*N-4) - pos.rows(3*N-3,3*N-1));
			double nvecN = norm(vecN);
			ten_N = (EA*(nvecN - 1.0)/(0.5*(roots(1)+1)*dL*nvecN)) * vecN;
			std::cout << "    WARNING: Towing line " << id << " tension is high. " << std::endl;
		}
		std::cout << "    Tension at anchor for line: " << id << " ; is: (" << ten_1(0) << " , " << ten_1(1) << " , " << ten_1(2) << " ) N" << std::endl;
		std::cout << "    Tension at fairlead for line: " << id << " ; is: (" << ten_N(0) << " , " << ten_N(1) << " , " << ten_N(2) << " ) N" << std::endl << std::endl;
	} else if (lineType == 3){
		xF=sqrt(pow((pos_N(0,0)-pos_1(0,0)),2)+pow((pos_N(1,0)-pos_1(1,0)),2));
		zF=(pos_N(2,0)-pos_1(2,0));
		double Li= sqrt(pow(xF,2)+pow(zF,2));
		if (Li<L) throw 4;
		if(xF==0){
			cosa=0; sina=0;		
		}else{
			cosa=(pos_N(0,0)-pos_1(0,0))/xF;
			sina=(pos_N(1,0)-pos_1(1,0))/xF;
		}
		double cost=xF/Li;
		double sint=zF/Li;
		for(int ii=0;ii<N;ii=ii+1) {
			xc(ii,0)=cost*(Li/L)*s(ii);
			zc(ii,0)=sint*(Li/L)*s(ii);
			pos(3*ii,0) = pos_1(0,0)+cosa*xc(ii,0);
			pos(3*ii+1,0) = pos_1(1,0)+sina*xc(ii,0);
			pos(3*ii+2,0) = pos_1(2,0)+zc(ii,0);
		}
		arma::mat vec1(pos.rows(3,5)-pos.rows(0,2));
		double nvec1 = norm(vec1);
		ten_1 = (EA*(nvec1 - 1.0)/(0.5*(this->roots(1)+1)*dL*nvec1)) * vec1;
		arma::mat vecN(pos.rows(3*N-6,3*N-4) - pos.rows(3*N-3,3*N-1));
		double nvecN = norm(vecN);
		ten_N = (EA*(nvecN - 1.0)/(0.5*(this->roots(1)+1)*dL*nvecN)) * vecN;
		std::cout << "     Tension at anchor for line: " << id << " ; is: (" << ten_1(0) << " , " << ten_1(1) << " , " << ten_1(2) << " ) N" << std::endl;
		std::cout << "     Tension at fairlead for line: " << id << " ; is: (" << ten_N(0) << " , " << ten_N(1) << " , " << ten_N(2) << " ) N" << std::endl << std::endl;
		floor_flag = -1;
	} else {
		throw std::invalid_argument( "Type of line not available." );
	}
}

void Line::OpenOutputFilesASCII (std::string path)
{
	char buffer1[50], buffer2[50], buffer3[50], buffer4[50], buffer5[50];

	int nn1 = sprintf(buffer1,"NodePosX_%d.txt", GetId());
	int nn2 = sprintf(buffer2,"NodePosY_%d.txt", GetId());
	int nn3 = sprintf(buffer3,"NodePosZ_%d.txt", GetId());
	int nn4 = sprintf(buffer4,"EndsTen_%d.txt", GetId());
	int nn5 = sprintf(buffer5,"LineTen_%d.txt", GetId());

	std::string file_path1 = JoinPath(path, buffer1);
	std::string file_path2 = JoinPath(path, buffer2);
	std::string file_path3 = JoinPath(path, buffer3);
	std::string file_path4 = JoinPath(path, buffer4);
	std::string file_path5 = JoinPath(path, buffer5);

	pfile_xpos = fopen (file_path1.c_str(),"w");
	pfile_ypos = fopen (file_path2.c_str(),"w");
	pfile_zpos = fopen (file_path3.c_str(),"w");
	pfile_ten =  fopen (file_path4.c_str(),"w");
	pfile_ten_line =  fopen (file_path5.c_str(),"w");
}

void Line::CloseOutputFilesASCII (void)
{
	fclose(pfile_xpos);
	fclose(pfile_ypos);
	fclose(pfile_zpos);
	fclose(pfile_ten);
	fclose(pfile_ten_line);
}

void Line::WriteOut (double t) 
{
	int ii;

	fprintf(pfile_xpos, "%f    ", t);
	for(ii=0;ii<this->N;ii=ii+1) fprintf(pfile_xpos, "%f    ", this->pos(ii,0));
	fprintf(pfile_xpos, "\n");

	fprintf(pfile_ypos, "%f    ", t);
	for(ii=0;ii<this->N;ii=ii+1) fprintf(pfile_ypos, "%f    ", this->pos(ii,1));
	fprintf(pfile_ypos, "\n");

	fprintf(pfile_zpos, "%f    ", t);
	for(ii=0;ii<this->N;ii=ii+1) fprintf(pfile_zpos, "%f    ", this->pos(ii,2));
	fprintf(pfile_zpos, "\n");

	fprintf(pfile_ten, "%f    %f    %f    %f    %f    %f    %f \n",t,ten_1(0,0),ten_1(1,0),ten_1(2,0),ten_N(0,0),ten_N(1,0),ten_N(2,0));
	
	if (t==0.0){
		fprintf(pfile_ten_line, "%f    ", t);
		for(ii=0;ii<this->N;ii=ii+1) fprintf(pfile_ten_line, "%f    ", this->s(ii,0));
		fprintf(pfile_ten_line, "\n");
	}
	fprintf(pfile_ten_line, "%f    ", t);
	for(ii=0;ii<this->N;ii=ii+1) fprintf(pfile_ten_line, "%f    ", this->T(ii,0));
	fprintf(pfile_ten_line, "\n");
}
