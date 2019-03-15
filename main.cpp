/*
MAIN DE LA IMPLEMENTACION DE NuevosFEM EN C++
*/

//LIBRERIAS Y OTROS COMANDOS
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include "classes.h"
#include <armadillo>

double PI;
double g;
double rhoW;
double fondo;


arma::mat fun(double t, arma::mat y, solver_data SD){

	arma::mat yprime = arma::zeros(size(y));
	int i0;

	// Copy info from y to the objects.
	int ini = 0;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii].N;jj=jj+1){
			SD.Lines[ii].pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii].N;jj=jj+1){
			SD.Lines[ii].vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}

	// Set boundary conditions on pos and vel
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		SD.Lines[ii].LineBCP[0]->getValues(t);
		SD.Lines[ii].LineBCP[1]->getValues(t);
		SD.Lines[ii].pos.row(0)                = SD.Lines[ii].LineBCP[0]->pos.t();
		SD.Lines[ii].pos.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->pos.t();
		SD.Lines[ii].vel.row(0)                = SD.Lines[ii].LineBCP[0]->vel.t();
		SD.Lines[ii].vel.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->vel.t();
	}



	// Compute forces vector for the different objects
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		SD.Lines[ii].SEM_computeF();
		SD.Lines[ii].F.row(0)                = SD.Lines[ii].LineBCP[0]->acc.t();
		SD.Lines[ii].F.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[0]->acc.t();
		SD.Lines[ii].acc = arma::solve(SD.Lines[ii].dL * SD.Lines[ii].MM, SD.Lines[ii].F);
	}



	// Copy info from the objects to yprime
	yprime.rows(0,SD.nSistema2-1) = y.rows(SD.nSistema2,SD.nSistema-1);
	ini = SD.nSistema2;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii].N;jj=jj+1){
			yprime.rows(ini,ini+2) = SD.Lines[ii].acc.row(jj).t();
			ini = ini + 3;
		}
	}


	return yprime;
}

arma::mat jac(double t, arma::mat y, solver_data SD){
	double dx = 1e-8, inv_dx = 1e8;
	arma::mat J = arma::zeros(SD.nSistema,SD.nSistema);
	arma::mat yprime0 = fun(t, y, SD);
	arma::mat yprime, e;
	for(int ii=0;ii<SD.nSistema;ii=ii+1){
		e = arma::zeros(SD.nSistema,1);
		e.row(ii) = dx;
		yprime = fun(t, y + e, SD);
		J.col(ii) = inv_dx * (yprime - yprime0);
	}
	return J;
}

void direction(arma::mat y, arma::mat& F, arma::mat& dy, arma::mat y0, double t, double h, solver_data SD){

	double dx = 1e-8, inv_dx = 1e8;

	arma::mat J = arma::zeros(SD.nSistema,SD.nSistema);
	arma::mat yprime0 = fun(t+h, y, SD);
	arma::mat yprime, e;

	for(int ii=0;ii<SD.nSistema;ii=ii+1){
		e = arma::zeros(SD.nSistema,1);
		e.row(ii) = dx;
		yprime = fun(t+h, y + e, SD);
		J.col(ii) = inv_dx * (yprime - yprime0);
	}

	arma::mat I = arma::eye(SD.nSistema,SD.nSistema);

	arma::mat M = I - h * J;
	F = y - y0 - h*yprime0;

	dy = arma::solve( M, F);
}

int main () {

	int flag_read_eq, flag_write_eq;
	int nLines;
	int nBCPs, nAnchBCPs, nFairBCPs;
	int nNodosTotal, nSistema;
	int solver_flag;
	solver_data SD;

	double t;
	double t_max;
	double dt;

	PI=acos(-1.0);
	t = 0.0;

	nNodosTotal=0;

 	std::ifstream datosProblema ("datosProblema.dat");
	datosProblema >> g;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> rhoW;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> fondo;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> dt;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> t_max;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> solver_flag;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_read_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_write_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema.close();


	//LEO DE FICHERO CUANTOS BCPs SE VAN A ESTUDIAR
	std::ifstream datosBCPs ("datosBCPs.dat");
	datosBCPs >> nFairBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nAnchBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs.close();

	nBCPs = nFairBCPs + nAnchBCPs;

	BCP * BCPs [nBCPs];

	BCP * F_BCPs = new FairleadBCP [nFairBCPs];
	for(int ii=0; ii<nFairBCPs; ii=ii+1){
		F_BCPs[ii].set_nBCP(ii+1);
		F_BCPs[ii].leer_datosBCPs();
		F_BCPs[ii].getValues(0.0);
		BCPs[ii] = &F_BCPs[ii];
	}

	BCP * A_BCPs = new AnchorBCP [nAnchBCPs];
	for(int ii=0; ii<nAnchBCPs; ii=ii+1){
		A_BCPs[ii].set_nBCP(ii+1+nFairBCPs);
		A_BCPs[ii].leer_datosBCPs();
		BCPs[ii+nFairBCPs] = &A_BCPs[ii];
	}

	std::cout << "BCPs[0].getValues(3.5)" << std::endl;
	BCPs[0]->getValues(3.5);
	std::cout << "BCPs[1].getValues(3.5)" << std::endl;
	BCPs[1]->getValues(3.5);


	//LEO DE FICHERO Y PINTO EN PANTALLA CUANTAS LINEAS SE VAN A ESTUDIAR
	std::ifstream datosLines ("datosLines.dat");
	datosLines >> nLines; datosLines.ignore(std::numeric_limits<int>::max(), '\n');
	datosLines.close();

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	Line * Lines = new Line[nLines];
	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nLines; ii=ii+1){
		Lines[ii].set_nLine(ii+1);
		try
		{
		Lines[ii].leer_datosLines();
		Lines[ii].LineBCP[0] = BCPs[Lines[ii].BCP_1-1];
		Lines[ii].LineBCP[1] = BCPs[Lines[ii].BCP_N-1];
		Lines[ii].pos_1 = BCPs[Lines[ii].BCP_1-1]->pos;
		Lines[ii].pos_N = BCPs[Lines[ii].BCP_N-1]->pos;
		nNodosTotal=nNodosTotal+Lines[ii].N;
		if(flag_read_eq==0) Lines[ii].initLine();
		Lines[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << Lines[ii].nLine << " is under the floor level." << std::endl << std::endl;
			if (e==1) std::cout<< "ERROR: Line " << Lines[ii].nLine << " touches the seafloor and it shouldn't. " << std::endl << std::endl;
			if (e==2) std::cout<< "ERROR: Line " << Lines[ii].nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl << std::endl;
			if (e==3) std::cout<< "ERROR: Line " << Lines[ii].nLine << " is not tense and vertical. It should be pretensed. " << std::endl << std::endl;
			if (e==4) std::cout<< "ERROR: Line " << Lines[ii].nLine << " is not tense and it should. " << std::endl << std::endl;
			if (e==5) std::cout<< "ERROR: Line " << Lines[ii].nLine << " initial shape can't be computed with QS method. " << std::endl;
			if (e==6) std::cout<< "ERROR: Line " << Lines[ii].nLine << " touches the seafloor althoug none of its ends are there. " << std::endl << std::endl;
			return 0;
		}
	}

	// Inicio el vector del sistema
	nSistema=3*2*nNodosTotal;
	SD.nSistema = 3*2*nNodosTotal;
	SD.nSistema2 = 3*nNodosTotal;

	arma::mat y = arma::zeros(nSistema,1);
	arma::mat yprime = arma::zeros(nSistema,1);

	int ini=0;
	if(flag_read_eq==0){
		for(int ii=0; ii<nLines; ii=ii+1){
			for(int jj=0; jj<Lines[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = Lines[ii].pos.row(jj).t();
					ini=ini+3;
			}
		}

	}else if (flag_read_eq==1){
	 	std::ifstream equi("Equilibrio.dat");
			for(int ii=0;ii<nSistema;ii=ii+1) {
				equi >> y(ii,0);    equi.ignore(std::numeric_limits<int>::max(), '\n');
			}
		equi.close();

		ini=0;
		for(int ii=0; ii<nLines; ii=ii+1){
			for(int jj=0; jj<Lines[ii].N; jj=jj+1){
				Lines[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}
		}

	}

	for(int ii=0; ii<nLines; ii=ii+1) Lines[ii].write_out(t);

	SD.nLines = nLines;
	SD.Lines = Lines;

	//BDF S (t,y,*fun,*jac,SD);

	if (solver_flag == 1){
		double dtrk = 1e-8;
		double t_old = 0.0;
		arma::mat K1, K2, K3;
		std::cout<< "The temporal integration begins." << std::endl << std::endl;
		std::cout<< "    t = " << t << std::endl;
		do{
			K1 = fun(t,y,SD);
			K2 = fun(t+dtrk,y+dtrk*K1,SD);
			K3 = fun(t+0.5*dtrk,y+0.25*dtrk*(K1+K2),SD);
			y = y + dtrk* ( (1.0/6.0)*(K1 + K2) + (2.0/3.0)*K3);
			t = t + dtrk;
			if (t >= t_old + dt){
				t_old = t;
				std::cout<< "    t = " << t << std::endl;
				for(int ii=0; ii<nLines; ii=ii+1) Lines[ii].write_out(t);
			}
		} while (t<=t_max);
	}
	if (solver_flag == 2){
		double dtei = 1e-4;
		double atol = 1e-6;
		double t_old = 0.0;
		arma::mat y0;
		arma::mat F = arma::zeros(nSistema,1);
		arma::mat dy = arma::zeros(nSistema,1);
		std::cout<< "The temporal integration begins." << std::endl << std::endl;
		std::cout<< "    t = " << t << std::endl;
		do{
			y0 = y;
			direction(y, F, dy, y0, t, dtei, SD);
			do{
				y = y - dy;
				direction(y, F, dy, y0, t, dtei, SD);
			} while (arma::norm(F,2)/nSistema > atol);
			t = t + dtei;
			if (t>0.02) dtei = 1e-3;
			if (t>0.04) dtei = 1e-2;
			if (t>2.9) dtei = 1e-4;
			std::cout<< "    t = " << t << std::endl;
			if (t >= t_old + dt){
				t_old = t;
				//std::cout<< "    t = " << t << std::endl;
				for(int ii=0; ii<nLines; ii=ii+1) Lines[ii].write_out(t);
			}
		} while (t<=t_max);
	}


	if (flag_write_eq == 1) {
		std::ofstream equi("Equilibrio.dat");
			for(int ii=0;ii<nSistema;ii=ii+1) equi << y(ii,0) << std::endl;
		equi.close();
	}

	return 0;
}