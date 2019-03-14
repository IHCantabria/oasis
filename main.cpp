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

struct solver_data {
	int nMoorLines;
	int nTowLines;
	int nTenLines;
	int nSistema;
	int nSistema2;
	MooringLine * MoorLine;
	TowingLine * TowLine;
	TensorLine * TenLine;
	int nAnchBCPs;
	int nFairBCPs;
	AnchorBCP * AnchBCP;
	FairleadBCP * FairBCP;
} Lines;

arma::mat fun(double t, arma::mat y, solver_data Lines){

	arma::mat yprime = arma::zeros(size(y));
	int i0;

	// Copy info from y to the objects.
	int ini = 0;
	for(int ii=0;ii<Lines.nMoorLines;ii=ii+1){
		for(int jj=0;jj<Lines.MoorLine[ii].N;jj=jj+1){
			Lines.MoorLine[ii].pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nTowLines;ii=ii+1){
		for(int jj=0;jj<Lines.TowLine[ii].N;jj=jj+1){
			Lines.TowLine[ii].pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nTenLines;ii=ii+1){
		for(int jj=0;jj<Lines.TenLine[ii].N;jj=jj+1){
			Lines.TenLine[ii].pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nMoorLines;ii=ii+1){
		for(int jj=0;jj<Lines.MoorLine[ii].N;jj=jj+1){
			Lines.MoorLine[ii].vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nTowLines;ii=ii+1){
		for(int jj=0;jj<Lines.TowLine[ii].N;jj=jj+1){
			Lines.TowLine[ii].vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nTenLines;ii=ii+1){
		for(int jj=0;jj<Lines.TenLine[ii].N;jj=jj+1){
			Lines.TenLine[ii].vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}

	// Set boundary conditions on pos and vel
	/*
	for(int ii=0;ii<Lines.nAnchBCPs;ii=ii+1){
		Lines.AnchBCP[ii].getValues();
		for(int jj=0;jj<Lines.AnchBCP[ii].nLinesBCP;jj=jj+1){
			if(Lines.AnchBCP[ii].BCPLineNode[jj] == 1) i0 = 0;
			if(Lines.AnchBCP[ii].BCPLineNode[jj] == 2) i0 = 0;
		}
	}

	for(int ii=0;ii<Lines.nFairBCPs;ii=ii+1){

	}
	*/

	// Compute forces vector for the different objects
	for(int ii=0;ii<Lines.nMoorLines;ii=ii+1){
		Lines.MoorLine[ii].SEM_computeF();
		Lines.MoorLine[ii].F.row(0) = arma::zeros(1,3);
		Lines.MoorLine[ii].F.row(Lines.MoorLine[ii].N-1) = arma::zeros(1,3);
		Lines.MoorLine[ii].acc = arma::solve(Lines.MoorLine[ii].dL * Lines.MoorLine[ii].MM,Lines.MoorLine[ii].F);
	}
	for(int ii=0;ii<Lines.nTowLines;ii=ii+1){
		Lines.TowLine[ii].SEM_computeF();
		Lines.TowLine[ii].F.row(0) = arma::zeros(1,3);
		Lines.TowLine[ii].F.row(Lines.TowLine[ii].N-1) = arma::zeros(1,3);
		Lines.TowLine[ii].acc = arma::solve(Lines.TowLine[ii].dL * Lines.TowLine[ii].MM,Lines.TowLine[ii].F);
	}
	for(int ii=0;ii<Lines.nTenLines;ii=ii+1){
		Lines.TenLine[ii].SEM_computeF();
		Lines.TenLine[ii].F.row(0) = arma::zeros(1,3);
		Lines.TenLine[ii].F.row(Lines.TenLine[ii].N-1) = arma::zeros(1,3);
		Lines.TenLine[ii].acc = arma::solve(Lines.TenLine[ii].dL * Lines.TenLine[ii].MM,Lines.TenLine[ii].F);
	}


	// Copy info from the objects to yprime
	yprime.rows(0,Lines.nSistema2-1) = y.rows(Lines.nSistema2,Lines.nSistema-1);
	ini = Lines.nSistema2;
	for(int ii=0;ii<Lines.nMoorLines;ii=ii+1){
		for(int jj=0;jj<Lines.MoorLine[ii].N;jj=jj+1){
			yprime.rows(ini,ini+2) = Lines.MoorLine[ii].acc.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nTowLines;ii=ii+1){
		for(int jj=0;jj<Lines.TowLine[ii].N;jj=jj+1){
			yprime.rows(ini,ini+2) = Lines.TowLine[ii].acc.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<Lines.nTenLines;ii=ii+1){
		for(int jj=0;jj<Lines.TenLine[ii].N;jj=jj+1){
			yprime.rows(ini,ini+2) = Lines.TenLine[ii].acc.row(jj).t();
			ini = ini + 3;
		}
	}


	return yprime;
}

void direction(arma::mat y, arma::mat& F, arma::mat& dy, arma::mat y0, double t, double h, solver_data Lines){

	double dx = 1e-8, inv_dx = 1e8;

	arma::mat J = arma::zeros(Lines.nSistema,Lines.nSistema);
	arma::mat yprime0 = fun(t+h, y, Lines);
	arma::mat yprime, e;

	for(int ii=0;ii<Lines.nSistema;ii=ii+1){
		e = arma::zeros(Lines.nSistema,1);
		e.row(ii) = dx;
		yprime = fun(t+h, y + e, Lines);
		J.col(ii) = inv_dx * (yprime - yprime0);
	}

	arma::mat I = arma::eye(Lines.nSistema,Lines.nSistema);

	arma::mat M = I - h * J;
	F = y - y0 - h*yprime0;

	dy = arma::solve( M, F);
}

int main () {

	int flag_read_eq, flag_write_eq;
	int nLines, nMoorLines, nTowLines, nTenLines;
	int nBCPs, nAnchBCPs, nFairBCPs;
	int nNodosTotal, nSistema;

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
	datosProblema >> flag_read_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_write_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema.close();


	//LEO DE FICHERO Y PINTO EN PANTALLA CUANTAS LINEAS SE VAN A ESTUDIAR
	std::ifstream datosMoorings ("datosMoorings.dat");
	datosMoorings >> nMoorLines; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> nTowLines; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> nTenLines; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings.close();

	nLines=nMoorLines+nTowLines+nTenLines;

	std::cout << std::endl << "El numero de lineas a estudiar es: " << nLines << std::endl << std::endl;

	Lines.nMoorLines = nMoorLines;
	Lines.nTenLines = nTenLines;
	Lines.nTowLines = nTowLines;

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	MooringLine * MoorLine = new MooringLine[nMoorLines];
	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nMoorLines; ii=ii+1){
		MoorLine[ii].set_nLine(ii+1);
		try
		{
		MoorLine[ii].leer_datosMoorings();
		//MoorLine[ii].print_out();
		nNodosTotal=nNodosTotal+MoorLine[ii].N;
		if(flag_read_eq==0) MoorLine[ii].initLine();
		MoorLine[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << MoorLine[ii].nLine << " is under the floor level." << std::endl << std::endl;
			if (e==2) std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl << std::endl;
			if (e==3) std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " is not tense and vertical. It should be pretensed. " << std::endl << std::endl;
			if (e==6) std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " touches the seafloor althoug none of its ends are there. " << std::endl << std::endl;
			if (e==5) {std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " initial shape can't be computed with QS method. " << std::endl;
			std::cout<< "Try reaching the desired initial condition with the dynamic method from a different initial condition " << std::endl << std::endl;}
			return 0;
		}

		//MoorLine[ii].SEM_computeA();	
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	TowingLine * TowLine = new TowingLine[nTowLines];
	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nTowLines; ii=ii+1){

		TowLine[ii].set_nLine(ii+1+nMoorLines);
		try
		{
		TowLine[ii].leer_datosMoorings();
		//TowLine[ii].print_out();
		nNodosTotal=nNodosTotal+TowLine[ii].N;
		if(flag_read_eq==0) TowLine[ii].initLine();
		TowLine[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << TowLine[ii].nLine << " is under the floor level. " << std::endl << std::endl;
			if (e==1) std::cout<< "ERROR: Towing line " << TowLine[ii].nLine << " touches the seafloor. " << std::endl << std::endl;
			if (e==5) {std::cout<< "ERROR: Towing line " << MoorLine[ii].nLine << " initial shape can't be computed with QS method. " << std::endl;
			std::cout<< "Try reaching the desired initial condition with the dynamic method from a different initial condition " << std::endl << std::endl;}
			return 0;
		}

		//TowLine[ii].SEM_computeA();		
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	TensorLine * TenLine = new TensorLine[nTenLines];
	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nTenLines; ii=ii+1){

		TenLine[ii].set_nLine(ii+1+nMoorLines+nTowLines);
		try
		{
		TenLine[ii].leer_datosMoorings();
		//TenLine[ii].print_out();
		nNodosTotal=nNodosTotal+TenLine[ii].N;
		if(flag_read_eq==0) TenLine[ii].initLine();
		TenLine[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << TenLine[ii].nLine << " is under the floor level. " << std::endl << std::endl;
			if (e==4) std::cout<< "ERROR: Tensor line " << TenLine[ii].nLine << " is not tense. " << std::endl << std::endl;
			return 0;
		}

		//TenLine[ii].SEM_computeA();	
	}

	//LEO DE FICHERO CUANTOS BCPs SE VAN A ESTUDIAR
	std::ifstream datosBCPs ("datosBCPs.dat");
	datosBCPs >> nFairBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nAnchBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs.close();

	nBCPs=nAnchBCPs+nFairBCPs;

	Lines.nAnchBCPs = nAnchBCPs;
	Lines.nFairBCPs = nFairBCPs;

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	FairleadBCP * FairBCP = new FairleadBCP[nFairBCPs];
	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nFairBCPs; ii=ii+1){
		FairBCP[ii].set_nBCP(ii+1);
		FairBCP[ii].leer_datosBCPs();
		FairBCP[ii].getValues(0.0);
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	AnchorBCP * AnchBCP = new AnchorBCP[nAnchBCPs];
	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nAnchBCPs; ii=ii+1){
		AnchBCP[ii].set_nBCP(ii+1+nFairBCPs);
		AnchBCP[ii].leer_datosBCPs();
		if(AnchBCP[ii].BCPLineType[0] == 1){
			if(AnchBCP[ii].BCPLineNode[0] == 1){
				AnchBCP[ii].pos = MoorLine[AnchBCP[ii].BCPLineIndex[0]].posAnch;
			}
			if(AnchBCP[ii].BCPLineNode[0] == 2){
				AnchBCP[ii].pos = MoorLine[AnchBCP[ii].BCPLineIndex[0]].posFair;
			}
		}
		if(AnchBCP[ii].BCPLineType[0] == 2){
			if(AnchBCP[ii].BCPLineNode[0] == 1){
				AnchBCP[ii].pos = TowLine[AnchBCP[ii].BCPLineIndex[0]].posAnch;
			}
			if(AnchBCP[ii].BCPLineNode[0] == 2){
				AnchBCP[ii].pos = TowLine[AnchBCP[ii].BCPLineIndex[0]].posFair;
			}
		}
		if(AnchBCP[ii].BCPLineType[0] == 3){
			if(AnchBCP[ii].BCPLineNode[0] == 1){
				AnchBCP[ii].pos = TenLine[AnchBCP[ii].BCPLineIndex[0]].posAnch;
			}
			if(AnchBCP[ii].BCPLineNode[0] == 2){
				AnchBCP[ii].pos = TenLine[AnchBCP[ii].BCPLineIndex[0]].posFair;
			}
		}
	}


	// Inicio el vector del sistema
	nSistema=3*2*nNodosTotal;
	Lines.nSistema = 3*2*nNodosTotal;
	Lines.nSistema2 = 3*nNodosTotal;

	arma::mat y = arma::zeros(nSistema,1);
	arma::mat yprime = arma::zeros(nSistema,1);

	int ini=0;
	if(flag_read_eq==0){
		for(int ii=0; ii<nMoorLines; ii=ii+1){
			for(int jj=0; jj<MoorLine[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = MoorLine[ii].pos.row(jj).t();
					ini=ini+3;
			}
		}
		for(int ii=0; ii<nTowLines; ii=ii+1){
			for(int jj=0; jj<TowLine[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = TowLine[ii].pos.row(jj).t();
					ini=ini+3;
			}		
		}
		for(int ii=0; ii<nTenLines; ii=ii+1){
			for(int jj=0; jj<TenLine[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = TenLine[ii].pos.row(jj).t();
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
		for(int ii=0; ii<nMoorLines; ii=ii+1){
			for(int jj=0; jj<MoorLine[ii].N; jj=jj+1){
				MoorLine[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}
		}
		for(int ii=0; ii<nTowLines; ii=ii+1){
			for(int jj=0; jj<TowLine[ii].N; jj=jj+1){
				TowLine[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}		
		}
		for(int ii=0; ii<nTenLines; ii=ii+1){
			for(int jj=0; jj<TenLine[ii].N; jj=jj+1){
				TenLine[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}		
		}
	}

	for(int ii=0; ii<nMoorLines; ii=ii+1) MoorLine[ii].write_out(t);
	for(int ii=0; ii<nTowLines; ii=ii+1) TowLine[ii].write_out(t);
	for(int ii=0; ii<nTenLines; ii=ii+1) TenLine[ii].write_out(t);

	Lines.MoorLine = MoorLine;
	Lines.TowLine = TowLine;
	Lines.TenLine = TenLine;

	//Lines.AnchBCP = AnchBCP;
	//Lines.FairBCP = FairBCP;

	int solver_flag = 2;

	if (solver_flag == 1){
		double dtrk = 1e-8;
		double t_old = 0.0;
		arma::mat K1, K2, K3;
		std::cout<< "The temporal integration begins." << std::endl << std::endl;
		std::cout<< "    t = " << t << std::endl;
		do{
			K1 = fun(t,y,Lines);
			K2 = fun(t+dtrk,y+dtrk*K1,Lines);
			K3 = fun(t+0.5*dtrk,y+0.25*dtrk*(K1+K2),Lines);
			y = y + dtrk* ( (1.0/6.0)*(K1 + K2) + (2.0/3.0)*K3);
			t = t + dtrk;
			//std::cout << "                       mean_acc = " << arma::norm(Lines.MoorLine[0].acc,2)/Lines.MoorLine[0].N << std::endl;
			//std::cout << "                              T = " << arma::norm(Lines.MoorLine[0].tenFair,2) << std::endl;
			if (t >= t_old + dt){
				t_old = t;
				std::cout<< "    t = " << t << std::endl;
				std::cout << "                       mean_acc = " << arma::norm(Lines.MoorLine[0].acc,2)/Lines.MoorLine[0].N << std::endl;
				for(int ii=0; ii<nMoorLines; ii=ii+1) Lines.MoorLine[ii].write_out(t);
				for(int ii=0; ii<nTowLines; ii=ii+1) Lines.TowLine[ii].write_out(t);
				for(int ii=0; ii<nTenLines; ii=ii+1) Lines.TenLine[ii].write_out(t);
			}
		} while (t<=t_max);
	}
	if (solver_flag == 2){
		double dtei = 1e-2;
		double t_old = 0.0;
		arma::mat y0;
		arma::mat F = arma::zeros(Lines.nSistema,1);
		arma::mat dy = arma::zeros(Lines.nSistema,1);
		std::cout<< "The temporal integration begins." << std::endl << std::endl;
		std::cout<< "    t = " << t << std::endl;
		do{
			y0 = y;
			//std::cout << "                              T = " << arma::norm(Lines.MoorLine[0].tenAnch,2) << std::endl;
			//std::cout << "                       mean_acc = " << arma::norm(Lines.MoorLine[0].acc,2)/Lines.MoorLine[0].N << std::endl;
			direction(y, F, dy, y0, t, dtei, Lines);
			do{
				y = y - dy;
				direction(y, F, dy, y0, t, dtei, Lines);
			} while (arma::norm(F,2)/nSistema > 1e-6);
			t = t + dtei;	
			if (t >= t_old + dt){
				t_old = t;
				std::cout<< "    t = " << t << std::endl;
				for(int ii=0; ii<nMoorLines; ii=ii+1) Lines.MoorLine[ii].write_out(t);
				for(int ii=0; ii<nTowLines; ii=ii+1) Lines.TowLine[ii].write_out(t);
				for(int ii=0; ii<nTenLines; ii=ii+1) Lines.TenLine[ii].write_out(t);
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