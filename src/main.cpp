/*
MAIN DE LA IMPLEMENTACION DE NuevosFEM EN C++
*/

//LIBRERIAS Y OTROS COMANDOS
#include <iostream>
#include <fstream>
#include <limits>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include <armadillo>
#include <ctime>
#include "ODE_solvers/ODE_solvers.hpp"


double PI;
double g;
double rhoW;
double fondo;
int nCalls = 0;


arma::mat fun(double t, arma::mat y, solver_data SD){

	nCalls = nCalls + 1;

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

	// Set boundary conditions on pos and vel if the BCP is not a joint
	if (SD.Lines[0].LineBCP[0]->tBCP != t){
		for(int ii=0;ii<SD.nLines;ii=ii+1){		
			if(SD.Lines[ii].LineBCP[0]->typeBCP != 3){
				SD.Lines[ii].LineBCP[0]->getValues(t);
			}
			if(SD.Lines[ii].LineBCP[1]->typeBCP != 3){
				SD.Lines[ii].LineBCP[1]->getValues(t);
			}
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii].LineBCP[0]->typeBCP != 3){
			SD.Lines[ii].pos.row(0)                = SD.Lines[ii].LineBCP[0]->pos.t();
			SD.Lines[ii].vel.row(0)                = SD.Lines[ii].LineBCP[0]->vel.t();
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP != 3){
			SD.Lines[ii].pos.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->pos.t();
			SD.Lines[ii].vel.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->vel.t();
		}
	}

	// Set boundary conditions on pos and vel if the BCP is a joint
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii].LineBCP[0]->typeBCP == 3){
			SD.Lines[ii].LineBCP[0]->posLines.row(SD.Lines[ii].LineBCP[0]->iLJ) = SD.Lines[ii].pos.row(0);
			SD.Lines[ii].LineBCP[0]->velLines.row(SD.Lines[ii].LineBCP[0]->iLJ) = SD.Lines[ii].vel.row(0);
			SD.Lines[ii].LineBCP[0]->iLJ = SD.Lines[ii].LineBCP[0]->iLJ + 1;
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP == 3){
			SD.Lines[ii].LineBCP[1]->posLines.row(SD.Lines[ii].LineBCP[1]->iLJ) = SD.Lines[ii].pos.row(SD.Lines[ii].N-1);
			SD.Lines[ii].LineBCP[1]->velLines.row(SD.Lines[ii].LineBCP[1]->iLJ) = SD.Lines[ii].vel.row(SD.Lines[ii].N-1);
			SD.Lines[ii].LineBCP[1]->iLJ = SD.Lines[ii].LineBCP[1]->iLJ + 1;
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){		
		if(SD.Lines[ii].LineBCP[0]->typeBCP == 3){
			SD.Lines[ii].LineBCP[0]->getValues(t);
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP == 3){
			SD.Lines[ii].LineBCP[1]->getValues(t);
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii].LineBCP[0]->typeBCP == 3){
			SD.Lines[ii].pos.row(0)                = SD.Lines[ii].LineBCP[0]->pos.t();
			SD.Lines[ii].vel.row(0)                = SD.Lines[ii].LineBCP[0]->vel.t();
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP == 3){
			SD.Lines[ii].pos.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->pos.t();
			SD.Lines[ii].vel.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->vel.t();
		}
	}

	// Compute forces vector for the different objects, setting boundary conditions on acc if the BCP is not a joint
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		SD.Lines[ii].SEM_computeF();
		if(SD.Lines[ii].LineBCP[0]->typeBCP != 3){
			SD.Lines[ii].F.row(0)                = SD.Lines[ii].LineBCP[0]->acc.t();
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP != 3){
			SD.Lines[ii].F.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->acc.t();
		}

		if((SD.Lines[ii].LineBCP[0]->typeBCP != 3)&&(SD.Lines[ii].LineBCP[1]->typeBCP != 3)){
			SD.Lines[ii].acc = SD.Lines[ii].inv_MM_1N * SD.Lines[ii].F / SD.Lines[ii].dL;
		} else if ((SD.Lines[ii].LineBCP[0]->typeBCP == 3)&&(SD.Lines[ii].LineBCP[1]->typeBCP != 3)){
			SD.Lines[ii].acc = SD.Lines[ii].inv_MM_N * SD.Lines[ii].F / SD.Lines[ii].dL;
		} else if ((SD.Lines[ii].LineBCP[0]->typeBCP != 3)&&(SD.Lines[ii].LineBCP[1]->typeBCP == 3)){
			SD.Lines[ii].acc = SD.Lines[ii].inv_MM_1 * SD.Lines[ii].F / SD.Lines[ii].dL;
		} else if ((SD.Lines[ii].LineBCP[0]->typeBCP == 3)&&(SD.Lines[ii].LineBCP[1]->typeBCP == 3)){
			SD.Lines[ii].acc = SD.Lines[ii].inv_MM * SD.Lines[ii].F / SD.Lines[ii].dL;
		}
	}

	// Set boundary conditions on acc if the BCP is a joint
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii].LineBCP[0]->typeBCP == 3){
			SD.Lines[ii].LineBCP[0]->accLines.row(SD.Lines[ii].LineBCP[0]->iLJ) = SD.Lines[ii].acc.row(0);
			SD.Lines[ii].LineBCP[0]->iLJ = SD.Lines[ii].LineBCP[0]->iLJ + 1;
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP == 3){
			SD.Lines[ii].LineBCP[1]->accLines.row(SD.Lines[ii].LineBCP[1]->iLJ) = SD.Lines[ii].acc.row(SD.Lines[ii].N-1);
			SD.Lines[ii].LineBCP[1]->iLJ = SD.Lines[ii].LineBCP[1]->iLJ + 1;
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){		
		if(SD.Lines[ii].LineBCP[0]->typeBCP == 3){
			SD.Lines[ii].LineBCP[0]->getValues(t);
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP == 3){
			SD.Lines[ii].LineBCP[1]->getValues(t);
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii].LineBCP[0]->typeBCP == 3){
			SD.Lines[ii].acc.row(0)                = SD.Lines[ii].LineBCP[0]->acc.t();
		}
		if(SD.Lines[ii].LineBCP[1]->typeBCP == 3){
			SD.Lines[ii].acc.row(SD.Lines[ii].N-1) = SD.Lines[ii].LineBCP[1]->acc.t();
		}
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

int main () {

	int flag_read_eq, flag_write_eq;
	int nLines;
	int nBCPs, nAnchBCPs, nFairBCPs, nJointBCPs;
	int nNodosTotal, nSistema;
	int solver_flag, nIterMax;
	double atol, rtol;
	solver_data SD;

	double t;
	double t_max;
	double dt;

	PI=acos(-1.0);
	t = 0.0;

	nNodosTotal=0;

 	std::ifstream datosProblema ("input/datosProblema.dat");
	datosProblema >> g;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> rhoW;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> fondo;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> dt;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> t_max;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> solver_flag;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> atol;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> rtol;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> nIterMax;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_read_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_write_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema.close();


	//LEO DE FICHERO CUANTOS BCPs SE VAN A ESTUDIAR
	std::ifstream datosBCPs ("input/datosBCPs.dat");
	datosBCPs >> nFairBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nAnchBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nJointBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs.close();


	nBCPs = nFairBCPs + nAnchBCPs + nJointBCPs;

	BCP * BCPs [nBCPs];
	int BCPcounter = 0;

	BCP * F_BCPs = new FairleadBCP [nFairBCPs];
	for(int ii=0; ii<nFairBCPs; ii=ii+1){
		F_BCPs[ii].set_nBCP(BCPcounter+1);
		F_BCPs[ii].leer_datosBCPs();
		F_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &F_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}

	BCP * A_BCPs = new AnchorBCP [nAnchBCPs];
	for(int ii=0; ii<nAnchBCPs; ii=ii+1){
		A_BCPs[ii].set_nBCP(BCPcounter+1);
		A_BCPs[ii].leer_datosBCPs();
		A_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &A_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}

	BCP * J_BCPs = new JointBCP [nJointBCPs];
	for(int ii=0; ii<nJointBCPs; ii=ii+1){
		J_BCPs[ii].set_nBCP(BCPcounter+1);
		J_BCPs[ii].leer_datosBCPs();
		J_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &J_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}

	//LEO DE FICHERO Y PINTO EN PANTALLA CUANTAS LINEAS SE VAN A ESTUDIAR
	std::ifstream datosLines ("input/datosLines.dat");
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
	 	std::ifstream equi("input/Equilibrio.dat");
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

	if (solver_flag == 1){
		BDF S (t, t_max, dt, y, *fun, SD);
		S.atol = atol;
		S.rtol = rtol;
		S.nIterMax = nIterMax;
		std::cout<< "The temporal integration begins." << std::endl << std::endl;
		time_t tstart, tend; 
 		tstart = time(0);
		std::cout<< "    t = " << t << std::endl;
		do{
			S.step();
			//std::cout<< "    t = " << S.t << std::endl;
			if (S.t >= t + dt){
				t = t + dt;
				std::cout<< "    t = " << t << std::endl;
				for(int ii=0; ii<nLines; ii=ii+1) Lines[ii].write_out(S.t);
			}
		} while (S.t<=t_max);

		tend = time(0); 
		std::cout << " Computational time  : " << difftime(tend, tstart) << " seconds" << std::endl;
		std::cout << " Total function calls: " << nCalls << std::endl;
		std::cout << " Total jac calls: " << S.iJ << std::endl;

	}


	if (flag_write_eq == 1) {
		std::ofstream equi("output/Equilibrio.dat");
			for(int ii=0;ii<nSistema;ii=ii+1) equi << y(ii,0) << std::endl;
		equi.close();
	}

	return 0;
}