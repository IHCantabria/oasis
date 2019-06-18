
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include <armadillo>
#include <ctime>

#include "Simulations/Simulation.hpp"
#include "Exceptions/Exception.hpp"
#include "BCPs/BCPs.hpp"
#include "BCPs/Winchies.hpp"
#include "BCPs/WinchiesController.hpp"
#include "Lines/Lines.hpp"
#include "Spring/Spring.hpp"
#include "Bodies/Bodies.hpp"
#include "Hydro/Hydro.hpp"
#include "ODE_solvers/ODE_solvers.hpp"
#include "os_tools.hpp"

// GLOBAL VARIABLES
double PI;
double g;
double rhoW;
double fondo;
int nCalls = 0;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////           fun           ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::mat fun(double t, arma::mat y, solver_data SD){

	nCalls = nCalls + 1;

	arma::mat yprime = arma::zeros(size(y));
	int i0;

	// Copy info from y to the objects.
	int ini = 0;
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii].pos = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii].N;jj=jj+1){
			SD.Lines[ii].pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	ini = ini + SD.nWinchies;
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii].vel = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii].N;jj=jj+1){
			SD.Lines[ii].vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<SD.nWinchies;ii=ii+1){
		SD.Winchies[ii].theta = arma::as_scalar(y.row(SD.nSistema2-(ii+1)));
		SD.Winchies[ii].omega = arma::as_scalar(y.row(SD.nSistema-(ii+1)));
	}
	

	// Update BodyBCP positions and velocities
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii].updateBCPs();
	}

	// Set boundary conditions on pos and vel of Lines if the BCP is not a joint
	if (SD.nLines >= 1){
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
	// Set boundary conditions on pos and vel of Lines if the BCP is a joint
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

	// Compute forces vector for the different Lines
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		SD.Lines[ii].SEM_computeF();
	}

	// Compute forces of Springs
	for(int ii=0;ii<SD.nSprings;ii=ii+1){
		SD.Springs[ii].computeSpringForces();
	}

	// Compute hydrostatic and hidrodynamic forces
	SD.Water->computeHydroForces(t);
	arma::mat Fb = SD.Water->HydroForces;
	// Compute forces on BCPs
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii].computeBCPForces();
		Fb(arma::span(6*ii,6*(ii+1)-1), 0) = Fb(arma::span(6*ii,6*(ii+1)-1), 0) + SD.Bodies[ii].BCPForces;
	}

	// Compute body acceleration
	arma::mat accB = (SD.Water->invM)*Fb;
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii].acc = accB(arma::span(6*ii,6*(ii+1)-1), 0);
	}

	// Update BodyBCP accelerations
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii].updateBCPs();
	}

	// Obtain Lines accelerations, imposing boundary conditions if the BCP is not a joint
	for(int ii=0;ii<SD.nLines;ii=ii+1){
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
	// Obtain Lines accelerations, imposing boundary conditions if the BCP is a joint
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

	// Compute Winchies
	for(int ii=0;ii<SD.nWinchies;ii=ii+1){
		SD.Winchies[ii].computeWinchie();
	}

	// Copy info from the objects to yprime
	yprime.rows(0,SD.nSistema2-1) = y.rows(SD.nSistema2,SD.nSistema-1);
	yprime.rows(SD.nSistema2,SD.nSistema2+6*SD.nBodies-1) = accB;
	ini = SD.nSistema2+6*SD.nBodies;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii].N;jj=jj+1){
			yprime.rows(ini,ini+2) = SD.Lines[ii].acc.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<SD.nWinchies;ii=ii+1){
		yprime.row(SD.nSistema2-(ii+1)) = SD.Winchies[ii].omega;
		yprime.row(SD.nSistema-(ii+1)) = SD.Winchies[ii].alpha;
	}


	if (yprime.has_nan()){
		std::cout << std::endl << "ERROR: NaN Detected!" << std::endl;
		throw std::exception();
	}
	return yprime;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////          MAIN           ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char* argv[]) 
{
	std::cout << std::endl << "----------------------------------------------"  << std::endl; //////////////////////////////////////////
	std::cout << "Starting OASIS: Offshore Advanced Simulation Software" << std::endl << std::endl; //////////////////////////////////////////

	// Declare variables
	int flag_read_eq, flag_write_eq;
	int numLines;
	int nSprings;
	int nWinchies;
	int nBodies;
	int nBCPs, nAnchBCPs, nFairBCPs, nJointBCPs, nBodyBCPs;
	int nNodosTotal, nSistema, nSistema2;
	int solver_flag, nIterMax;
	double atol, rtol;
	solver_data SD;
	double t;
	double t_max;
	double dt;
	std::string project_path;
	std::string inputs_path;
	std::string outputs_path;
	std::string file_path;

	// Read input arguments
	if (argc < 2)
	{
		printf("Not enought input arguments. First argument must be the project root path.");
		return 1;
	}
	else
	{
		project_path = argv[1];
		printf("PROJECT ROOT PATH: %s\n\n", project_path.c_str());
	}

	PI=acos(-1.0);
	t = 0.0;
	nNodosTotal=0;

	try
	{
		Simulation* mySim = new Simulation(project_path, "ASCII");
		mySim->ReadProperties();
		printf("Water Depth: %f\n", mySim->waterDepth);
	}
	catch(Exception& error)
	{
		error.PrintDebug();
	}
	
	

	/**
	std::cout << "  Reading datosProblema.dat ..." << std::endl << std::endl; //////////////////////////////////////////

	file_path = JoinPath(project_path, "input/datosProblema.dat");
 	std::ifstream datosProblema(file_path.c_str());
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

	std::cout<< "  Reading and starting all BCPs ..." << std::endl << std::endl; //////////////////////////////////////////

	//LEO DE FICHERO CUANTOS BCPs SE VAN A ESTUDIAR
	file_path = JoinPath(project_path, "input/datosBCPs.dat");
	std::ifstream datosBCPs (file_path.c_str());
	datosBCPs >> nFairBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nAnchBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nJointBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> nBodyBCPs; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs.close();
	nBCPs = nFairBCPs + nAnchBCPs + nJointBCPs + nBodyBCPs;
	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA BCP
	BCP * BCPs [nBCPs];
	int BCPcounter = 0;
	//INICIO LOS BCPs
	BCP * F_BCPs = new FairleadBCP [nFairBCPs];
	for(int ii=0; ii<nFairBCPs; ii=ii+1){
		F_BCPs[ii].set_nBCP(BCPcounter+1);
		F_BCPs[ii].leer_datosBCPs(file_path);
		F_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &F_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}
	BCP * A_BCPs = new AnchorBCP [nAnchBCPs];
	for(int ii=0; ii<nAnchBCPs; ii=ii+1){
		A_BCPs[ii].set_nBCP(BCPcounter+1);
		A_BCPs[ii].leer_datosBCPs(file_path);
		A_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &A_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}
	BCP * J_BCPs = new JointBCP [nJointBCPs];
	for(int ii=0; ii<nJointBCPs; ii=ii+1){
		J_BCPs[ii].set_nBCP(BCPcounter+1);
		J_BCPs[ii].leer_datosBCPs(file_path);
		J_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &J_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}
	BCP * B_BCPs = new BodyBCP [nBodyBCPs];
	for(int ii=0; ii<nBodyBCPs; ii=ii+1){
		B_BCPs[ii].set_nBCP(BCPcounter+1);
		B_BCPs[ii].leer_datosBCPs(file_path);
		BCPs[BCPcounter] = &B_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}

	std::cout<< "  Reading and starting all Bodies ..." << std::endl << std::endl; //////////////////////////////////////////

	//LEO DE FICHERO CUANTOS CUERPOS SE VAN A ESTUDIAR
	file_path = JoinPath(project_path, "input/datosBodies.dat");
	std::ifstream datosBodies (file_path.c_str());
	datosBodies >> nBodies; datosBodies.ignore(std::numeric_limits<int>::max(), '\n');
	datosBodies.close();
	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA CUERPO
	Body * Bodies = new Body[nBodies];
	//INICIO LOS CUERPOS
	for(int ii=0; ii<nBodies; ii=ii+1){
		Bodies[ii].set_nBody(ii+1);
		Bodies[ii].leer_datosBody(project_path);
		for(int jj=0; jj<Bodies[ii].nBCPs; jj=jj+1){
			Bodies[ii].BodyBCPs[jj] = BCPs[Bodies[ii].index_BCPs[jj]-1];
		}
		Bodies[ii].updateBCPs();
	}

	for(int ii=0; ii<nBCPs; ii=ii+1){
		BCPs[ii]->getValues(0.0);
	}

	std::cout<< "  Reading and starting all Lines ..." << std::endl << std::endl; //////////////////////////////////////////

	// Read Lines properties
	FILE* fid_lines;
	file_path = JoinPath(project_path, "input/datosLines.dat");
	printf("Lines files path: %s\n", file_path.c_str());
	fid_lines = fopen(file_path.c_str(), "r");
	std::ifstream datosLines (file_path.c_str());
	datosLines >> numLines; datosLines.ignore(std::numeric_limits<int>::max(), '\n');
	datosLines.close();
	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA LINEA
	Line * Lines = new Line[numLines];
	//INICIO LAS LINEAS
	for(int ii=0; ii<numLines; ii=ii+1){
		Lines[ii].set_nLine(ii+1);
		printf("aquiii\n");
		try
		{
		Lines[ii].leer_datosLines(file_path);
		Lines[ii].print_out();
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

	std::cout<< "  Reading and starting all Winchies ..." << std::endl << std::endl; //////////////////////////////////////////

	//LEO DE FICHERO CUANTOS WINCHIES SE VAN A ESTUDIAR
	file_path = JoinPath(project_path, "input/datosWinchies.dat");
	std::ifstream datosWinchies (file_path);
	datosWinchies >> nWinchies; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');
	datosWinchies.close();
	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA MUELLE
	Winchie * Winchies = new Winchie[nWinchies];
	//INICIO LLOS MUELLES
	for(int ii=0; ii<nWinchies; ii=ii+1){
		Winchies[ii].set_nWinchie(ii+1);
		Winchies[ii].leer_datosWinchies();
		Winchies[ii].LineW = &Lines[Winchies[ii].nLine - 1];
	}
	fclose(fid_lines);
	std::cout<< "  Reading and starting Winchies controllers ... " << std::endl << std::endl; //////////////////////////////////////////

	// DEFINO UN UN OBJETO DE LA CLASE WINCHIE CONTROLLER
	WinchieController CW;
	//INICIO EL OBJETO QUE CONTIENE EL CONTROLADOR
	CW.set_WinchieController(nWinchies,Winchies);
	CW.leer_datosWinchieController();


	std::cout<< "  Reading and starting all Springs ..." << std::endl << std::endl; //////////////////////////////////////////

	//LEO DE FICHERO CUANTOS MUELLES SE VAN A ESTUDIAR
	file_path = JoinPath(project_path, "input/datosSprings.dat");
	std::ifstream datosSprings (file_path);
	datosSprings >> nSprings; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings.close();
	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA MUELLE
	Spring * Springs = new Spring[nSprings];
	//INICIO LLOS MUELLES
	for(int ii=0; ii<nSprings; ii=ii+1){
		Springs[ii].set_nSpring(ii+1);
		Springs[ii].leer_datosSprings(file_path);
		Springs[ii].SpringBCP[0] = BCPs[Springs[ii].BCP_1-1];
		Springs[ii].SpringBCP[1] = BCPs[Springs[ii].BCP_2-1];
	}

	std::cout<< "  Reading and starting Hyrodynamics ... " << std::endl << std::endl; //////////////////////////////////////////

	// DEFINO UN POINTER A UN OBJETO DE LA CLASE HYDRO
	file_path = JoinPath(project_path, "input/flotante.h5");
	Hydro * Water = new Hydro;
	//INICIO EL OBJETO QUE CONTIENE LA HIDRODINAMICA
	Water->set_Hydro(nBodies,Bodies);
	Water->leer_datosHydro(file_path);
	Water->computeIRF();
	Water->computeWaveSpectrum();
	Water->computeFe();

	std::cout<< "End of input reading" << std::endl << std::endl; //////////////////////////////////////////


	std::cout<< "  Initializing ODE system vector ..." << std::endl << std::endl; //////////////////////////////////////////
	// Inicio el vector del sistema
	nSistema = 2*(3*nNodosTotal + 6*nBodies + nWinchies);
	nSistema2 = 3*nNodosTotal + 6*nBodies + nWinchies;
	SD.nSistema = nSistema;
	SD.nSistema2 = nSistema2;

	arma::mat y = arma::zeros(nSistema,1);
	arma::mat yprime = arma::zeros(nSistema,1);

	int ini=0;
	if(flag_read_eq==0){
		for(int ii=0; ii<nBodies; ii=ii+1){
				y.rows(ini,ini+5) = Bodies[ii].pos;
				ini=ini+6;
		}
		for(int ii=0; ii<numLines; ii=ii+1){
			for(int jj=0; jj<Lines[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = Lines[ii].pos.row(jj).t();
					ini=ini+3;
			}
		}
	}else if (flag_read_eq==1){
		file_path = JoinPath(project_path, "input/Equilibrio.dat");
	 	std::ifstream equi(file_path);
			for(int ii=6*nBodies;ii<nSistema2;ii=ii+1) {
				equi >> y(ii,0);    equi.ignore(std::numeric_limits<int>::max(), '\n');
			}
		equi.close();

		ini=6*nBodies;
		for(int ii=0; ii<numLines; ii=ii+1){
			for(int jj=0; jj<Lines[ii].N; jj=jj+1){
				Lines[ii].pos.row(jj) = y.rows(ini,ini+2).t();
				ini=ini+3;
			}
		}
	}

	// ESCIBIENDO CONDICION INICIAL A FICHEROS
	for(int ii=0; ii<numLines; ii=ii+1) Lines[ii].write_out(t);
	for(int ii=0; ii<nBodies; ii=ii+1) Bodies[ii].write_out(t);

	// GUARDANDO DATOS LEIDOS EN ESTRUCTURA DEL SOLVER TEMPORAL
	SD.nLines = numLines;
	SD.Lines = Lines;
	SD.nSprings = nSprings;
	SD.Springs = Springs;
	SD.nBodies = nBodies;
	SD.Bodies = Bodies;
	SD.Water = Water;
	SD.nWinchies = nWinchies;
	SD.Winchies = Winchies;

	std::cout<< "  Starting temporal integration ..." << std::endl << std::endl; //////////////////////////////////////////

	if (solver_flag == 1){
		BDF S (t, t_max, dt, y, *fun, SD);
		S.atol = atol;
		S.rtol = rtol;
		S.nIterMax = nIterMax;
		time_t tstart, tend; 
 		tstart = time(0);
		std::cout<< "    t = " << t << " s" << std::endl;
		do{
			S.step();
			if (S.t >= t + dt){
				t = t + dt;
				std::cout<< "    t = " << t << " s"  << std::endl;
				if (nWinchies>0) CW.controlWinchies();
				for(int ii=0; ii<numLines; ii=ii+1) Lines[ii].write_out(S.t);
				for(int ii=0; ii<nBodies; ii=ii+1) Bodies[ii].write_out(S.t);
			}
		} while (S.t<=t_max);

		tend = time(0); 
		std::cout << std::endl << "    Computational time  : " << difftime(tend, tstart) << " seconds" << std::endl;
		std::cout << "    Total function calls: " << nCalls << std::endl;
		std::cout << "    Total jac calls: " << S.iJ << std::endl << std::endl;
	}


	if (flag_write_eq == 1) {

		std::cout << "  Writting data to Equilibrio.dat ..." << std::endl << std::endl; //////////////////////////////////////////

		std::ofstream equi("output/Equilibrio.dat");
			for(int ii=6*nBodies;ii<nSistema2;ii=ii+1) equi << y(ii,0) << std::endl;
		equi.close();
	}
	**/
	std::cout << "End of the program." << std::endl; //////////////////////////////////////////
	std::cout << "----------------------------------------------" << std::endl << std::endl; //////////////////////////////////////////
	return 0;
}