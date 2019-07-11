
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
#include "Hydro/HydroDatabase.hpp"
#include "ODE_solvers/ODE_solvers.hpp"
#include "os_tools.hpp"

// GLOBAL VARIABLES
double PI;
int nCalls = 0;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////           fun           ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
arma::mat fun(double t, arma::mat y, solver_data& SD){

	nCalls = nCalls + 1;
	//std::cout << "Main::fun - At first" << std::endl; 
	arma::mat yprime = arma::zeros(size(y));
	int i0;
	//std::cout << "Main::fun - Before copy y to objects" << std::endl;
	// Copy info from y to the objects.
	int ini = 0;
	for(int ii=0;ii<SD.nBodies;ii=ii+1)
	{
		SD.Bodies[ii]->pos = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1)
	{
		for(int jj=0;jj<SD.Lines[ii]->N;jj=jj+1){
			SD.Lines[ii]->pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	ini = ini + SD.nWinchies;
	for(int ii=0;ii<SD.nBodies;ii=ii+1)
	{
		SD.Bodies[ii]->vel = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1)
	{
		for(int jj=0;jj<SD.Lines[ii]->N;jj=jj+1)
		{
			SD.Lines[ii]->vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<SD.nWinchies;ii=ii+1){
		SD.Winchies[ii]->theta = arma::as_scalar(y.row(SD.nSistema2-(ii+1)));
		SD.Winchies[ii]->omega = arma::as_scalar(y.row(SD.nSistema-(ii+1)));
	}
	
	// Update BodyBCP positions and velocities
	//std::cout << "Main::fun - Update BCP positions and velocities" << std::endl;
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii]->UpdateBcps();
	}
	// Set boundary conditions on pos and vel of Lines if the BCP is not a joint
	//std::cout << "Main::fun - Set Boundary conditios" << std::endl;
	if (SD.nLines >= 1){
		if (SD.Lines[0]->pLineBcps[0]->tBCP != t){
			for(int ii=0;ii<SD.nLines;ii=ii+1){		
				if(SD.Lines[ii]->pLineBcps[0]->GetType() != 3){
					SD.Lines[ii]->pLineBcps[0]->GetValues(t);
				}
				if(SD.Lines[ii]->pLineBcps[1]->GetType() != 3){
					SD.Lines[ii]->pLineBcps[1]->GetValues(t);
				}
			}
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() != 3){
			SD.Lines[ii]->pos.row(0)                = SD.Lines[ii]->pLineBcps[0]->pos.t();
			SD.Lines[ii]->vel.row(0)                = SD.Lines[ii]->pLineBcps[0]->vel.t();
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() != 3){
			SD.Lines[ii]->pos.row(SD.Lines[ii]->N-1) = SD.Lines[ii]->pLineBcps[1]->pos.t();
			SD.Lines[ii]->vel.row(SD.Lines[ii]->N-1) = SD.Lines[ii]->pLineBcps[1]->vel.t();
		}
	}
	// Set boundary conditions on pos and vel of Lines if the BCP is a joint
	//std::cout << "Main::fun - Set Boundary conditios if BCP is a Joint" << std::endl;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[0]->posLines.row(SD.Lines[ii]->pLineBcps[0]->iLJ) = SD.Lines[ii]->pos.row(0);
			SD.Lines[ii]->pLineBcps[0]->velLines.row(SD.Lines[ii]->pLineBcps[0]->iLJ) = SD.Lines[ii]->vel.row(0);
			SD.Lines[ii]->pLineBcps[0]->iLJ = SD.Lines[ii]->pLineBcps[0]->iLJ + 1;
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[1]->posLines.row(SD.Lines[ii]->pLineBcps[1]->iLJ) = SD.Lines[ii]->pos.row(SD.Lines[ii]->N-1);
			SD.Lines[ii]->pLineBcps[1]->velLines.row(SD.Lines[ii]->pLineBcps[1]->iLJ) = SD.Lines[ii]->vel.row(SD.Lines[ii]->N-1);
			SD.Lines[ii]->pLineBcps[1]->iLJ = SD.Lines[ii]->pLineBcps[1]->iLJ + 1;
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[0]->GetValues(t);
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[1]->GetValues(t);
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() == 3){
			SD.Lines[ii]->pos.row(0)                = SD.Lines[ii]->pLineBcps[0]->pos.t();
			SD.Lines[ii]->vel.row(0)                = SD.Lines[ii]->pLineBcps[0]->vel.t();
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() == 3){
			SD.Lines[ii]->pos.row(SD.Lines[ii]->N-1) = SD.Lines[ii]->pLineBcps[1]->pos.t();
			SD.Lines[ii]->vel.row(SD.Lines[ii]->N-1) = SD.Lines[ii]->pLineBcps[1]->vel.t();
		}
	}
	// Compute forces vector for the different Lines
	//std::cout << "Main::fun - Compute forces vector for different lines" << std::endl;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		SD.Lines[ii]->SEM_computeF();
	}
	// Compute forces of Springs
	//std::cout << "Main::fun - Compute spring" << std::endl;
	for(int ii=0;ii<SD.nSprings;ii=ii+1){
		SD.Springs[ii]->computeSpringForces();
	}
	// Compute hydrostatic and hidrodynamic forces
	//std::cout << "Main::fun - Compute hydrodynamic and hydrostatic forces" << std::endl;
	arma::mat Fb = SD.Water->ComputeHydrostaticForces();
	arma::mat Fr = SD.Water->ComputeRadiationForces(t, SD);
	arma::mat Fe = SD.Water->ComputeFirstWaveExcForce(t);
	// Compute forces on BCPs
	//std::cout << "Main::fun - Compute forces on BCPs" << std::endl;
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii]->ComputeBcpForces();
		Fb(arma::span(6*ii,6*(ii+1)-1), 0) =  Fe - (Fr + Fb(arma::span(6*ii,6*(ii+1)-1), 0)) + SD.Bodies[ii]->bcpForces;
	}
	WriteASCII("output/bcpForces.dat", SD.Bodies[0]->bcpForces, true);
	// Compute body acceleration
	//(**SD.Water->pAddedMassHf).print();
	//arma::mat accB = (*SD.Water->pStructuralMass+**SD.Water->pAddedMassHf)*Fb;
	//std::cout << "Main::fun - Compute accelerations" << std::endl;
	arma::mat accB;
	arma::mat total_mass;
	WriteASCII("output/accB.dat", accB, true);
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		total_mass = *SD.Water->pStructuralMass+(*SD.Water->pAddedMassHf[ii]);
		//std::cout << "Total mass matrix(2,2): " << mmm(2,2) << std::endl;
		accB = arma::solve(total_mass,Fb);
		//accB = arma::solve(*SD.Water->pStructuralMass,Fb);
		SD.Bodies[ii]->acc = accB(arma::span(6*ii,6*(ii+1)-1), 0);
	}

	// Update BodyBCP accelerations
	//std::cout << "Main::fun - Compute BCP accelerations" << std::endl;
	for(int ii=0;ii<SD.nBodies;ii=ii+1){
		SD.Bodies[ii]->UpdateBcps();
	}
	//WriteASCII("output/bcpAcceleration.dat", SD.Bodies[0]->pBodyBcps[0]->accG_BCP, true);
	// Obtain Lines accelerations, imposing boundary conditions if the BCP is not a joint
	//std::cout << "Main::fun - Compute lines accelerations" << std::endl;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() != 3){
			SD.Lines[ii]->F.row(0)                = SD.Lines[ii]->pLineBcps[0]->acc.t();
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() != 3){
			SD.Lines[ii]->F.row(SD.Lines[ii]->N-1) = SD.Lines[ii]->pLineBcps[1]->acc.t();
		}
		if((SD.Lines[ii]->pLineBcps[0]->GetType() != 3)&&(SD.Lines[ii]->pLineBcps[1]->GetType() != 3)){
			SD.Lines[ii]->acc = SD.Lines[ii]->inv_MM_1N * SD.Lines[ii]->F / SD.Lines[ii]->dL;
		} else if ((SD.Lines[ii]->pLineBcps[0]->GetType() == 3)&&(SD.Lines[ii]->pLineBcps[1]->GetType() != 3)){
			SD.Lines[ii]->acc = SD.Lines[ii]->inv_MM_N * SD.Lines[ii]->F / SD.Lines[ii]->dL;
		} else if ((SD.Lines[ii]->pLineBcps[0]->GetType() != 3)&&(SD.Lines[ii]->pLineBcps[1]->GetType() == 3)){
			SD.Lines[ii]->acc = SD.Lines[ii]->inv_MM_1 * SD.Lines[ii]->F / SD.Lines[ii]->dL;
		} else if ((SD.Lines[ii]->pLineBcps[0]->GetType() == 3)&&(SD.Lines[ii]->pLineBcps[1]->GetType() == 3)){
			SD.Lines[ii]->acc = SD.Lines[ii]->inv_MM * SD.Lines[ii]->F / SD.Lines[ii]->dL;
		}
	}
	// Obtain Lines accelerations, imposing boundary conditions if the BCP is a joint
	//std::cout << "Main::fun - Lines accelerations if the Line is a Joint" << std::endl;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[0]->accLines.row(SD.Lines[ii]->pLineBcps[0]->iLJ) = SD.Lines[ii]->acc.row(0);
			SD.Lines[ii]->pLineBcps[0]->iLJ = SD.Lines[ii]->pLineBcps[0]->iLJ + 1;
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[1]->accLines.row(SD.Lines[ii]->pLineBcps[1]->iLJ) = SD.Lines[ii]->acc.row(SD.Lines[ii]->N-1);
			SD.Lines[ii]->pLineBcps[1]->iLJ = SD.Lines[ii]->pLineBcps[1]->iLJ + 1;
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[0]->GetValues(t);
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() == 3){
			SD.Lines[ii]->pLineBcps[1]->GetValues(t);
		}
	}
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		if(SD.Lines[ii]->pLineBcps[0]->GetType() == 3){
			SD.Lines[ii]->acc.row(0)                = SD.Lines[ii]->pLineBcps[0]->acc.t();
		}
		if(SD.Lines[ii]->pLineBcps[1]->GetType() == 3){
			SD.Lines[ii]->acc.row(SD.Lines[ii]->N-1) = SD.Lines[ii]->pLineBcps[1]->acc.t();
		}
	}

	// Compute Winchies
	//std::cout << "Main::fun - Compute Winchies" << std::endl;
	for(int ii=0;ii<SD.nWinchies;ii=ii+1){
		SD.Winchies[ii]->computeWinchie();
	}

	// Copy info from the objects to yprime
	//std::cout << "Main::fun - Copy info to yprime" << std::endl;
	yprime.rows(0,SD.nSistema2-1) = y.rows(SD.nSistema2,SD.nSistema-1);
	yprime.rows(SD.nSistema2,SD.nSistema2+6*SD.nBodies-1) = accB;
	ini = SD.nSistema2+6*SD.nBodies;
	for(int ii=0;ii<SD.nLines;ii=ii+1){
		for(int jj=0;jj<SD.Lines[ii]->N;jj=jj+1){
			yprime.rows(ini,ini+2) = SD.Lines[ii]->acc.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<SD.nWinchies;ii=ii+1){
		yprime.row(SD.nSistema2-(ii+1)) = SD.Winchies[ii]->omega;
		yprime.row(SD.nSistema-(ii+1)) = SD.Winchies[ii]->alpha;
	}
	//std::cout << "Main::fun - Check if yprime has a NaN" << std::endl;
	WriteASCII("output/yprime.dat", yprime, true);
	if (yprime.has_nan()){
		std::cout << std::endl << "ERROR: NaN Detected!" << std::endl;
		throw std::exception();
	}
	//std::cout << "Main::fun - End of fcn" << std::endl;
	return yprime;
	
}
**/


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

	// PI
	PI=acos(-1.0);
	t = 0.0;

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
	printf("Before initializing...\n");
	PI=acos(-1.0);
	t = 0.0;
	nNodosTotal=0;
	printf("Before initializing...\n");
	try
	{
		printf("Before initializing...\n");
		Simulation* mySim = new Simulation(project_path, "ASCII");
		printf("Before initializing...\n");
		mySim->LoadCase();
		mySim->Initialize();
		mySim->Run();
		printf("Water Depth: %f\n", mySim->waterDepth);
	//}
	//catch(Exception& error)
	//{
//		error.PrintDebug();
	//}
	
	

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
	BCP* BCPs [nBCPs];
	int BCPcounter = 0;

	//INICIO LOS BCPs
	BCP* F_BCPs = new FairleadBCP [nFairBCPs];
	for(int ii=0; ii<nFairBCPs; ii=ii+1){
		F_BCPs[ii].set_nBCP(BCPcounter+1);
		F_BCPs[ii].leer_datosBCPs(file_path);
		F_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &F_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}
	BCP* A_BCPs = new AnchorBCP [nAnchBCPs];
	for(int ii=0; ii<nAnchBCPs; ii=ii+1){
		A_BCPs[ii].set_nBCP(BCPcounter+1);
		A_BCPs[ii].leer_datosBCPs(file_path);
		A_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &A_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}
	BCP* J_BCPs = new JointBCP [nJointBCPs];
	for(int ii=0; ii<nJointBCPs; ii=ii+1){
		J_BCPs[ii].set_nBCP(BCPcounter+1);
		J_BCPs[ii].leer_datosBCPs(file_path);
		J_BCPs[ii].getValues(0.0);
		BCPs[BCPcounter] = &J_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}
	BCP* B_BCPs = new BodyBCP [nBodyBCPs];
	for(int ii=0; ii<nBodyBCPs; ii=ii+1){
		B_BCPs[ii].set_nBCP(BCPcounter+1);
		B_BCPs[ii].leer_datosBCPs(file_path);
		BCPs[BCPcounter] = &B_BCPs[ii];
		BCPcounter = BCPcounter + 1;
	}

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
	**/
		// DEFINO UN POINTER A UN OBJETO DE LA CLASE HYDRO
		
		//Water->computeIRF();
		//Water->computeWaveSpectrum();
		//Water->computeFe();

		//std::cout<< "End of input reading" << std::endl << std::endl; //////////////////////////////////////////

		
		//std::cout<< "  Initializing ODE system vector ..." << std::endl << std::endl; //////////////////////////////////////////
		// Inicio el vector del sistema
		/**
		std::cout << "Num DOFs Total: " << mySim->numDofTotal << std::endl;
		std::cout << "Num NumBodies: " << mySim->numBodies << std::endl;
		std::cout << "Num NumWinchies: " << mySim->numWinches << std::endl;
		nSistema2 = 3*mySim->numDofTotal + 6*mySim->numBodies + mySim->numWinches;
		nSistema = 2*nSistema2;
		SD.nSistema = nSistema;
		SD.nSistema2 = nSistema2;

		printf("Sistema size: %d\n", nSistema);
		printf("Sistema2 size: %d\n", nSistema2);

		arma::mat y = arma::zeros(nSistema,1);
		arma::mat yprime = arma::zeros(nSistema,1);

		int ini=0;
		if(flag_read_eq==0){
			for(int ii=0; ii<mySim->numBodies; ii=ii+1){
					y.rows(ini,ini+5) = mySim->pBodies[ii]->pos;
					ini=ini+6;
			}
			for(int ii=0; ii<mySim->numLines; ii=ii+1){
				for(int jj=0; jj<mySim->pLines[ii]->N; jj=jj+1){
						y.rows(ini,ini+2) = mySim->pLines[ii]->pos.row(jj).t();
						ini=ini+3;
				}
			}
		}else if (mySim->readEquilibrium){
			file_path = JoinPath(mySim->inputFolderPath, "Equilibrio.dat");
			std::ifstream equi(file_path);
				for(int ii=6*mySim->numBodies;ii<nSistema2;ii=ii+1) {
					equi >> y(ii,0);    equi.ignore(std::numeric_limits<int>::max(), '\n');
				}
			equi.close();

			ini=6*mySim->numBodies;
			for(int ii=0; ii<mySim->numLines; ii=ii+1){
				for(int jj=0; jj<mySim->pLines[ii]->N; jj=jj+1){
					mySim->pLines[ii]->pos.row(jj) = y.rows(ini,ini+2).t();
					ini=ini+3;
				}
			}
		}
		

		// ESCIBIENDO CONDICION INICIAL A FICHEROS
		for(int ii=0; ii<mySim->numLines; ii=ii+1) mySim->pLines[ii]->write_out(t);
		for(int ii=0; ii<mySim->numBodies; ii=ii+1) mySim->pBodies[ii]->WriteOut(t);
		**/

		// GUARDANDO DATOS LEIDOS EN ESTRUCTURA DEL SOLVER TEMPORAL
		/**
		SD.nLines = mySim->numLines;
		SD.Lines = mySim->pLines;
		SD.nSprings = mySim->numSprings;
		SD.Springs = mySim->pSprings;
		SD.nBodies = mySim->numBodies;
		SD.Bodies = mySim->pBodies;
		SD.Water = mySim->pBodies[0]->hydro;
		SD.nWinchies = mySim->numWinches;
		SD.Winchies = mySim->pWinches;
		arma::mat* timeBufferNew;
		SD.timeBufferCount = new int(0);
		SD.timeBuffer = new arma::mat(1, SD.timeBufferSize, arma::fill::zeros);
		(*SD.timeBuffer)(0,0) = 1;

		// Save first data
		std::cout << "Antes de update system" << std::endl;
		mySim->UpdateSystem(0.0, y, SD);
		// Update time vector if any
		if ((*SD.timeBufferCount) < SD.timeBufferSize)
		{
			(*SD.timeBuffer)(0, (*SD.timeBufferCount)) = 0.0;
		}
		else
		{
			timeBufferNew = new arma::mat(1, SD.timeBufferSize, arma::fill::zeros);
			(*timeBufferNew).cols(0, mySim->pBodies[0]->hydro->numPointsIRF-1) = (*SD.timeBuffer).cols(SD.timeBufferSize-mySim->pBodies[0]->hydro->numPointsIRF, SD.timeBufferSize-1);
			SD.timeBuffer = timeBufferNew;
			(*SD.timeBufferCount) = mySim->pBodies[0]->hydro->numPointsIRF-1;
			delete timeBufferNew;
		}
		(*SD.timeBufferCount)++;

		std::cout<< "  Starting temporal integration ..." << std::endl << std::endl; //////////////////////////////////////////
		std::cout<< "  Starting temporal integration ... " << mySim->timeIntMethod << std::endl;
		std::cout<< "  Starting temporal integration ... " << mySim->simulationTime << std::endl;
		//if (mySim->timeIntMethod == 1)
		if (mySim->timeIntMethod == 1)
		{
			std::cout << "Initializing temporal solver..." << std::endl;
			BDF S (t, mySim->simulationTime, mySim->maxTimeStep, y, *fun, SD);
			S.atol = mySim->timeIntAbsTol;
			S.rtol = mySim->timeIntRelTol;
			S.nIterMax = mySim->maxIterStep;
			time_t tstart, tend;
			
			tstart = time(0);
			std::cout<< "    t = " << t << " s" << std::endl;
			do{
				std::cout << "Antes de Step" << std::endl;
				S.step();
				
				// Update time vector if any
				std::cout << "Antes de Update Time" << std::endl;
				if ((*SD.timeBufferCount) < SD.timeBufferSize)
				{
					std::cout << "New time step" << std::endl;
					(*SD.timeBuffer)(0, (*SD.timeBufferCount)) = S.t;
					std::cout << "New time step --> done" << std::endl;
				}
				else
				{
					timeBufferNew = new arma::mat(1, SD.timeBufferSize, arma::fill::zeros);
					(*timeBufferNew).cols(0, mySim->pBodies[0]->hydro->numPointsIRF-1) = (*SD.timeBuffer).cols(SD.timeBufferSize-mySim->pBodies[0]->hydro->numPointsIRF, SD.timeBufferSize-1);
					SD.timeBuffer = timeBufferNew;
					(*SD.timeBufferCount) = mySim->pBodies[0]->hydro->numPointsIRF-1;
					delete timeBufferNew;
				}
				(*SD.timeBufferCount)++;
				std::cout << "Antes de Update System" << std::endl;
				mySim->UpdateSystem(S.t, S.y, SD);
				std::cout << "Despues de Update System" << std::endl;
				
				// Print out time if any
				if (S.t >= t + mySim->maxTimeStep)
				{
					t = t + mySim->maxTimeStep;
					std::cout<< "    t = " << t << " s"  << std::endl;
					//if (nWinchies>0) CW.controlWinchies();
					for(int ii=0; ii<mySim->numLines; ii=ii+1) mySim->pLines[ii]->WriteOut(S.t);
					for(int ii=0; ii<mySim->numBodies; ii=ii+1) mySim->pBodies[ii]->WriteOut(S.t);
				}
			} while (S.t<=mySim->simulationTime);
			
			tend = time(0); 
			std::cout << std::endl << "    Computational time  : " << difftime(tend, tstart) << " seconds" << std::endl;
			std::cout << "    Total function calls: " << numCallsSysFun << std::endl;
			std::cout << "    Total jac calls: " << S.iJ << std::endl << std::endl;
			
			
		}

		if (flag_write_eq == 1) {

			std::cout << "  Writting data to Equilibrio.dat ..." << std::endl << std::endl; //////////////////////////////////////////

			std::ofstream equi("output/Equilibrio.dat");
				for(int ii=6*nBodies;ii<nSistema2;ii=ii+1) equi << y(ii,0) << std::endl;
			equi.close();
		}
		**/
	}
	catch(Exception& error)
	{
		error.PrintDebug();
	}
	
	std::cout << "End of the program." << std::endl; //////////////////////////////////////////
	std::cout << "----------------------------------------------" << std::endl << std::endl; //////////////////////////////////////////
	return 0;
}