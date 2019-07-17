
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
	
	try
	{
		printf("Before initializing...\n");
		Simulation* mySim = new Simulation(project_path, "ASCII");
		printf("Before initializing...\n");
		mySim->LoadCase();
		mySim->Initialize();
		mySim->Run();
		printf("Water Depth: %f\n", mySim->waterDepth);
	}
	catch(Exception& error)
	{
		error.PrintDebug();
	}
	
	std::cout << "End of the program.\n"; //////////////////////////////////////////
	std::cout << "----------------------------------------------\n\n"; //////////////////////////////////////////
	return 0;
}