
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
#include <filesystem>

int main(int argc, char *argv[])
{
	std::cout << std::endl
			  << "----------------------------------------------" << std::endl;
	std::cout << "OASIS: Offshore Advanced Simulation Software" << std::endl
			  << std::endl;

	// Declare variables
	std::string project_path;

	// Read input arguments
	if (argc < 2)
	{
		std::cout << "Not enought input arguments. First argument must be the project root path." << std::endl;
		return 1;
	}
	else
	{
		project_path = argv[1];
		std::cout << "Simulation directory path:\n";
		if (project_path.length() > 0)
		{
			std::cout << "  " << project_path << std::endl;
		}
		else
		{
			std::cout << "  /.  (pwd)" << std::endl;
		}
	}

	try
	{
		// Auto-detect input format
		std::string data_format;
		std::filesystem::path yaml_path = std::filesystem::path(project_path) / "input" / "dataProblem.yaml";
		std::filesystem::path dat_path = std::filesystem::path(project_path) / "input" / "dataProblem.dat";
		if (std::filesystem::exists(yaml_path))
		{
			data_format = "YAML";
			std::cout << "  Input format: YAML (detected dataProblem.yaml)" << std::endl;
		}
		else if (std::filesystem::exists(dat_path))
		{
			data_format = "ASCII";
			std::cout << "  Input format: ASCII (detected dataProblem.dat)" << std::endl;
		}
		else
		{
			std::cout << "ERROR: No input file found. Expected input/dataProblem.yaml or input/dataProblem.dat" << std::endl;
			return 1;
		}

		std::cout << std::endl
				  << "Creating simulation..." << std::endl;
		Simulation *mySim = new Simulation(project_path, data_format);
		std::cout << std::endl
				  << "Loading..." << std::endl;
		mySim->LoadCase();
		std::cout << std::endl
				  << "Initializing..." << std::endl;
		mySim->Initialize();
		std::cout << std::endl
				  << "Running..." << std::endl;
		mySim->Run();
		std::cout << std::endl
				  << "Closing simulation..." << std::endl;
		mySim->CloseCase();
	}
	catch (Exception &error)
	{
		error.PrintDebug();
	}

	std::cout << "\n";
	std::cout << "End of OASIS.\n";
	std::cout << "----------------------------------------------\n\n";

	return 0;
}