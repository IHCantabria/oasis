
#include <armadillo>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include "WindTurbine.hpp"
#include "../MathTools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"
#include <OpenFAST.H>


void WindTurbine::ReadPropertiesASCII(FILE*& pFile){

    // Declare variables
	char buffer_line [1000];
    char cFastFileName [1000];
    int body_id;

    //Ignoro las tres primeras lineas, donde pone "New Turbine"
	for(int ii=0; ii<3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile);
	}

    //Leo todo
	fscanf(pFile, "%d %[^\n]\n", &body_id, buffer_line);
    pBody = pSim->pBodies[body_id-1];

    if (fscanf(pFile, "%s %[^\n]\n", cFastFileName, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the FAST input file name in Turbine: " << nTurbine << "\n";
		throw ValueError(ss.str());
	}

    std::string FASTInputFileName = JoinPath(pSim->inputFolderPath, "FAST");
    FASTInputFileName = JoinPath(FASTInputFileName, cFastFileName);


    fi.nTurbinesGlob = 1;
    fi.dryRun = false;
    fi.debug = false;
    fi.simStart = fast::init;
    fi.tStart = 0.0;
    fi.nEveryCheckPoint = 160;
    fi.dtFAST = pSim->maxTimeStep;
    fi.tMax = pSim->simulationTime;
    fi.scStatus = false;
    fi.scLibFile = "banana";
    fi.globTurbineData.resize(fi.nTurbinesGlob);

    fi.globTurbineData[0].TurbID = 1;
    fi.globTurbineData[0].FASTInputFileName = FASTInputFileName;
    fi.globTurbineData[0].FASTRestartFileName = "banana";
    fi.globTurbineData[0].TurbineBasePos = {0.0, 0.0, 0.0};
    fi.globTurbineData[0].TurbineHubPos = {0.0, 0.0, 0.0};
    fi.globTurbineData[0].numForcePtsBlade = 10;
    fi.globTurbineData[0].numForcePtsTwr = 10;
    fi.globTurbineData[0].nacelle_cd = 0.1;
    fi.globTurbineData[0].nacelle_area = 10;
    fi.globTurbineData[0].air_density = 1.225;


}


void WindTurbine::SetUp(void){

    FAST.setInputs(fi);
    FAST.allocateTurbinesToProcsSimple();
    FAST.init();
    if (FAST.isTimeZero()) FAST.solution0();

}


void WindTurbine::Finalize(void){

    FAST.end();

}

void WindTurbine::ComputeForceOnBase(void){

    forceOnBase = arma::zeros(6,1);
    // FAST.calc_nacelle_force(u,v,w,dc,area,rho,fx,fy,fz);
    // FAST.getForce(currentForce, iNode, iTurbGlob, nSize)

}

void WindTurbine::UpdateNodesVelocities(void){

}