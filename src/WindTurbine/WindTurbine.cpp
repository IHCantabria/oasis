
#include <armadillo>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include "WindTurbine.hpp"
#include "../MathTools.hpp"
#include "../Exceptions/Exception.hpp"
#include <OpenFAST.H>


void WindTurbine::ReadPropertiesASCII(std::string file_path){

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

    if (fscanf(pFile, "%s %[^\n]\n", cHydroDatabaseName, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the FAST input file name in Turbine: " << nTurbine << "\n";
		throw ValueError(ss.str());
	}


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
    fi.globTurbineData[0].FASTInputFileName = cHydroDatabaseName;
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