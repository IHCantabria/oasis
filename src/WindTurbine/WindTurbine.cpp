
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
#include <mpi.h>


void WindTurbine::ReadPropertiesASCII(FILE*& pFile){

    // Declare variables
	char buffer_line [1000];
    char cFastFileName [1000];

    fi.nTurbinesGlob = numWindTurbines;
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

    pBodies = new Body* [numWindTurbines];

    for(int iTurbLoc=0; iTurbLoc<numWindTurbines; iTurbLoc++)
    {
        // Ignoro las tres primeras lineas, donde pone "New Turbine"
        for(int ii=0; ii<3; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }

        // Leo y asigno el cuerpo correspondiente a la turbina
        int tmp_body_id;
        fscanf(pFile, "%d %[^\n]\n", &tmp_body_id, buffer_line); tmp_body_id--;
        pBodies[iTurbLoc] = pSim->pBodies[tmp_body_id];

        // Leo y uno al path el nombre del archivo de input de fast
        if (fscanf(pFile, "%s %[^\n]\n", cFastFileName, buffer_line) != 2)
        {
            std::stringstream ss;
            ss << "An error ocurred when trying to read the FAST input file name in Turbine: " << iTurbLoc+1 << "\n";
            throw ValueError(ss.str());
        }
        std::string FASTInputFileName = JoinPath(pSim->inputFolderPath, "FAST");
        FASTInputFileName = JoinPath(FASTInputFileName, cFastFileName);

        // Leo el numero de puntos en palas y torre        
        int tmp_numForcePtsBlade,tmp_numForcePtsTwr;
        fscanf(pFile, "%d %[^\n]\n", &tmp_numForcePtsBlade, buffer_line);
        fscanf(pFile, "%d %[^\n]\n", &tmp_numForcePtsTwr, buffer_line);

        // Leo la posicion de la base de la turbina
        double TurbineBasePosX,TurbineBasePosY,TurbineBasePosZ;
        fscanf(pFile, "%lf %lf %lf %[^\n]\n", &TurbineBasePosX, &TurbineBasePosY, &TurbineBasePosZ, buffer_line);

        // Guardo en estructura de fast los datos leidos
        fi.globTurbineData[iTurbLoc].TurbID = iTurbLoc+1;
        fi.globTurbineData[iTurbLoc].FASTInputFileName = FASTInputFileName;
        fi.globTurbineData[iTurbLoc].FASTRestartFileName = "banana";
        fi.globTurbineData[iTurbLoc].TurbineBasePos = {TurbineBasePosX,TurbineBasePosY,TurbineBasePosZ};
        // fi.globTurbineData[iTurbLoc].TurbineHubPos = {0.0, 0.0, 0.0};
        fi.globTurbineData[iTurbLoc].numForcePtsBlade = tmp_numForcePtsBlade;
        fi.globTurbineData[iTurbLoc].numForcePtsTwr = tmp_numForcePtsTwr;
        // fi.globTurbineData[iTurbLoc].nacelle_cd = 0.1;
        // fi.globTurbineData[iTurbLoc].nacelle_area = 10;
        // fi.globTurbineData[iTurbLoc].air_density = 1.225;
    }

}


void WindTurbine::Initialize(void){

    int iErr;
    int nProcs;
    int rank;
    std::cout << "        --> MPI_Init" << std::endl;
    iErr = MPI_Init(NULL, NULL);
    std::cout << "        --> MPI_Comm_size" << std::endl;
    iErr = MPI_Comm_size( MPI_COMM_WORLD, &nProcs);
    std::cout << "        --> MPI_Comm_rank" << std::endl;
    iErr = MPI_Comm_rank( MPI_COMM_WORLD, &rank);

    std::cout << "        --> fi.comm = MPI_COMM_WORLD;" << std::endl;
    fi.comm = MPI_COMM_WORLD;

    std::cout << "        --> FAST.setInputs(fi);" << std::endl;
    FAST.setInputs(fi);
    FAST.allocateTurbinesToProcsSimple();
    FAST.init();
    if (FAST.isTimeZero()) FAST.solution0();

}


void WindTurbine::Step(void){
    SetBaseMovements();
    FAST.step();
    GetBaseForces();
}


void WindTurbine::Finalize(void){
    FAST.end();
    MPI_Finalize();

    char* s1 = new char[1000]; char* s2 = new char[1000]; char* s3 = new char[1000];
    #if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
        sprintf(s1,"move %s\FAST\*.out %s",(pSim->inputFolderPath).c_str(),(pSim->outputFolderPath).c_str());
        sprintf(s2,"move %s\FAST\*.outb %s",(pSim->inputFolderPath).c_str(),(pSim->outputFolderPath).c_str());
        sprintf(s3,"move %s\FAST\*.sum %s",(pSim->inputFolderPath).c_str(),(pSim->outputFolderPath).c_str());
        system("del *.h5");
    #elif defined(__unix__)
        sprintf(s1,"mv %s/FAST/*.out %s",(pSim->inputFolderPath).c_str(),(pSim->outputFolderPath).c_str());
        sprintf(s2,"mv %s/FAST/*.outb %s",(pSim->inputFolderPath).c_str(),(pSim->outputFolderPath).c_str());
        sprintf(s3,"mv %s/FAST/*.sum %s",(pSim->inputFolderPath).c_str(),(pSim->outputFolderPath).c_str());
        system("rm *.h5");
    #endif
    system(s1); system(s2); system(s3);
}

void WindTurbine::GetBaseForces(void){

}

void WindTurbine::SetBaseMovements(void){

}