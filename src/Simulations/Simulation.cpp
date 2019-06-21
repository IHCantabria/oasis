
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include "Simulation.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Bodies/Bodies.hpp"


void Simulation::Initialize()
{
    // Read Simulation Properties
    this->ReadProperties();

    // Read Components Data
    this->ReadBodies();
    this->ReadLines();
    this->ReadBcps();
    
}


void Simulation::ReadProperties()
{
    (this->*pReadProperties)();
}


void Simulation::ReadBcps()
{
    (this->*pReadBcps)();
}


void Simulation::ReadBcpsASCII()
{
    // Declare local variables
    char bufferLine [1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosBCPs.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosBCPs.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
    
    // Read data
    fscanf(file_pointer, "%d %[^\n]", &numFairBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]", &numAnchorBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]", &numJointBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]", &numBodyBcps, bufferLine);
	numBcps = numFairBcps + numAnchorBcps + numJointBcps + numBodyBcps;

	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA BCP
	pBcps = new BCP*[numBcps];
	int bcp_count = 0;

    // Se leen los BCPs
	pFairleadBcps = new FairleadBCP*[numFairBcps];
	for(int ii=0; ii<numFairBcps; ii++)
    {
		pFairleadBcps[ii] = new FairleadBCP(bcp_count);
		pFairleadBcps[ii]->ReadPropertiesASCII(file_path);
		pFairleadBcps[ii]->GetValues(0.0);
		pBcps[bcp_count] = pFairleadBcps[ii];
		bcp_count++;
	}
	pAnchorBcps = new AnchorBCP*[numAnchorBcps];
	for(int ii=0; ii<numAnchorBcps; ii++)
    {
		pAnchorBcps[ii] = new AnchorBCP(bcp_count);
		pAnchorBcps[ii]->ReadPropertiesASCII(file_path);
		pAnchorBcps[ii]->GetValues(0.0);
		pAnchorBcps[bcp_count] = pAnchorBcps[ii];
		bcp_count++;
	}
	pJointBcps = new JointBCP*[numJointBcps];
	for(int ii=0; ii<numJointBcps; ii++)
    {
		pJointBcps[ii] = JointBCP(bcp_count);
		pJointBcps[ii]->ReadPropertiesASCII(file_path);
		pJointBcps[ii]->GetValues(0.0);
		pBcps[bcp_count] = pJointBcps[ii];
		bcp_count++;
	}
	pBodyBcps = new BodyBCP*[numBodyBcps];
	for(int ii=0; ii<numBodyBcps; ii++)
    {
		pBodyBcps[ii] = BodyBCP(bcp_count);
		pBodyBcps[ii].ReadPropertiesASCII(file_path);
		pBodyBcps[bcp_count] = pBodyBcps[ii];
		bcp_count++;
	}

    // Close file
    fclose(file_pointer);
}


void Simulation::ReadBcpsHDF5()
{
    std::stringstream ss;
    ss << "Method ReadBcpsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
}


void Simulation::ReadBodiesASCII()
{
    // Declare local variables
    char bufferLine [1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosBodies.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosBodies.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

    // Read total number of bodies to read 
    fscanf(file_pointer, "%d %[^\n]", &numBodies, bufferLine);

	//INICIO LOS CUERPOS
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii] = Body(ii);
		pBodies[ii].ReadPropertiesASCII(file_pointer);
		for(int jj=0; jj<pBodies[ii].numBcps; jj++)
        {
			pBodies[ii].pBodyBcps[jj] = BCPs[pBodies[ii].pIndexBcps[jj]];
		}
		pBodies[ii].UpdateBcps();
	}

	for(int ii=0; ii<nBCPs; ii=ii+1)
    {
		BCPs[ii]->getValues(0.0);
	}

    fclose(file_pointer);

}


void Simulation::ReadBodiesHDF5()
{
    std::stringstream ss;
    ss << "Method ReadBodiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
}


void Simulation::ReadPropertiesASCII()
{
    std::string file_path = JoinPath(inputFolderPath, "datosProblema.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosProblema.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
	
    char bufferLine [1000];
    int dummyBool;
    fscanf(file_pointer, "%lf %[^\n]", &gravity, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]", &waterDensity, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]", &waterDepth, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]", &maxTimeStep, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]", &simulationTime, bufferLine);
	fscanf(file_pointer, "%s %[^\n]", &timeIntMethod, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]", &timeIntAbsTol, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]", &timeIntRelTol, bufferLine);
	fscanf(file_pointer, "%d %[^\n]", &maxIterStep, bufferLine);
	fscanf(file_pointer, "%d %[^\n]", &dummyBool, bufferLine);
    readEquilibrium = dummyBool;
	fscanf(file_pointer, "%d %[^\n]", &dummyBool, bufferLine);
    writeEquilibrium = dummyBool;

    fclose(file_pointer);
}


void Simulation::ReadPropertiesHDF5()
{
    std::stringstream ss;
    ss << "Method ReadPropertiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
}


Simulation::Simulation(std::string incProjectPath, std::string incDataFormat)
{
    // Assign direct variables
    projectFolderPath = incProjectPath;

    // Process indirect variables
    if (!incDataFormat.compare("ASCII"))
    {
        dataFormat = 0;
        dataFormatStr = "ASCII";
        inputFolderPath = JoinPath(incProjectPath, "input");
        outputFolderPath = JoinPath(incProjectPath, "output");
        pReadProperties = &this->ReadPropertiesASCII;
        pReadBcps = &this->ReadBcpsASCII;
        pReadBodies = &this->ReadBodiesASCII;
    }
    else if (!incDataFormat.compare("HDF5"))
    {
        dataFormat = 1;
        dataFormatStr = "HDF5";
        inputFolderPath = incProjectPath;
        outputFolderPath = incProjectPath;
        pReadProperties = &this->ReadPropertiesHDF5;
        pReadBcps = &this->ReadBcpsHDF5;
        pReadBodies = &this->ReadBodiesHDF5;
    }
    else
    {
        std::stringstream ss;
        ss << "Simulation data format --> " << incDataFormat << " is not available.\n    Available formats: ASCII | HDF5.";
        throw ValueError(ss.str());
    }

    // 

}