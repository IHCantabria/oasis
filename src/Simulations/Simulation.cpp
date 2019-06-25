
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include "Simulation.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"


void Simulation::Initialize()
{
    // Read Simulation Properties
    this->ReadProperties();

    // Read Components Data
    this->ReadBodies();
    this->ReadLines();
    this->ReadBcps();
    this->ReadWinches();
    this->ReadSprings();

    // Setup case
    this->SetupCase();
    
}


void Simulation::ReadBcps()
{
    (this->*pReadBcps)();
}


void Simulation::ReadBcpsASCII()
{
    std::cout << "--> Reading BCPs Properties (ASCII format)" << std::endl;
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
    fscanf(file_pointer, "%d %[^\n]\n", &numFairBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numAnchorBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numJointBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numBodyBcps, bufferLine);
	numBcps = numFairBcps + numAnchorBcps + numJointBcps + numBodyBcps;
    printf("Number of bcps read!!\n");
	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA BCP
	pBcps = new BCP* [numBcps];
	int bcp_count = 0;
    printf("Allocated BCPs\n");
    // Se leen los BCPs
	pFairleadBcps = new FairleadBCP* [numFairBcps];
	for(int ii=0; ii<numFairBcps; ii++)
    {
		pFairleadBcps[ii] = new FairleadBCP(bcp_count);
		pFairleadBcps[ii]->ReadPropertiesASCII(file_pointer);
		dynamic_cast<FairleadBCP*>(pFairleadBcps[ii])->Initialize(inputFolderPath);
		pBcps[bcp_count] = pFairleadBcps[ii];
		bcp_count++;
	}
	pAnchorBcps = new AnchorBCP* [numAnchorBcps];
	for(int ii=0; ii<numAnchorBcps; ii++)
    {
		pAnchorBcps[ii] = new AnchorBCP(bcp_count);
		pAnchorBcps[ii]->ReadPropertiesASCII(file_pointer);
		pBcps[bcp_count] = pAnchorBcps[ii];
		bcp_count++;
	}
	pJointBcps = new JointBCP* [numJointBcps];
	for(int ii=0; ii<numJointBcps; ii++)
    {
		pJointBcps[ii] = new JointBCP(bcp_count);
		pJointBcps[ii]->ReadPropertiesASCII(file_pointer);
		pBcps[bcp_count] = pJointBcps[ii];
		bcp_count++;
	}
	pBodyBcps = new BodyBCP* [numBodyBcps];
	for(int ii=0; ii<numBodyBcps; ii++)
    {
		pBodyBcps[ii] = new BodyBCP(bcp_count);
		pBodyBcps[ii]->ReadPropertiesASCII(file_pointer);
		pBcps[bcp_count] = pBodyBcps[ii];
		bcp_count++;
	}

    // Close file
    fclose(file_pointer);
    std::cout << "----> BCPs Properties Read" << std::endl;
}


void Simulation::ReadBcpsHDF5()
{
    std::cout << "--> Reading BCPs Properties (ASCII format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBcpsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> BCPs Properties Read" << std::endl;
}


void Simulation::ReadBodies()
{
    (this->*pReadBodies)();
}


void Simulation::ReadBodiesASCII()
{
    std::cout << "--> Reading Bodies Properties (ASCII format)" << std::endl;
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
    fscanf(file_pointer, "%d %[^\n]\n", &numBodies, bufferLine);

	//INICIO LOS CUERPOS
    pBodies = new Body* [numBodies];
    printf("Total number of bodies: %d\n", numBodies);
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii] = new Body(ii);
		pBodies[ii]->ReadPropertiesASCII(file_pointer);
        /**
		for(int jj=0; jj<pBodies[ii]->numBcps; jj++)
        {
			pBodies[ii]->pBodyBcps[jj] = BCPs[pBodies[ii]->pIndexBcps[jj]];
		}
		pBodies[ii]->UpdateBcps();
        **/
	}
    /**
	for(int ii=0; ii<nBCPs; ii=ii+1)
    {
		BCPs[ii]->getValues(0.0);
	}
    **/
    fclose(file_pointer);
    std::cout << "----> Bodies Properties Read" << std::endl;
}


void Simulation::ReadBodiesHDF5()
{
    std::cout << "--> Reading Bodies Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBodiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Bodies Properties Read" << std::endl;
}


void Simulation::ReadLines()
{
    (this->*pReadLines)();
}


void Simulation::ReadLinesASCII()
{
    std::cout << "--> Reading Lines Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
    FILE* fid_lines;
	std::string file_path = JoinPath(inputFolderPath, "datosLines.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosLines.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

   // Read total number of springs to read 
    fscanf(file_pointer, "%d %[^\n]\n", &numLines, bufferLine);

	pLines = new Line*[numLines];
	//INICIO LAS LINEAS
	for(int ii=0; ii<numLines; ii++)
    {
		pLines[ii] = new Line(ii);
		try
		{
            pLines[ii]->ReadPropertiesASCII(file_pointer);
            pLines[ii]->print_out();
            /**
            pLines[ii]->LineBCP[0] = BCPs[pLines[ii]->BCP_1-1];
            pLines[ii]->LineBCP[1] = BCPs[pLines[ii]->BCP_N-1];
            pLines[ii]->pos_1 = BCPs[pLines[ii].BCP_1-1]->pos;
            pLines[ii]->pos_N = BCPs[pLines[ii].BCP_N-1]->pos;
            nNodosTotal=nNodosTotal+pLines[ii].N;
            if(flag_read_eq==0) pLines[ii]->initLine();
            Lines[ii].SEM_getBaseFunctions();
            **/
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is under the floor level." << std::endl << std::endl;
			if (e==1) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor and it shouldn't. " << std::endl << std::endl;
			if (e==2) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl << std::endl;
			if (e==3) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and vertical. It should be pretensed. " << std::endl << std::endl;
			if (e==4) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and it should. " << std::endl << std::endl;
			if (e==5) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " initial shape can't be computed with QS method. " << std::endl;
			if (e==6) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor althoug none of its ends are there. " << std::endl << std::endl;
		}
        
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Lines Properties Read" << std::endl;
}


void Simulation::ReadLinesHDF5()
{
    std::cout << "--> Reading Lines Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadLinesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Lines Properties Read" << std::endl;
}


void Simulation::ReadSprings()
{
    (this->*pReadSprings)();
}


void Simulation::ReadSpringsASCII()
{
    std::cout << "--> Reading Spring Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
    FILE* fid_lines;
	std::string file_path = JoinPath(inputFolderPath, "datosSprings.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosSprings.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

    // Read file contents
    fscanf(file_pointer, "%d %[^\n]", &numSprings, bufferLine);

	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA MUELLE
	pSprings = new Spring*[numSprings];

	//INICIO LLOS MUELLES
	for(int ii=0; ii<numSprings; ii++)
    {
		pSprings[ii] = new Spring(ii);
		pSprings[ii]->ReadPropertiesASCII(file_path);
		//pSprings[ii].SpringBCP[0] = BCPs[Springs[ii].BCP_1-1];
		//pSprings[ii].SpringBCP[1] = BCPs[Springs[ii].BCP_2-1];
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Spring Properties Read" << std::endl;
}


void Simulation::ReadSpringsHDF5()
{
    std::cout << "--> Reading Spring Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadSpringsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Spring Properties Read" << std::endl;
}


void Simulation::ReadProperties()
{
    (this->*pReadProperties)();
}


void Simulation::ReadPropertiesASCII()
{
    std::cout << "--> Reading Simulation Properties (ASCII format)" << std::endl;
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
    printf("MaxIterStep: %d\n", maxIterStep);

    // Close file
    fclose(file_pointer);

    std::cout << "----> Simulation Properties Read" << std::endl;
}


void Simulation::ReadPropertiesHDF5()
{
    std::cout << "--> Reading Simulation Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadPropertiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Simulation Properties Read" << std::endl;
}


void Simulation::ReadWinches()
{
    (this->*pReadWinches)();
}


void Simulation::ReadWinchesASCII()
{
    std::cout << "--> Reading Winches Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
    FILE* fid_lines;
	std::string file_path = JoinPath(inputFolderPath, "datosWinchies.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosWinchies.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

    // Read file contents
	fscanf(file_pointer, "%d %[^\n]\n", &numWinches, bufferLine);

	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA MUELLE
	pWinches = new Winchie*[numWinches];
	//INICIO LLOS MUELLES
	for(int ii=0; ii<numWinches; ii++)
    {
		pWinches[ii] = new Winchie(ii);
		pWinches[ii]->ReadPropertiesASCII(file_pointer);
		pWinches[ii]->LineW = pLines[pWinches[ii]->nLine - 1];
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Winches Properties Read" << std::endl;
}


void Simulation::ReadWinchesHDF5()
{
    std::cout << "--> Reading Winches Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWinchesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Winches Properties Read" << std::endl;
}


void Simulation::SetupCase()
{
    // Setup Bodies
    double a = 0.0;
    /**
    // Setup Lines
    for (int ii=0; ii<numLines; ii++)
    {

    }

    // Setup Springs
	for(int ii=0; ii<numSprings; ii++)
    {
		pSprings[ii].SpringBCP[0] = BCPs[Springs[ii].BCP_1-1];
		pSprings[ii].SpringBCP[1] = BCPs[Springs[ii].BCP_2-1];
	}
    **/
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
        pReadLines = &this->ReadLinesASCII;
        pReadSprings = &this->ReadSpringsASCII;
        pReadWinches = &this->ReadWinchesASCII;
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
        pReadLines = &this->ReadLinesHDF5;
        pReadSprings = &this->ReadSpringsHDF5;
        pReadWinches = &this->ReadWinchesHDF5;
    }
    else
    {
        std::stringstream ss;
        ss << "Simulation data format --> " << incDataFormat << " is not available.\n    Available formats: ASCII | HDF5.";
        throw ValueError(ss.str());
    }

    // 

}