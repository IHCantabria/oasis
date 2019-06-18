
#include <cstdio>
#include <string>
#include <sstream>
#include "Simulation.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"


void Simulation::ReadProperties()
{
    (this->*pReadProperties)();
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
    double a = 0.0;
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
    }
    else if (!incDataFormat.compare("HDF5"))
    {
        dataFormat = 1;
        dataFormatStr = "HDF5";
        inputFolderPath = incProjectPath;
        outputFolderPath = incProjectPath;
        pReadProperties = &this->ReadPropertiesHDF5;
    }
    else
    {
        std::stringstream ss;
        ss << "Simulation data format --> " << incDataFormat << " is not available.\n    Available formats: ASCII | HDF5.";
        throw ValueError(ss.str());
    }

    // 

}