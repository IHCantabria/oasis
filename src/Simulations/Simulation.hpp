
#include <string>


class Simulation
{
private:
    int dataFormat;
    std::string dataFormatStr;
public:
    // Declare file system attributes
    std::string inputFolderPath;
    std::string projectFolderPath;
    std::string outputFolderPath;

    // Declare simulation properties attributes
    double gravity;
    int maxIterStep;
    double maxTimeStep;
    double simulationTime;
    bool readEquilibrium;
    double timeIntAbsTol;
    std::string timeIntMethod;
    double timeIntRelTol;
    double waterDensity;
    double waterDepth;
    bool writeEquilibrium;

    // Declare class functions
    void (Simulation::*pReadProperties)(void);
    void ReadProperties(void);
    void ReadPropertiesASCII(void);
    void ReadPropertiesHDF5(void);
    Simulation(std::string projectPath, std::string incDataFormat);
    
};