
#ifndef simulation_hpp__
#define simulation_hpp__

#include <string>
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"


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

    // Declare Components Setup Attributes
    BCP** pAnchorBcps;
    BCP** pFairleadBcps;
    BCP** pJointBcps;
    BCP** pBodyBcps;
    BCP** pBcps;
    Body** pBodies;
    int numAnchorBcps;
    int numBcps;
    int numBodies;
    int numBodyBcps;
    int numFairBcps;
    int numJointBcps;

    // Declare constructors
    Simulation(std::string projectPath, std::string incDataFormat);

    // Declare IO methods
    void Initialize(void);
    void (Simulation::*pReadProperties)(void);
    void (Simulation::*pReadBcps)(void);
    void (Simulation::*pReadBodies)(void);
    void ReadBcps(void);
    void ReadBcpsASCII(void);
    void ReadBcpsHDF5(void);
    void ReadBodies(void);
    void ReadBodiesASCII(void);
    void ReadBodiesHDF5(void);
    void ReadLines(void);
    void ReadLinesASCII(void);
    void ReadLinesHDF5(void);
    void ReadProperties(void);
    void ReadPropertiesASCII(void);
    void ReadPropertiesHDF5(void);
    

};

#endif // simulation_hpp__