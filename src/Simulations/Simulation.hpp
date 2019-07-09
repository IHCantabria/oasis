
#ifndef simulationdef_hpp__
#define simulationdef_hpp__

#include <string>
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../BCPs/Winchies.hpp"
#include "../Lines/Lines.hpp"
#include "../Spring/Spring.hpp"


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

    // Declare IO attributes
    bool useWinches = false;

    // Declare simulation properties attributes
    double gravity;
    int maxIterStep;
    double maxTimeStep;
    double simulationTime;
    bool readEquilibrium;
    double timeIntAbsTol;
    int timeIntMethod;
    double timeIntRelTol;
    double waterDensity;
    double waterDepth;
    bool writeEquilibrium;

    // Declare time simulation attributes
    int timeBufferSize=1e4;
	arma::mat timeBuffer = arma::zeros(1, timeBufferSize);
    int timeBufferCount=0;

    // Declare Components Setup Attributes
    AnchorBCP** pAnchorBcps;
    BCP** pBcps;
    Body** pBodies;
    BodyBCP** pBodyBcps;
    FairleadBCP** pFairleadBcps;
    JointBCP** pJointBcps;
    Line** pLines;
    Spring** pSprings;
    Winchie** pWinches;
    int numAnchorBcps;
    int numBcps;
    int numBodies;
    int numBodyBcps;
    int numDofTotal;
    int numFairBcps;
    int numJointBcps;
    int numLines;
    int numSprings;
    int numWinches;

    // Declare constructors
    Simulation(std::string projectPath, std::string incDataFormat);

    // Declare IO methods
    void Initialize(void);
    void (Simulation::*pReadBcps)(void);
    void (Simulation::*pReadBodies)(void);
    void (Simulation::*pReadLines)(void);
    void (Simulation::*pReadSprings)(void);
    void (Simulation::*pReadProperties)(void);
    void (Simulation::*pReadWinches)(void);
    void ReadBcps(void);
    void ReadBcpsASCII(void);
    void ReadBcpsHDF5(void);
    void ReadBodies(void);
    void ReadBodiesASCII(void);
    void ReadBodiesHDF5(void);
    void ReadHydrodynamicsHDF5(void);
    void ReadLines(void);
    void ReadLinesASCII(void);
    void ReadLinesHDF5(void);
    void ReadSprings(void);
    void ReadSpringsASCII(void);
    void ReadSpringsHDF5(void);
    void ReadWinches(void);
    void ReadWinchesASCII(void);
    void ReadWinchesHDF5(void);
    void ReadProperties(void);
    void ReadPropertiesASCII(void);
    void ReadPropertiesHDF5(void);

    // Declare general purpose class methods
    void SetupCase();
    void UpdateSystem(arma::mat y);

};

#endif // simulation_hpp__