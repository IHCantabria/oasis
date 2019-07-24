
#ifndef simulationdef_hpp__
#define simulationdef_hpp__

#include <string>
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../BCPs/Winchies.hpp"
#include "../Lines/Lines.hpp"
#include "../Spring/Spring.hpp"

class BDF;


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
    int numSystem;
    int numSystem2;
    double simulationTime;
    bool readEquilibrium;
    double timeIntAbsTol;
    int timeIntMethod;
    double timeIntRelTol;
    double waterDensity;
    double waterDepth;
    bool writeEquilibrium;

    // Declare time simulation attributes
    int numCallsSysFun=0;
    BDF* pTimeSolver;
    int timeBufferSize=1e2;
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
    int numAnchorBcps=0;
    int numBcps=0;
    int numBodies=0;
    int numBodyBcps=0;
    int numDofTotal=0;
    int numFairBcps=0;
    int numJointBcps=0;
    int numLines=0;
    int numSprings=0;
    int numWinches=0;

    // Declare constructors
    Simulation(std::string projectPath, std::string incDataFormat);

    // Declare IO methods
    void LoadCase(void);
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
    arma::mat CalculateSystemDynamics(double time, arma::mat y);
    void Initialize(void);
    void Run();
    void SetupCase(void);
    void UpdateSystem(void);

};

#endif // simulation_hpp__