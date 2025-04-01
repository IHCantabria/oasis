
#ifndef simulationdef_hpp__
#define simulationdef_hpp__

#include <string>
#include "ISimulation.hpp"
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../BCPs/Winchies.hpp"
#include "../BCPs/WinchiesController.hpp"
#include "../Lines/Lines.hpp"
#include "../Spring/Spring.hpp"
#include "../Waves/Wave.hpp"
#include "../Sinking/Sinking.hpp"
#include "../SeaFloor/SeaFloor.hpp"
#include "../WindTurbine/WindTurbine.hpp"

class ODE_solver;

class Simulation : public ISimulation
{
private:
    int dataFormat;
    std::string dataFormatStr;
    bool status;

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
    double writeTimeStep;
    double maxTimeStep;
    double hydroTimeStep;
    double fastTimeStep;
    double fastControllerTimeStep;
    double timeIRF;
    double sinkingTimeStep;
    double controllerTimeStep;
    int numSystem;
    int numSystem2;
    double simulationTime;
    bool readEquilibrium;
    bool rotSimpFlag;
    double timeIntAbsTol;
    int timeIntMethod;
    int timeIntOrder;
    bool timeIntAdaptivity;
    double timeIntRelTol;
    double waterDensity;
    double waterDepth;
    bool writeEquilibrium;
    int flagStatic;

    // Declare time simulation attributes
    int numCallsSysFun = 0;
    ODE_solver *pTimeSolver;
    int timeBufferSize = 0;
    arma::mat timeBuffer = arma::zeros(1, timeBufferSize);
    int timeBufferCount = 0;

    // Declare system properties
    arma::mat *pSystemMatrix;
    arma::mat *pSystemMatrixInv;
    arma::mat *pSystemMatrixFF;
    arma::mat *pSystemMatrixFFInv;
    arma::mat *pSystemMatrixFL;
    arma::uvec sysMatIndFree, sysMatIndLock;

    int numAllLinesNodes = 0;
    arma::sp_mat *pLinesCouplingMatrix_sp;
    arma::mat *pLinesCouplingMatrix;
    arma::mat *pLinesCouplingMatrixInv;
    arma::uvec indexesNoFairNoAnchor;
    arma::uvec indexesFairAnchor;
    arma::mat removed_rows;

    // Declare Components Setup Attributes
    AnchorBCP **pAnchorBcps;
    BCP **pBcps;
    Body **pBodies;
    Body **pBodiesFree;
    Body **pBodiesLock;
    BodyBCP **pBodyBcps;
    FairleadBCP **pFairleadBcps;
    JointBCP **pJointBcps;
    Line **pLines;
    Spring **pSprings;
    Winchie **pWinches;
    WinchieController WinchesController;
    Wave *pWave;
    Sinking **pSinking;
    SeaFloor **pSeaFloor;
    Bathymetry **pBathymetry;
    Inclined **pInclined;
    Flat **pFlat;
    WindTurbine **pWindTurbines;
    int numAnchorBcps = 0;
    int numBcps = 0;
    int numBodies = 0;
    int numBodiesFree = 0;
    int numBodiesLock = 0;
    int numBodyBcps = 0;
    int numDofTotal = 0;
    int numFairBcps = 0;
    int numJointBcps = 0;
    int numLines = 0;
    int numSprings = 0;
    int numWinches = 0;
    int numSinking = 0;
    int numFloor = 0;
    int numBathymetry = 0;
    int numInclined = 0;
    int numFlat = 0;
    int numWindTurbines = 0;

    // Declare constructors
    Simulation(std::string projectPath, std::string incDataFormat);

    // Declare IO methods
    void LoadCase(void);
    void (Simulation::*pReadBcps)(void);
    void (Simulation::*pReadBodies)(void);
    void (Simulation::*pReadLines)(void);
    void (Simulation::*pReadSprings)(void);
    void (Simulation::*pReadProperties)(void);
    void (Simulation::*pReadWaves)(void);
    void (Simulation::*pReadWind)(void);
    void (Simulation::*pReadWinches)(void);
    void (Simulation::*pReadSinking)(void);
    void (Simulation::*pReadSeaFloor)(void);
    void (Simulation::*pReadWindTurbines)(void);
    void PrintSetup(void);
    void ReadBcps(void);
    void ReadBcpsASCII(void);
    void ReadBcpsHDF5(void);
    void ReadBodies(void);
    void ReadBodiesASCII(void);
    void ReadBodiesHDF5(void);
    void ReadLines(void);
    void ReadLinesASCII(void);
    void ReadLinesHDF5(void);
    void ReadSinking(void);
    void ReadSinkingASCII(void);
    void ReadSinkingHDF5(void);
    void ReadSprings(void);
    void ReadSpringsASCII(void);
    void ReadSpringsHDF5(void);
    void ReadWinches(void);
    void ReadWinchesASCII(void);
    void ReadWinchesHDF5(void);
    void ReadProperties(void);
    void ReadPropertiesASCII(void);
    void ReadPropertiesHDF5(void);
    void ReadWaves(void);
    void ReadWavesASCII(void);
    void ReadWavesHDF5(void);
    void ReadSeaFloor(void);
    void ReadSeaFloorASCII(void);
    void ReadSeaFloorHDF5(void);
    void ReadWindTurbines(void);
    void ReadWindTurbinesASCII(void);
    void ReadWindTurbinesHDF5(void);

    // Declare general purpose class methods
    arma::mat CalculateSystemDynamics(double time, arma::mat y);
    void Initialize(void);
    void Run();
    void SetupCase(void);
    void CloseCase(void);
    void UpdateSystem(void);
    void UpdateSystemMatrix(void);
    void ComputeLinesCouplingMatrix(void);
    arma::mat ComputeLinesInitialPoint(void);
    arma::mat ComputeLinesForces(arma::mat position);
    arma::mat ComputeLinesJacobian(arma::mat position, arma::mat fuerzaEnPosInicial);
    void ComputeLinesEquilibrium(arma::mat x);
};

#endif // simulation_hpp__
