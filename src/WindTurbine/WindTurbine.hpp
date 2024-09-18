
#ifndef WINDTURBINE_FLAG
#define WINDTURBINE_FLAG

#include <armadillo>
#include <string>
#include "FASTurbW_Library.h"
#include "FAST_Library.h"

class Simulation;
class Body;

class WindTurbine
{
public:
    int idWindTurbine; // Turbine number
    Simulation *pSim;  // Pointer to simulation instance
    Body *pBody;       // Pointer to corresponding body

    char InputFileName_AD[INTERFACE_STRING_LENGTH];
    char InputFileName_IW[INTERFACE_STRING_LENGTH];
    char InputFileName_SD[INTERFACE_STRING_LENGTH];
    char InputFileName_ED[INTERFACE_STRING_LENGTH];
    char OutputPathName[INTERFACE_STRING_LENGTH];
    FSTW_InitInputType_t FSTW_InitInput;
    FSTW_InputType_t FSTW_Input;
    FSTW_OutputType_t FSTW_Output;
    int ErrStat;
    char ErrMsg[INTERFACE_STRING_LENGTH];
    int AbortErrLev = ErrID_Fatal;

    double rotIner = 0.0;
    arma::mat bodyInerMat = arma::zeros(6, 6); // Inertia matrix of the platform without WT in the platform COG frame
    arma::mat towrInerMat = arma::zeros(6, 6); // Inertia matrix of the current tower without WT in the platform COG frame
    arma::mat turbInerMat = arma::zeros(6, 6); // Inertia matrix of the current turbine (nacelle+hub+blades) without WT in the platform COG frame

    double rotPos = 0.0;
    double rotSpeed = 0.0;
    double rotAcc = 0.0;
    int YCMode = 0;
    double yaw = 0.0;
    double yaw_ini = 0.0;
    double yawSpeed = 0.0;

    arma::mat forceBodyCOG = arma::zeros(6, 1);
    double airTrq = 0.0;
    double genTrq = 0.0;
    double yawTrq = 0.0;

    WindTurbine(int n, Simulation *pSimInp)
    {
        idWindTurbine = n + 1;
        pSim = pSimInp;
    }                                      // Inicializa un objeto de clase turbina dandole el indice
    void ReadPropertiesASCII(FILE *pFile); // Lee inputs de las turbinas
    void Initialize(void);                 // Configura el objeto turbina
    void Finalize(void);                   // Cierra el caso
    void WriteOut(double t);               // Escribir datos a fichero

    void ComputeForces(double time);
    void SetInputsFAST(void);
    void ComputeControler(double time);
    void ComputeRotorAcc(void);

    void CheckError(void);
};

#endif