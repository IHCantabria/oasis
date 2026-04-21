// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef WINDTURBINE_FLAG
#define WINDTURBINE_FLAG

#include <armadillo>
#include <string>
#include <cstdio>
#include <iostream>
#include <yaml-cpp/yaml.h>

class Simulation;
class Body;

#ifdef OASIS_USE_OPENFAST
// ============================================================================
// Full WindTurbine class with OpenFAST/FASTurbine coupling
// ============================================================================
#include "FASTurbW_Library.h"
#include "FAST_Library.h"

class WindTurbine
{
public:
    int idWindTurbine; // Turbine number
    Simulation* pSim;  // Pointer to simulation instance
    Body* pBody;       // Pointer to corresponding body

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
    arma::mat bodyInerMat = arma::zeros(6, 6);
    arma::mat towrInerMat = arma::zeros(6, 6);
    arma::mat turbInerMat = arma::zeros(6, 6);

    double rotPos = 0.0;
    double rotSpeed = 0.0;
    double rotAcc = 0.0;
    int YCMode = 0;
    double yaw = 0.0;
    double yaw_ini = 0.0;
    double yawSpeed = 0.0;
    bool isRotorBlocked = false;

    arma::mat forceBodyCOG = arma::zeros(6, 1);
    double airTrq = 0.0;
    double genTrq = 0.0;
    double yawTrq = 0.0;

    WindTurbine(int n, Simulation* pSimInp)
    {
        idWindTurbine = n + 1;
        pSim = pSimInp;
    }
    void ReadPropertiesASCII(FILE* pFile);
    void ReadPropertiesYAML(YAML::Node node);
    void Initialize(void);
    void Finalize(void);
    void WriteOut(double t);

    void ComputeForces(double time);
    void SetInputsFAST(void);
    void ComputeControler(double time);
    void ComputeRotorAcc(void);

    void CheckError(void);
};

#else
// ============================================================================
// Stub WindTurbine class (OpenFAST not available)
// ============================================================================
class WindTurbine
{
public:
    int idWindTurbine;
    Simulation* pSim;
    Body* pBody;

    double rotIner = 0.0;
    arma::mat bodyInerMat = arma::zeros(6, 6);
    arma::mat towrInerMat = arma::zeros(6, 6);
    arma::mat turbInerMat = arma::zeros(6, 6);

    double rotPos = 0.0;
    double rotSpeed = 0.0;
    double rotAcc = 0.0;
    int YCMode = 0;
    double yaw = 0.0;
    double yaw_ini = 0.0;
    double yawSpeed = 0.0;
    bool isRotorBlocked = false;

    arma::mat forceBodyCOG = arma::zeros(6, 1);
    double airTrq = 0.0;
    double genTrq = 0.0;
    double yawTrq = 0.0;

    WindTurbine(int n, Simulation* pSimInp)
    {
        idWindTurbine = n + 1;
        pSim = pSimInp;
    }
    void ReadPropertiesASCII(FILE* pFile)
    {
        std::cerr << "ERROR: WindTurbine requires compilation with -DOASIS_USE_OPENFAST=ON" << std::endl;
    }
    void ReadPropertiesYAML(YAML::Node node)
    {
        std::cerr << "ERROR: WindTurbine requires compilation with -DOASIS_USE_OPENFAST=ON" << std::endl;
    }
    void Initialize(void)
    {
    }
    void Finalize(void)
    {
    }
    void WriteOut(double t)
    {
    }
    void ComputeForces(double time)
    {
    }
    void SetInputsFAST(void)
    {
    }
    void ComputeControler(double time)
    {
    }
    void ComputeRotorAcc(void)
    {
    }
    void CheckError(void)
    {
    }
};

#endif // OASIS_USE_OPENFAST

#endif // WINDTURBINE_FLAG