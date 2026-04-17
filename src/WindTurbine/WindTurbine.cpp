
// WindTurbine.cpp - Only compiled when OASIS_USE_OPENFAST is enabled
#ifdef OASIS_USE_OPENFAST

#include <armadillo>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include "WindTurbine.hpp"
#include "../MathTools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"
#include "FASTurbW_Library.h"
#include "FAST_Library.h"

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#if (__cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)) && __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#elif __has_include(<boost/filesystem.hpp>)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif
#endif

void WindTurbine::ReadPropertiesASCII(FILE *pFile)
{

    // Declare variables
    char buffer_line[1000];
    char cADFileName[1000];
    char cIWFileName[1000];
    char cSDFileName[1000];
    char cEDFileName[1000];

    // Ignoro las tres primeras lineas, donde pone "New Turbine"
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), pFile);
    }

    // Leo y asigno el cuerpo correspondiente a la turbina
    int tmp_body_id;
    fscanf(pFile, "%d %[^\n]\n", &tmp_body_id, buffer_line);
    tmp_body_id--;
    pBody = pSim->pBodies[tmp_body_id];

    // Leo y uno al path el nombre del archivo de input de AeroDyn
    if (fscanf(pFile, "%s %[^\n]\n", cADFileName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the AeroDyn input file name in Turbine: " << idWindTurbine << "\n";
        throw ValueError(ss.str());
    }
    std::string ADFileName = JoinPath(pSim->inputFolderPath, cADFileName);
    std::copy(ADFileName.data(),
              ADFileName.data() + (ADFileName.size() + 1),
              InputFileName_AD);

    // Leo y uno al path el nombre del archivo de input de InflowWind
    if (fscanf(pFile, "%s %[^\n]\n", cIWFileName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the InflowWind input file name in Turbine: " << idWindTurbine << "\n";
        throw ValueError(ss.str());
    }
    std::string IWFileName = JoinPath(pSim->inputFolderPath, cIWFileName);
    std::copy(IWFileName.data(),
              IWFileName.data() + (IWFileName.size() + 1),
              InputFileName_IW);

    // Leo y uno al path el nombre del archivo de input de ServoDyn
    if (fscanf(pFile, "%s %[^\n]\n", cSDFileName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the ServoDyn input file name in Turbine: " << idWindTurbine << "\n";
        throw ValueError(ss.str());
    }
    std::string SDFileName = JoinPath(pSim->inputFolderPath, cSDFileName);
    std::copy(SDFileName.data(),
              SDFileName.data() + (SDFileName.size() + 1),
              InputFileName_SD);

    // Leo y uno al path el nombre del archivo de input de ElastoDyn
    if (fscanf(pFile, "%s %[^\n]\n", cEDFileName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the ElastoDyn input file name in Turbine: " << idWindTurbine << "\n";
        throw ValueError(ss.str());
    }
    std::string EDFileName = JoinPath(pSim->inputFolderPath, cEDFileName);
    std::copy(EDFileName.data(),
              EDFileName.data() + (EDFileName.size() + 1),
              InputFileName_ED);
}

void WindTurbine::ReadPropertiesYAML(YAML::Node node)
{
    int tmp_body_id = node["body_index"].as<int>() - 1;
    pBody = pSim->pBodies[tmp_body_id];

    std::string ADFileName = JoinPath(pSim->inputFolderPath, node["aerodyn_file"].as<std::string>());
    std::copy(ADFileName.data(), ADFileName.data() + (ADFileName.size() + 1), InputFileName_AD);

    std::string IWFileName = JoinPath(pSim->inputFolderPath, node["inflowwind_file"].as<std::string>());
    std::copy(IWFileName.data(), IWFileName.data() + (IWFileName.size() + 1), InputFileName_IW);

    std::string SDFileName = JoinPath(pSim->inputFolderPath, node["servodyn_file"].as<std::string>());
    std::copy(SDFileName.data(), SDFileName.data() + (SDFileName.size() + 1), InputFileName_SD);

    std::string EDFileName = JoinPath(pSim->inputFolderPath, node["elastodyn_file"].as<std::string>());
    std::copy(EDFileName.data(), EDFileName.data() + (EDFileName.size() + 1), InputFileName_ED);
}

void WindTurbine::Initialize(void)
{

    std::cout << "    WindTurbine::Initialize" << std::endl;
    std::string OutputPath = JoinPath(pSim->outputFolderPath, "FAST");
    std::copy(OutputPath.data(),
              OutputPath.data() + (OutputPath.size() + 1),
              OutputPathName);
    fs::create_directory(OutputPath);

    std::cout << "        --> Introduce inputs in FSTW_InitInput" << std::endl;
    FSTW_InitInput.Tmax = pSim->simulationTime;
    FSTW_InitInput.TimeInterval = pSim->fastTimeStep;
    FSTW_InitInput.TimeInterval_SrvD = pSim->fastControllerTimeStep;
    // Set motion position input [1: Platform, 2: Hub]
    // Only platform is implemented
    // TODO: Implement hub motion option
    FSTW_InitInput.flag_input_pos = 1;

    std::cout << "        --> FSTW_Init..." << std::endl;
    FSTW_Init(&idWindTurbine,
              InputFileName_AD,
              InputFileName_IW,
              InputFileName_SD,
              InputFileName_ED,
              OutputPathName,
              &FSTW_InitInput,
              &FSTW_Input,
              &FSTW_Output,
              &ErrStat,
              ErrMsg);
    CheckError();
    std::cout << "        --> ... done!" << std::endl;

    std::cout << "        --> Compute inertia matrices" << std::endl;
    rotIner = FSTW_InitInput.rotorInertia;
    for (int ii = 0; ii < 6; ii++)
    {
        for (int jj = 0; jj < 6; jj++)
        {
            this->bodyInerMat(ii, jj) = FSTW_InitInput.platInerMat[ii][jj];
            this->towrInerMat(ii, jj) = FSTW_InitInput.towrInerMat[ii][jj];
            this->turbInerMat(ii, jj) = FSTW_InitInput.turbInerMat[ii][jj];
        }
    }
    std::cout << "        --> ... done!" << std::endl;

    std::cout << "        --> Save inertia matrices" << std::endl;
    char buffer1[50], buffer2[50], buffer3[50];
    int nn1 = sprintf(buffer1, "FASTurbW_%d_bodyInerMat.dat", idWindTurbine);
    int nn2 = sprintf(buffer2, "FASTurbW_%d_towrInerMat.dat", idWindTurbine);
    int nn3 = sprintf(buffer3, "FASTurbW_%d_turbInerMat.dat", idWindTurbine);
    this->bodyInerMat.save(JoinPath(OutputPath, buffer1), arma::raw_ascii);
    this->towrInerMat.save(JoinPath(OutputPath, buffer2), arma::raw_ascii);
    this->turbInerMat.save(JoinPath(OutputPath, buffer3), arma::raw_ascii);
    std::cout << "        --> ... done!" << std::endl;

    std::cout << "        --> Get initial conditions" << std::endl;
    rotSpeed = FSTW_InitInput.turbIniRotSpeed;
    yaw = FSTW_InitInput.turbIniYaw;
    yaw_ini = yaw;
    YCMode = FSTW_InitInput.YCMode;
    if (FSTW_InitInput.isFixed_GenDOF > 0)
    {
        isRotorBlocked = true;
    }
    std::cout << "        --> ... done!" << std::endl;

    if (YCMode > 0)
    {
        std::stringstream ss;
        ss << "Option with YAW control not implemented yet, fix ServoDyn Input for wind turbine " << idWindTurbine << "\n";
        throw ValueError(ss.str());
    }

    // Calculate turbine at time 0.0...
    std::cout << "        --> Calculate turbine at time 0.0..." << std::endl;
    SetInputsFAST();
    std::cout << "            --> SetInputsFAST... done!" << std::endl;
    ComputeForces(0.0);
    std::cout << "            --> ComputeForces... done!" << std::endl;
    ComputeControler(0.0);
    std::cout << "            --> ComputeControler... done!" << std::endl;
    WriteOut(0.0);
    std::cout << "            --> WriteOut... done!" << std::endl;
    std::cout << "        --> ... done!" << std::endl;
}

void WindTurbine::Finalize(void)
{
    std::cout << "    WindTurbine::Finalize" << std::endl;
    FSTW_End(&ErrStat, ErrMsg);
    CheckError();
}

void WindTurbine::ComputeForces(double time)
{
    FSTW_CalcWind(&time, &ErrStat, ErrMsg);
    CheckError();
    FSTW_CalcForces(&time, &ErrStat, ErrMsg);
    CheckError();

    forceBodyCOG(0, 0) = FSTW_Output.plat_forces[0];
    forceBodyCOG(1, 0) = FSTW_Output.plat_forces[1];
    forceBodyCOG(2, 0) = FSTW_Output.plat_forces[2];
    forceBodyCOG(3, 0) = FSTW_Output.plat_forces[3];
    forceBodyCOG(4, 0) = FSTW_Output.plat_forces[4];
    forceBodyCOG(5, 0) = FSTW_Output.plat_forces[5];

    airTrq = FSTW_Output.aero_torque[0];

    // this->pBody->windTurbForces += forceBodyCOG;
}

void WindTurbine::SetInputsFAST(void)
{
    // std::cout << "    WindTurbine::SetBaseMovements" << std::endl;

    FSTW_Input.plat_pos[0] = pBody->pos(0, 0);
    FSTW_Input.plat_pos[1] = pBody->pos(1, 0);
    FSTW_Input.plat_pos[2] = pBody->pos(2, 0);
    FSTW_Input.plat_pos[3] = pBody->pos(3, 0);
    FSTW_Input.plat_pos[4] = pBody->pos(4, 0);
    FSTW_Input.plat_pos[5] = pBody->pos(5, 0);

    FSTW_Input.plat_vel[0] = pBody->vel(0, 0);
    FSTW_Input.plat_vel[1] = pBody->vel(1, 0);
    FSTW_Input.plat_vel[2] = pBody->vel(2, 0);
    FSTW_Input.plat_vel[3] = pBody->vel(3, 0);
    FSTW_Input.plat_vel[4] = pBody->vel(4, 0);
    FSTW_Input.plat_vel[5] = pBody->vel(5, 0);

    FSTW_Input.plat_acc[0] = pBody->acc(0, 0);
    FSTW_Input.plat_acc[1] = pBody->acc(1, 0);
    FSTW_Input.plat_acc[2] = pBody->acc(2, 0);
    FSTW_Input.plat_acc[3] = pBody->acc(3, 0);
    FSTW_Input.plat_acc[4] = pBody->acc(4, 0);
    FSTW_Input.plat_acc[5] = pBody->acc(5, 0);

    FSTW_Input.RotPos[0] = rotPos;
    FSTW_Input.RotSpeed[0] = rotSpeed;

    if (YCMode > 0)
    {
        FSTW_Input.Yaw[0] = yaw;
        FSTW_Input.YawSpeed[0] = yawSpeed;
    }
    else
    {
        FSTW_Input.Yaw[0] = yaw_ini;
        FSTW_Input.YawSpeed[0] = 0.0;
    }
}

void WindTurbine::ComputeControler(double time)
{
    // std::cout << "    WindTurbine::ComputeControler" << std::endl;
    FSTW_CalcController(&time, &ErrStat, ErrMsg);
    CheckError();

    genTrq = FSTW_Output.gen_torque[0];
    yawTrq = FSTW_Output.YawMom[0];
}

void WindTurbine::WriteOut(double time)
{
    // std::cout << "    WindTurbine::ComputeControler" << std::endl;
    FSTW_WriteOutput(&time, &ErrStat, ErrMsg);
    CheckError();
}

void WindTurbine::ComputeRotorAcc(void)
{
    // std::cout << "    WindTurbine::ComputeRotorAcc" << std::endl;

    if (isRotorBlocked)
    {
        rotAcc = 0.0;
    }
    else
    {
        rotAcc = (airTrq - genTrq) / rotIner;
    }
    // std::cout << "WindTurbine::ComputeRotorAcc - airTrq = " << airTrq << std::endl;
    // std::cout << "WindTurbine::ComputeRotorAcc - genTrq = " << genTrq << std::endl;
}

void WindTurbine::CheckError(void)
{
    // std::cout << "    WindTurbine::CheckError" << std::endl;
    int ErrStat2;
    char ErrMsg2[INTERFACE_STRING_LENGTH];
    if (ErrStat != ErrID_None)
    {
        if (ErrStat >= AbortErrLev)
        {
            FSTW_End(&ErrStat2, ErrMsg2);
            if (ErrStat != ErrStat2)
            {
                std::cout << "ERROR: " << ErrMsg2 << std::endl;
            }
            throw std::runtime_error(ErrMsg);
        }
    }
}

#endif // OASIS_USE_OPENFAST