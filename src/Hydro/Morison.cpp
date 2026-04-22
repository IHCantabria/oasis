// SPDX-License-Identifier: GPL-3.0-or-later
#include <armadillo>
#include <string>
#include <fstream>
#include <limits>
#include <sstream>
#include "Morison.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"
#include "../Logger.hpp"

Morison::Morison(int numBodies_inp, Simulation* pSim_inp)
{
    pi = arma::datum::pi;
    numBodies = numBodies_inp;
    pSim = pSim_inp;

    flag_wind = new bool[numBodies]();  // value-initialised to false
    flag_curr = new bool[numBodies]();

    wind_spd_body = arma::zeros<arma::vec>(numBodies);
    wind_dir_body = arma::zeros<arma::vec>(numBodies);
    curr_spd_body = arma::zeros<arma::vec>(numBodies);
    curr_dir_body = arma::zeros<arma::vec>(numBodies);

    pHeadings = new arma::mat[numBodies];
    pWindFKCoeff   = new arma::cube*[numBodies];
    pWindDragCoeff = new arma::cube*[numBodies];
    pCurrFKCoeff   = new arma::cube*[numBodies];
    pCurrDragCoeff = new arma::cube*[numBodies];
    for (int ii = 0; ii < numBodies; ii++)
    {
        pWindFKCoeff[ii]   = nullptr;
        pWindDragCoeff[ii] = nullptr;
        pCurrFKCoeff[ii]   = nullptr;
        pCurrDragCoeff[ii] = nullptr;
    }
}

void Morison::ReadMorisonData(void)
{
    if (pSim->GetDataFormat() == 0)
        ReadMorisonDataASCII();
    else
        ReadMorisonDataYAML();
}

void Morison::ReadMorisonDataASCII(void)
{
    std::string file_path = JoinPath(pSim->inputFolderPath, "dataMorison.dat");
    std::ifstream f(file_path);
    if (!f.is_open())
    {
        Logger::warning("dataMorison.dat was not found! Setting wind and currents to zero for all bodies.");
        return;
    }

    // Helpers that read one token + rest-of-line, or skip a full line
    auto skipLine = [&]() { f.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); };
    auto readInt  = [&](int& v)    { f >> v;  f.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); };
    auto readDbl  = [&](double& v) { f >> v;  f.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); };
    auto readStr  = [&](std::string& v) { f >> v; f.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); };
    auto readMat6x2 = [&](arma::mat& M) {
        skipLine(); // header comment line
        for (int r = 0; r < 6; r++)
        {
            f >> M(r, 0) >> M(r, 1);
            f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    };

    // Read number of body sections
    int numSections;
    readInt(numSections);

    for (int s = 0; s < numSections; s++)
    {
        // Skip 3 header lines  (//, // Body [N], //)
        skipLine(); skipLine(); skipLine();

        // Body index (1-based → 0-based)
        int bodyIndex1;
        readInt(bodyIndex1);
        int idx = bodyIndex1 - 1;
        if (idx < 0 || idx >= numBodies)
        {
            std::stringstream ss;
            ss << "dataMorison.dat: body index " << bodyIndex1
               << " is out of range (numBodies=" << numBodies << ").\n";
            throw ValueError(ss.str());
        }

        // --- WIND AND CURRENT section ---
        // 3 header lines (// ----, // WIND AND CURRENT DEFINITION, // ----)
        skipLine(); skipLine(); skipLine();

        int flowType;
        readInt(flowType);
        if (flowType != 1)
        {
            std::stringstream ss;
            ss << "dataMorison.dat: flow type " << flowType
               << " is not implemented. Only FlowType=1 (constant) is supported.\n";
            throw NotImplementedError(ss.str());
        }

        // // Constant flow data case
        skipLine();
        double wspd, wdir, cspd, cdir;
        readDbl(wspd);
        readDbl(wdir);
        readDbl(cspd);
        readDbl(cdir);

        // // Variable flow data case  +  placeholder filename (unused)
        skipLine();
        std::string dummy;
        readStr(dummy);

        // --- MORISON COEFFICIENTS section ---
        // 4 header lines (// ----, // MORISON COEFFICIENTS DEFINITION, // Must include..., // ----)
        skipLine(); skipLine(); skipLine(); skipLine();

        int symFlag;
        readInt(symFlag);
        if (symFlag != 2)
        {
            std::stringstream ss;
            ss << "dataMorison.dat: Sym_flag=" << symFlag
               << " is not implemented. Only Sym_flag=2 (symmetric) is supported.\n";
            throw NotImplementedError(ss.str());
        }

        // // Non-symmetric case comment  +  placeholder filename
        skipLine();
        readStr(dummy);

        // // Symmetric case comment
        skipLine();

        int symOrder;
        readInt(symOrder);
        if (symOrder != 2)
        {
            std::stringstream ss;
            ss << "dataMorison.dat: SymOrder=" << symOrder
               << " is not implemented. Only SymOrder=2 is supported.\n";
            throw NotImplementedError(ss.str());
        }

        arma::mat windFK_X  = arma::zeros(6, 2);
        arma::mat windDrag_X = arma::zeros(6, 2);
        arma::mat windFK_Y  = arma::zeros(6, 2);
        arma::mat windDrag_Y = arma::zeros(6, 2);
        arma::mat currFK_X  = arma::zeros(6, 2);
        arma::mat currDrag_X = arma::zeros(6, 2);
        arma::mat currFK_Y  = arma::zeros(6, 2);
        arma::mat currDrag_Y = arma::zeros(6, 2);

        readMat6x2(windFK_X);
        readMat6x2(windDrag_X);
        readMat6x2(windFK_Y);
        readMat6x2(windDrag_Y);
        readMat6x2(currFK_X);
        readMat6x2(currDrag_X);
        readMat6x2(currFK_Y);
        readMat6x2(currDrag_Y);

        StoreBodyMorisonData(idx, wspd, wdir, cspd, cdir,
                             windFK_X, windDrag_X, windFK_Y, windDrag_Y,
                             currFK_X, currDrag_X, currFK_Y, currDrag_Y);
    }

    f.close();
}

void Morison::ReadMorisonDataYAML(void)
{
    YAML::Node yamlRoot = pSim->GetYamlRoot();
    if (!yamlRoot["morison"])
    {
        Logger::warning("'morison' section not found in YAML. Setting wind and currents to zero for all bodies.");
        return;
    }

    YAML::Node morNode = yamlRoot["morison"];

    for (std::size_t s = 0; s < morNode.size(); s++)
    {
        YAML::Node entry = morNode[s];

        int bodyIndex1 = entry["body_index"].as<int>();
        int idx = bodyIndex1 - 1;
        if (idx < 0 || idx >= numBodies)
        {
            std::stringstream ss;
            ss << "dataProblem.yaml morison: body_index " << bodyIndex1
               << " is out of range (numBodies=" << numBodies << ").\n";
            throw ValueError(ss.str());
        }

        int flowType = entry["flow_type"].as<int>();
        if (flowType != 1)
        {
            std::stringstream ss;
            ss << "dataProblem.yaml morison: flow_type=" << flowType
               << " is not implemented. Only flow_type=1 (constant) is supported.\n";
            throw NotImplementedError(ss.str());
        }

        int symOrder = entry["symmetry_order"].as<int>();
        if (symOrder != 2)
        {
            std::stringstream ss;
            ss << "dataProblem.yaml morison: symmetry_order=" << symOrder
               << " is not implemented. Only symmetry_order=2 is supported.\n";
            throw NotImplementedError(ss.str());
        }

        double wspd = entry["wind_speed"].as<double>();
        double wdir = entry["wind_direction"].as<double>();
        double cspd = entry["current_speed"].as<double>();
        double cdir = entry["current_direction"].as<double>();

        // Read each 6x2 matrix from a nested YAML sequence [[a,b],[c,d],...]
        auto readYAMLMat6x2 = [](YAML::Node node) -> arma::mat {
            arma::mat M = arma::zeros(6, 2);
            for (int r = 0; r < 6; r++)
            {
                M(r, 0) = node[r][0].as<double>();
                M(r, 1) = node[r][1].as<double>();
            }
            return M;
        };

        arma::mat windFK_X   = readYAMLMat6x2(entry["wind_fk_x"]);
        arma::mat windDrag_X = readYAMLMat6x2(entry["wind_drag_x"]);
        arma::mat windFK_Y   = readYAMLMat6x2(entry["wind_fk_y"]);
        arma::mat windDrag_Y = readYAMLMat6x2(entry["wind_drag_y"]);
        arma::mat currFK_X   = readYAMLMat6x2(entry["current_fk_x"]);
        arma::mat currDrag_X = readYAMLMat6x2(entry["current_drag_x"]);
        arma::mat currFK_Y   = readYAMLMat6x2(entry["current_fk_y"]);
        arma::mat currDrag_Y = readYAMLMat6x2(entry["current_drag_y"]);

        StoreBodyMorisonData(idx, wspd, wdir, cspd, cdir,
                             windFK_X, windDrag_X, windFK_Y, windDrag_Y,
                             currFK_X, currDrag_X, currFK_Y, currDrag_Y);
    }
}

void Morison::StoreBodyMorisonData(int bodyIdx, double wspd, double wdir, double cspd, double cdir,
                                   const arma::mat& windFK_X,  const arma::mat& windDrag_X,
                                   const arma::mat& windFK_Y,  const arma::mat& windDrag_Y,
                                   const arma::mat& currFK_X,  const arma::mat& currDrag_X,
                                   const arma::mat& currFK_Y,  const arma::mat& currDrag_Y)
{
    wind_spd_body(bodyIdx) = wspd;
    wind_dir_body(bodyIdx) = wdir;
    curr_spd_body(bodyIdx) = cspd;
    curr_dir_body(bodyIdx) = cdir;

    // SymOrder=2: headings at 0, 90, 180, 270, 360 (5 points)
    pHeadings[bodyIdx] = arma::linspace(0.0, 360.0, 5);

    // Build 5×6×2 cubes by assigning X-heading data to even indices, Y-heading to odd indices
    arma::cube WindFKCoeff   = arma::zeros(5, 6, 2);
    arma::cube WindDragCoeff = arma::zeros(5, 6, 2);
    arma::cube CurrFKCoeff   = arma::zeros(5, 6, 2);
    arma::cube CurrDragCoeff = arma::zeros(5, 6, 2);

    for (int ii = 0; ii < 5; ii += 2)
    {
        WindFKCoeff  (arma::span(ii), arma::span::all, arma::span::all) = windFK_X;
        WindDragCoeff(arma::span(ii), arma::span::all, arma::span::all) = windDrag_X;
        CurrFKCoeff  (arma::span(ii), arma::span::all, arma::span::all) = currFK_X;
        CurrDragCoeff(arma::span(ii), arma::span::all, arma::span::all) = currDrag_X;
    }
    for (int ii = 1; ii < 5; ii += 2)
    {
        WindFKCoeff  (arma::span(ii), arma::span::all, arma::span::all) = windFK_Y;
        WindDragCoeff(arma::span(ii), arma::span::all, arma::span::all) = windDrag_Y;
        CurrFKCoeff  (arma::span(ii), arma::span::all, arma::span::all) = currFK_Y;
        CurrDragCoeff(arma::span(ii), arma::span::all, arma::span::all) = currDrag_Y;
    }

    pWindFKCoeff[bodyIdx]   = new arma::cube(WindFKCoeff);
    pWindDragCoeff[bodyIdx] = new arma::cube(WindDragCoeff);
    pCurrFKCoeff[bodyIdx]   = new arma::cube(CurrFKCoeff);
    pCurrDragCoeff[bodyIdx] = new arma::cube(CurrDragCoeff);

    if (wspd > 0.0)
    {
        Logger::info("    --> Morison wind forces activated for body " + std::to_string(bodyIdx + 1) + ".");
        flag_wind[bodyIdx] = true;
    }
    if (cspd > 0.0)
    {
        Logger::info("    --> Morison current forces activated for body " + std::to_string(bodyIdx + 1) + ".");
        flag_curr[bodyIdx] = true;
    }
}

void Morison::ReadFlowData_HDF5(void)
{
    Logger::info("--> Reading Flow (HDF5 format)");
    std::stringstream ss;
    ss << "Method ReadFlowData_HDF5 in class Morison not implemented yet. \n";
    throw NotImplementedError(ss.str());
}

void Morison::ReadMorCoeffData_HDF5(void)
{
    Logger::info("--> Reading Morison Coefficients (HDF5 format)");
    std::stringstream ss;
    ss << "Method ReadMorCoeffData_HDF5 in class Morison not implemented yet. \n";
    throw NotImplementedError(ss.str());
}

arma::mat Morison::ComputeWindForce(int idBody, double yaw, double t)
{
    arma::mat F = arma::zeros(6, 1);

    double wspd = wind_spd_body(idBody);
    double wdir = wind_dir_body(idBody);

    arma::mat vel = arma::zeros(2, 1);
    vel(0, 0) = wspd * cos(pi * wdir / 180.0);
    vel(1, 0) = wspd * sin(pi * wdir / 180.0);

    arma::mat h = arma::zeros(1, 1) + yaw + wdir;

    arma::cube temp_B = interp1(pHeadings[idBody], *(pWindDragCoeff[idBody]), h);
    arma::mat B = temp_B(arma::span(0), arma::span::all, arma::span::all);

    F = F + B * (vel % arma::abs(vel));
    return F;
}

arma::mat Morison::ComputeCurrForce(int idBody, double yaw, double t)
{
    arma::mat F = arma::zeros(6, 1);

    double cspd = curr_spd_body(idBody);
    double cdir = curr_dir_body(idBody);

    arma::mat vel = arma::zeros(2, 1);
    vel(0, 0) = cspd * cos(pi * cdir / 180.0);
    vel(1, 0) = cspd * sin(pi * cdir / 180.0);

    arma::mat h = arma::zeros(1, 1) + cdir;

    arma::cube temp_B = interp1(pHeadings[idBody], *(pCurrDragCoeff[idBody]), h);
    arma::mat B = temp_B(arma::span(0), arma::span::all, arma::span::all);

    F = F + B * (vel % arma::abs(vel));
    return F;
}