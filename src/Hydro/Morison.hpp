// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef morisondef_hpp__
#define morisondef_hpp__
#include <armadillo>
#include <string>
#include <cstdio>

class Simulation;

class Morison
{
public:
    // Pointer to simulation — provides access to inputFolderPath, dataFormat, yamlRoot
    Simulation* pSim;
    int numBodies;
    double pi;

    // Per-body flags: flag_wind[i] / flag_curr[i] are true only if body i has active flow
    bool* flag_wind;
    bool* flag_curr;

    // Per-body constant flow scalars (one value per body)
    arma::vec wind_spd_body;
    arma::vec wind_dir_body;
    arma::vec curr_spd_body;
    arma::vec curr_dir_body;

    // Per-body heading arrays and drag/FK coefficient cubes
    // pHeadings[i]: (5,) linspace(0,360,5) when SymOrder=2
    arma::mat* pHeadings;
    arma::cube** pWindFKCoeff;   // [numBodies][5 headings x 6 dofs x 2 vel_comps]
    arma::cube** pWindDragCoeff;
    arma::cube** pCurrFKCoeff;
    arma::cube** pCurrDragCoeff;

    // Constructor
    Morison(int numBodies_inp, Simulation* pSim_inp);

    // Top-level read dispatcher (routes to ASCII or YAML based on pSim->GetDataFormat())
    void ReadMorisonData(void);
    void ReadMorisonDataASCII(void);
    void ReadMorisonDataYAML(void);

    // Shared helper: store parsed data for one body index (0-based)
    void StoreBodyMorisonData(int bodyIdx, double wspd, double wdir, double cspd, double cdir,
                              const arma::mat& windFK_X, const arma::mat& windDrag_X,
                              const arma::mat& windFK_Y, const arma::mat& windDrag_Y,
                              const arma::mat& currFK_X, const arma::mat& currDrag_X,
                              const arma::mat& currFK_Y, const arma::mat& currDrag_Y);

    void ReadFlowData_HDF5(void);
    void ReadMorCoeffData_HDF5(void);
    arma::mat ComputeWindForce(int idBody, double yaw, double t);
    arma::mat ComputeCurrForce(int idBody, double yaw, double t);
};

#endif // morisondef_hpp__