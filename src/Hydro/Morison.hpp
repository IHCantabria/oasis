#ifndef morisondef_hpp__
#define morisondef_hpp__
#include <armadillo>
#include <string>
#include <cstdio>

class Simulation;

class Morison
{
public:
    // Declare class variables
    Simulation *pSim; // Pointer to the simulation instance. It gives fast access to the necessary simulation variables
    int numBodies;
    double pi;

    std::string FlowDataFile;
    std::string MorCoeffDataFile;

    bool flag_wind = false;
    bool flag_curr = false;

    int FlowType_flag;
    arma::mat time, wind_spd, wind_dir, curr_spd, curr_dir, wind_acc, curr_acc;
    arma::mat headings;
    arma::cube **pWindFKCoeff;   // Matrix components: [body, headings, dofs, vel_comps];
    arma::cube **pWindDragCoeff; // Matrix components: [body, headings, dofs, vel_comps];
    arma::cube **pCurrFKCoeff;   // Matrix components: [body, headings, dofs, vel_comps];
    arma::cube **pCurrDragCoeff; // Matrix components: [body, headings, dofs, vel_comps];

    // Declare class constructors
    Morison(int numBodies_inp, Simulation *pSim_inp);

    // Methods
    void ReadMorisonData(void);
    void ReadFlowData_HDF5(void);
    void ReadMorCoeffData_HDF5(void);
    arma::mat ComputeCurrForce(int idBody, double yaw, double t);
    arma::mat ComputeWindForce(int idBody, double yaw, double t);
};

#endif // morisondef_hpp__