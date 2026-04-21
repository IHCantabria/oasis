// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SPRING_FLAG
#define SPRING_FLAG
#include <armadillo>
#include <string>
#include <cstdio>
#include <yaml-cpp/yaml.h>
#include "../BCPs/BCPs.hpp"
#include "SpringType.hpp"

class Simulation;

class Spring
{
public:
    int nSpring;   // Spring number
    int typeIndex; // 0-based index into pSpringTypes array

    int stressModelFlag; // Flag that determines the stress model: [1: linear matrix, 2: stress-strain curve]
    int dampingFlag;     // Flag that determines whether damping is considered: [0: no, 1: yes]
    int frictionFlag;    // Flag that determines whether friction is considered: [0: no, 1: yes]
    int frameFlag; // Flag that indicates which frame is used for strain/force calculation: [0: mean, 1: BCP1, 2: BCP2]

    int BCP_1; // Indices of the BCPs connected to this spring
    int BCP_2;
    int BCP_1_type; // Flags indicating the type of each BCP
    int BCP_2_type;
    BCP* SpringBCP[2]; // Boundary condition points attached to this spring

    arma::field<arma::mat> SpringVectors; // Unit vectors defining the spring orientation in local frame for each BCP
                                          // (order: normal (x), tangent 1 (y), tangent 2 (z))
    arma::mat SpringMatrix_K = arma::zeros(6, 6); // Spring stiffness matrix
    double mu_d, mu_s, vt,
        Dt;                // Dynamic friction coefficient, static friction coefficient, and maximum friction velocity
    int flagStickSlip = 0; // [0: Slip, 1: Stick]
    arma::mat SpringStrains_Stick;
    arma::mat a_1 = arma::zeros(4, 1),
              a_2 = arma::zeros(4,
                                1);               // Polynomial coefficients for the friction coefficient curve
    arma::mat SpringMatrix_M = arma::zeros(6, 6); // Spring friction matrix
    arma::mat SpringMatrix_D = arma::zeros(6, 1); // Spring damping coefficients
    // For 6 degrees of freedom: number of data points, n, in the stress-strain curves
    int n_StressStrain[6];
    // For 6 degrees of freedom: nx1 matrix of strain states and nx6 matrix of corresponding stresses
    arma::field<arma::mat> data_StressStrain;

    Spring(int n)
    {
        nSpring = n;
    } // Initialises a spring object with the given index
    void ReadPropertiesASCII(FILE* pFilePointer, SpringType** pTypes,
                             int numTypes); // Read spring inputs from ASCII file
    void ReadPropertiesYAML(YAML::Node node, SpringType** pTypes, int numTypes);
    void computeSpringForces(void); // Computes spring forces at each BCP and stores them in the BCP objects

    // Output
    Simulation* pSim = nullptr;
    FILE* pfile_spring_csv = nullptr;
    arma::mat lastForceG = arma::zeros(6, 1); // Last computed global force for output
    void OpenOutputFiles(std::string path);
    void CloseOutputFiles(void);
    void WriteOut(double t);
};

#endif