// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef SPRINGTYPE_HPP__
#define SPRINGTYPE_HPP__

#include <armadillo>
#include <stdio.h>
#include <yaml-cpp/yaml.h>

struct SpringType
{
    int stressModelFlag;
    int dampingFlag;
    int frictionFlag;
    int frameFlag;

    arma::mat SpringMatrix_K = arma::zeros(6, 6);
    arma::mat SpringMatrix_D = arma::zeros(6, 1);
    arma::mat SpringMatrix_M = arma::zeros(6, 6);

    double mu_d, mu_s, vt, Dt;

    // Friction polynomial coefficients (computed from mu_s and vt)
    arma::mat a_1 = arma::zeros(4, 1);
    arma::mat a_2 = arma::zeros(4, 1);

    int n_StressStrain[6];
    arma::field<arma::mat> data_StressStrain;

    void ReadPropertiesASCII(FILE* fp);
    void ReadPropertiesYAML(YAML::Node node);
};

#endif // SPRINGTYPE_HPP__
