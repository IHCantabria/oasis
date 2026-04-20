#ifndef LINETYPE_HPP__
#define LINETYPE_HPP__

#include <armadillo>
#include <stdio.h>
#include <yaml-cpp/yaml.h>

struct LineType
{
    int flag_tension;
    double rho0, d;
    int flag_stiffness;

    // flag_stiffness == 0: constant axial stiffness
    double EA, beta;

    // flag_stiffness == 1: viscoelastic
    int num_kernel_coef, num_elastic_coef;
    arma::mat kernel_lin_coef;
    arma::mat kernel_exp_coef;
    arma::mat elastic_coef;

    // flag_stiffness > 1: tabulated (value = number of data points)
    arma::mat strain_data, stress_data;

    // Hydrodynamic / seabed coefficients
    double CB, Cmn, Cdn, Cdt, GK, GC;

    // Contact / friction model
    int smoothstep;
    int frictionModel;
    double vth, ust, usn, ud, deltamax;

    // Friction polynomial coefficients (computed when frictionModel == 1)
    arma::mat a_1 = arma::zeros(4, 1);

    void ReadPropertiesASCII(FILE* fp);
    void ReadPropertiesYAML(YAML::Node node);
};

#endif // LINETYPE_HPP__
