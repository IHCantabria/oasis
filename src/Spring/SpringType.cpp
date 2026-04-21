// SPDX-License-Identifier: GPL-3.0-or-later
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <armadillo>
#include "SpringType.hpp"
#include "../Exceptions/Exception.hpp"

void SpringType::ReadPropertiesASCII(FILE* fp)
{
    char buffer_line[1000];
    int jj, kk, ll, temp_N;

    data_StressStrain.set_size(6, 2);

    // Skip 3 header lines ("New spring type [N]")
    for (int ii = 0; ii < 3; ii++)
        fgets(buffer_line, sizeof(buffer_line), fp);

    fscanf(fp, "%d %[^\n]\n", &stressModelFlag, buffer_line);
    fscanf(fp, "%d %[^\n]\n", &dampingFlag, buffer_line);
    fscanf(fp, "%d %[^\n]\n", &frictionFlag, buffer_line);
    fscanf(fp, "%d %[^\n]\n", &frameFlag, buffer_line);

    // Skip stiffness matrix comment line
    fgets(buffer_line, sizeof(buffer_line), fp);
    for (jj = 0; jj < 6; jj++)
    {
        for (kk = 0; kk < 6; kk++)
            fscanf(fp, "%lf", &SpringMatrix_K(jj, kk));
        fscanf(fp, "%[^\n]\n", buffer_line);
    }

    // Skip friction coefficients comment line
    fgets(buffer_line, sizeof(buffer_line), fp);
    fscanf(fp, "%lf %[^\n]\n", &mu_d, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &mu_s, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &vt, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &Dt, buffer_line);

    // Compute friction polynomial coefficients
    arma::mat tmpA, tmpB = arma::zeros(4, 1);
    tmpA = {{pow(vt / 2, 3), pow(vt / 2, 2), vt / 2, 1},
            {3 * pow(vt / 2, 2), 2 * (vt / 2), 1, 0},
            {pow(vt, 3), pow(vt, 2), vt, 1},
            {3 * pow(vt, 2), 2 * vt, 1, 0}};
    tmpB(0, 0) = mu_s / 2;
    tmpB(1, 0) = mu_s / vt;
    tmpB(2, 0) = mu_s;
    tmpB(3, 0) = 0;
    a_1 = arma::solve(tmpA, tmpB);

    tmpA = {{pow(vt, 3), pow(vt, 2), vt, 1},
            {3 * pow(vt, 2), 2 * vt, 1, 0},
            {pow(vt * 2, 3), pow(vt * 2, 2), vt * 2, 1},
            {3 * pow(vt * 2, 2), 2 * (vt * 2), 1, 0}};
    tmpB(0, 0) = mu_s;
    tmpB(1, 0) = 0;
    tmpB(2, 0) = mu_d;
    tmpB(3, 0) = 0;
    a_2 = arma::solve(tmpA, tmpB);

    // Skip friction matrix comment line
    fgets(buffer_line, sizeof(buffer_line), fp);
    for (jj = 0; jj < 6; jj++)
    {
        for (kk = 0; kk < 6; kk++)
            fscanf(fp, "%lf", &SpringMatrix_M(jj, kk));
        fscanf(fp, "%[^\n]\n", buffer_line);
    }

    // Skip damping coefficients comment line
    fgets(buffer_line, sizeof(buffer_line), fp);
    for (jj = 0; jj < 6; jj++)
        fscanf(fp, "%lf %[^\n]\n", &SpringMatrix_D(jj, 0), buffer_line);

    // Read stress-strain curves for 6 DOFs
    for (jj = 0; jj < 6; jj++)
    {
        // Skip DOF comment line
        fgets(buffer_line, sizeof(buffer_line), fp);
        fscanf(fp, "%d %[^\n]\n", &temp_N, buffer_line);
        n_StressStrain[jj] = temp_N;

        arma::mat temp_vec2 = arma::zeros(temp_N, 1);
        for (ll = 0; ll < temp_N; ll++)
            fscanf(fp, "%lf", &temp_vec2(ll, 0));
        fscanf(fp, "%[^\n]\n", buffer_line);
        data_StressStrain(jj, 0) = temp_vec2;

        arma::mat temp_mat = arma::zeros(temp_N, 6);
        for (kk = 0; kk < 6; kk++)
        {
            for (ll = 0; ll < temp_N; ll++)
                fscanf(fp, "%lf", &temp_mat(ll, kk));
            fscanf(fp, "%[^\n]\n", buffer_line);
        }
        data_StressStrain(jj, 1) = temp_mat;
    }
}

void SpringType::ReadPropertiesYAML(YAML::Node node)
{
    data_StressStrain.set_size(6, 2);

    stressModelFlag = node["stress_model_flag"].as<int>();
    dampingFlag = node["damping_flag"].as<int>();
    frictionFlag = node["friction_flag"].as<int>();
    frameFlag = node["frame_flag"].as<int>();

    YAML::Node kNode = node["stiffness_matrix"];
    for (int jj = 0; jj < 6; jj++)
        for (int kk = 0; kk < 6; kk++)
            SpringMatrix_K(jj, kk) = kNode[jj][kk].as<double>();

    mu_d = node["mu_d"].as<double>();
    mu_s = node["mu_s"].as<double>();
    vt = node["vt"].as<double>();
    Dt = node["Dt"].as<double>();

    arma::mat tmpA, tmpB = arma::zeros(4, 1);
    tmpA = {{pow(vt / 2, 3), pow(vt / 2, 2), vt / 2, 1},
            {3 * pow(vt / 2, 2), 2 * (vt / 2), 1, 0},
            {pow(vt, 3), pow(vt, 2), vt, 1},
            {3 * pow(vt, 2), 2 * vt, 1, 0}};
    tmpB(0, 0) = mu_s / 2;
    tmpB(1, 0) = mu_s / vt;
    tmpB(2, 0) = mu_s;
    tmpB(3, 0) = 0;
    a_1 = arma::solve(tmpA, tmpB);

    tmpA = {{pow(vt, 3), pow(vt, 2), vt, 1},
            {3 * pow(vt, 2), 2 * vt, 1, 0},
            {pow(vt * 2, 3), pow(vt * 2, 2), vt * 2, 1},
            {3 * pow(vt * 2, 2), 2 * (vt * 2), 1, 0}};
    tmpB(0, 0) = mu_s;
    tmpB(1, 0) = 0;
    tmpB(2, 0) = mu_d;
    tmpB(3, 0) = 0;
    a_2 = arma::solve(tmpA, tmpB);

    YAML::Node mNode = node["mass_matrix"];
    for (int jj = 0; jj < 6; jj++)
        for (int kk = 0; kk < 6; kk++)
            SpringMatrix_M(jj, kk) = mNode[jj][kk].as<double>();

    YAML::Node dNode = node["damping_vector"];
    for (int jj = 0; jj < 6; jj++)
        SpringMatrix_D(jj, 0) = dNode[jj].as<double>();

    YAML::Node ssNode = node["stress_strain"];
    for (int jj = 0; jj < 6; jj++)
    {
        YAML::Node dofNode = ssNode[jj];
        int temp_N = (int)dofNode["displacements"].size();
        n_StressStrain[jj] = temp_N;

        arma::mat temp_vec2 = arma::zeros(temp_N, 1);
        for (int ll = 0; ll < temp_N; ll++)
            temp_vec2(ll, 0) = dofNode["displacements"][ll].as<double>();
        data_StressStrain(jj, 0) = temp_vec2;

        arma::mat temp_mat = arma::zeros(temp_N, 6);
        YAML::Node forcesNode = dofNode["forces"];
        for (int kk = 0; kk < 6; kk++)
            for (int ll = 0; ll < temp_N; ll++)
                temp_mat(ll, kk) = forcesNode[kk][ll].as<double>();
        data_StressStrain(jj, 1) = temp_mat;
    }
}
