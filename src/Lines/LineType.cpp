#include <cmath>
#include <sstream>
#include <armadillo>
#include "LineType.hpp"
#include "../Exceptions/Exception.hpp"

void LineType::ReadPropertiesASCII(FILE* fp)
{
    char buffer_line[1000];
    double dtemp;
    fpos_t carriage_init;

    // Skip 3 header lines ("New line type [N]")
    for (int ii = 0; ii < 3; ii++)
        fgets(buffer_line, sizeof(buffer_line), fp);

    fscanf(fp, "%d %[^\n]\n", &flag_tension, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &rho0, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &d, buffer_line);
    fscanf(fp, "%d %[^\n]\n", &flag_stiffness, buffer_line);

    if (flag_stiffness == 0)
    {
        fscanf(fp, "%lf %[^\n]\n", &EA, buffer_line);
        fscanf(fp, "%lf %[^\n]\n", &beta, buffer_line);
    }
    else if (flag_stiffness == 1)
    {
        num_kernel_coef = 0;
        fgetpos(fp, &carriage_init);
        while (fscanf(fp, "%lf", &dtemp) == 1)
            num_kernel_coef++;
        fsetpos(fp, &carriage_init);
        kernel_lin_coef = arma::zeros(num_kernel_coef, 1);
        for (int ii = 0; ii < num_kernel_coef; ii++)
        {
            if (fscanf(fp, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "Error reading kernel linear coefficients for line type.\n";
                throw ValueError(ss.str());
            }
            kernel_lin_coef(ii, 0) = dtemp;
        }
        fscanf(fp, "%[^\n]\n", buffer_line);

        int num_kernel_coef_tmp = 0;
        fgetpos(fp, &carriage_init);
        while (fscanf(fp, "%lf", &dtemp) == 1)
            num_kernel_coef_tmp++;
        if (num_kernel_coef_tmp != num_kernel_coef)
        {
            std::stringstream ss;
            ss << "Number of linear and exponential kernel coefficients must be equal in line type.\n";
            throw ValueError(ss.str());
        }
        fsetpos(fp, &carriage_init);
        kernel_exp_coef = arma::zeros(num_kernel_coef, 1);
        for (int ii = 0; ii < num_kernel_coef; ii++)
        {
            if (fscanf(fp, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "Error reading kernel exponential coefficients for line type.\n";
                throw ValueError(ss.str());
            }
            kernel_exp_coef(ii, 0) = dtemp;
        }
        fscanf(fp, "%[^\n]\n", buffer_line);

        num_elastic_coef = 0;
        fgetpos(fp, &carriage_init);
        while (fscanf(fp, "%lf", &dtemp) == 1)
            num_elastic_coef++;
        fsetpos(fp, &carriage_init);
        elastic_coef = arma::zeros(num_elastic_coef, 1);
        for (int ii = 0; ii < num_elastic_coef; ii++)
        {
            if (fscanf(fp, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "Error reading elastic polynomial coefficients for line type.\n";
                throw ValueError(ss.str());
            }
            elastic_coef(ii, 0) = dtemp;
        }
        fscanf(fp, "%[^\n]\n", buffer_line);

        if (elastic_coef(0, 0) > 0.0)
            EA = elastic_coef(0, 0);
        else
        {
            double strain_EA = 0.1;
            EA = elastic_coef(0, 0) + 2 * elastic_coef(1, 0) * strain_EA +
                 3 * elastic_coef(2, 0) * strain_EA * strain_EA;
        }
    }
    else if (flag_stiffness > 1)
    {
        strain_data = arma::zeros(flag_stiffness, 1);
        stress_data = arma::zeros(flag_stiffness, 1);
        for (int ii = 0; ii < flag_stiffness; ii++)
        {
            if (fscanf(fp, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "Error reading strain data for line type.\n";
                throw ValueError(ss.str());
            }
            strain_data(ii, 0) = dtemp;
        }
        fscanf(fp, "%[^\n]\n", buffer_line);
        for (int ii = 0; ii < flag_stiffness; ii++)
        {
            if (fscanf(fp, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "Error reading stress data for line type.\n";
                throw ValueError(ss.str());
            }
            stress_data(ii, 0) = dtemp;
        }
        fscanf(fp, "%[^\n]\n", buffer_line);
        fscanf(fp, "%lf %[^\n]\n", &beta, buffer_line);

        EA = 0;
        int i_EA = 0;
        while (EA <= 0)
        {
            EA = (stress_data(i_EA + 1, 0) - stress_data(i_EA, 0)) / (strain_data(i_EA + 1, 0) - strain_data(i_EA, 0));
            i_EA++;
            if (i_EA > flag_stiffness - 1)
            {
                std::stringstream ss;
                ss << "Strain-Stress curve should have positive slope for line type.\n";
                throw ValueError(ss.str());
            }
        }
    }

    fscanf(fp, "%lf %[^\n]\n", &CB, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &Cmn, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &Cdn, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &Cdt, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &GK, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &GC, buffer_line);
    fscanf(fp, "%d %[^\n]\n", &smoothstep, buffer_line);
    fscanf(fp, "%d %[^\n]\n", &frictionModel, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &vth, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &ust, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &usn, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &ud, buffer_line);
    fscanf(fp, "%lf %[^\n]\n", &deltamax, buffer_line);

    if (frictionModel == 1)
    {
        arma::mat tmpA = arma::zeros(4, 4);
        arma::mat tmpB = arma::zeros(4, 1);
        double x1 = 0;
        double x2 = 1e-04;
        tmpA = {{pow(x1, 3), pow(x1, 2), 1 * x1, 1},
                {3 * pow(x1, 2), 2 * (x1), 1, 0},
                {pow(x2, 3), pow(x2, 2), 1 * x2, 1},
                {3 * pow(x2, 2), 2 * x2, 1, 0}};
        tmpB(0, 0) = 1e-07;
        tmpB(1, 0) = 0.0;
        tmpB(2, 0) = x2;
        tmpB(3, 0) = 1;
        a_1 = arma::solve(tmpA, tmpB);
    }
}

void LineType::ReadPropertiesYAML(YAML::Node node)
{
    flag_tension = node["flag_tension"].as<int>();
    rho0 = node["density"].as<double>();
    d = node["diameter"].as<double>();
    flag_stiffness = node["flag_stiffness"].as<int>();

    if (flag_stiffness == 0)
    {
        EA = node["EA"].as<double>();
        beta = node["beta"].as<double>();
    }
    else if (flag_stiffness == 1)
    {
        YAML::Node klNode = node["kernel_lin_coef"];
        num_kernel_coef = (int)klNode.size();
        kernel_lin_coef = arma::zeros(num_kernel_coef, 1);
        for (int ii = 0; ii < num_kernel_coef; ii++)
            kernel_lin_coef(ii, 0) = klNode[ii].as<double>();

        YAML::Node keNode = node["kernel_exp_coef"];
        if ((int)keNode.size() != num_kernel_coef)
        {
            std::stringstream ss;
            ss << "Number of linear and exponential kernel coefficients must be equal in line type.\n";
            throw ValueError(ss.str());
        }
        kernel_exp_coef = arma::zeros(num_kernel_coef, 1);
        for (int ii = 0; ii < num_kernel_coef; ii++)
            kernel_exp_coef(ii, 0) = keNode[ii].as<double>();

        YAML::Node ecNode = node["elastic_coef"];
        num_elastic_coef = (int)ecNode.size();
        elastic_coef = arma::zeros(num_elastic_coef, 1);
        for (int ii = 0; ii < num_elastic_coef; ii++)
            elastic_coef(ii, 0) = ecNode[ii].as<double>();

        if (elastic_coef(0, 0) > 0.0)
            EA = elastic_coef(0, 0);
        else
        {
            double strain_EA = 0.1;
            EA = elastic_coef(0, 0) + 2 * elastic_coef(1, 0) * strain_EA +
                 3 * elastic_coef(2, 0) * strain_EA * strain_EA;
        }
    }
    else if (flag_stiffness > 1)
    {
        YAML::Node strainNode = node["strain_data"];
        YAML::Node stressNode = node["stress_data"];
        strain_data = arma::zeros(flag_stiffness, 1);
        stress_data = arma::zeros(flag_stiffness, 1);
        for (int ii = 0; ii < flag_stiffness; ii++)
            strain_data(ii, 0) = strainNode[ii].as<double>();
        for (int ii = 0; ii < flag_stiffness; ii++)
            stress_data(ii, 0) = stressNode[ii].as<double>();
        beta = node["beta"].as<double>();

        EA = 0;
        int i_EA = 0;
        while (EA <= 0)
        {
            EA = (stress_data(i_EA + 1, 0) - stress_data(i_EA, 0)) / (strain_data(i_EA + 1, 0) - strain_data(i_EA, 0));
            i_EA++;
            if (i_EA > flag_stiffness - 1)
            {
                std::stringstream ss;
                ss << "Strain-Stress curve should have positive slope for line type.\n";
                throw ValueError(ss.str());
            }
        }
    }

    CB = node["CB"].as<double>();
    Cmn = node["Cmn"].as<double>();
    Cdn = node["Cdn"].as<double>();
    Cdt = node["Cdt"].as<double>();
    GK = node["GK"].as<double>();
    GC = node["GC"].as<double>();
    smoothstep = node["smoothstep"].as<int>();
    frictionModel = node["friction_model"].as<int>();
    vth = node["vth"].as<double>();
    ust = node["ust"].as<double>();
    usn = node["usn"].as<double>();
    ud = node["ud"].as<double>();
    deltamax = node["deltamax"].as<double>();

    if (frictionModel == 1)
    {
        arma::mat tmpA = arma::zeros(4, 4);
        arma::mat tmpB = arma::zeros(4, 1);
        double x1 = 0;
        double x2 = 1e-04;
        tmpA = {{pow(x1, 3), pow(x1, 2), 1 * x1, 1},
                {3 * pow(x1, 2), 2 * (x1), 1, 0},
                {pow(x2, 3), pow(x2, 2), 1 * x2, 1},
                {3 * pow(x2, 2), 2 * x2, 1, 0}};
        tmpB(0, 0) = 1e-07;
        tmpB(1, 0) = 0.0;
        tmpB(2, 0) = x2;
        tmpB(3, 0) = 1;
        a_1 = arma::solve(tmpA, tmpB);
    }
}
