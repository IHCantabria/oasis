
#include <armadillo>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include "OWC.hpp"
#include "../MathTools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"

void OWC::ReadPropertiesASCII(FILE *pFile)
{
    // Declare variables
    char buffer_line[1000];

    // Ignore th first three lines, where it says "OWC #"
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), pFile);
    }

    // Read and assign the body corresponding to the OWC
    fscanf(pFile, "%d %[^\n]\n", &idBodyOWC, buffer_line);
    fscanf(pFile, "%d %[^\n]\n", &idBodyFloater, buffer_line);
    idBodyOWC--;
    idBodyFloater--;
    pBodyOWC = pSim->pBodies[idBodyOWC];
    pBodyFloater = pSim->pBodies[idBodyFloater];

    // Check if the OWC body has a single degree of freedom
    if (pBodyOWC->numDofs != 1)
    {
        std::stringstream ss;
        ss << "The OWC body must have a single degree of freedom. The body " << idBodyOWC + 1 << " has " << pBodyOWC->numDofs << " degrees of freedom.\n";
        throw ValueError(ss.str());
    }
    // Check if the single degree of freedom is vertical
    if (pBodyOWC->pDofs[0] != 2)
    {
        std::stringstream ss;
        ss << "The OWC body must have a single degree of freedom in the vertical direction. The body " << idBodyOWC + 1 << " has the degree of freedom " << pBodyOWC->pDofs[0] + 1 << ".\n";
        throw ValueError(ss.str());
    }
    // Check that the OWC body has no BCPs
    if (pBodyOWC->numBcps != 0)
    {
        std::stringstream ss;
        ss << "The OWC body must not have any BCPs. The body " << idBodyOWC + 1 << " has " << pBodyOWC->numBcps << " BCPs.\n";
        throw ValueError(ss.str());
    }
    // Check that the OWC body has no wind turbines
    if (pBodyOWC->numWindTurbs != 0)
    {
        std::stringstream ss;
        ss << "The OWC body must not have any wind turbines. The body " << idBodyOWC + 1 << " has " << pBodyOWC->numWindTurbs << " wind turbines.\n";
        throw ValueError(ss.str());
    }
    // Check that the OWC body has linear hydrostatics
    if (pBodyOWC->flag_hydrostatics != 0)
    {
        std::stringstream ss;
        ss << "The OWC body must have linear hydrostatics. The body " << idBodyOWC + 1 << " has flag_hydrostatics = " << pBodyOWC->flag_hydrostatics << ".\n";
        throw ValueError(ss.str());
    }

    fscanf(pFile, "%lf %[^\n]\n", &waterplane_area, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &reference_air_volume, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &hole_area, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &hole_discharge_coef, buffer_line);
    fscanf(pFile, "%d %[^\n]\n", &turbine_type, buffer_line);
    double dtemp;
    for (int ii = 0; ii < 3; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the local position of the OWC: " << idOWC << "\n";
            throw ValueError(ss.str());
        }
        pos_local(ii, 0) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
}

void OWC::Initialize(FILE *pFile)
{
    ReadPropertiesASCII(pFile);

    // Initialize variables
    displacement = pBodyOWC->pos(2, 0); // TODO: Review this. Ref frames?
    velocity = 0.0;
    floater_displacement = pBodyFloater->pos(2, 0); // TODO: Review this. Ref frames?
    floater_velocity = 0.0;
    rel_pressure = 0.0;
    rel_pressure_dot = 0.0;
    air_pressure = pSim->airAtmPres;

    // Open output file
    char buffer[50];
    int nn = sprintf(buffer, "OWC_%d.txt", idOWC);
    std::string file_path = JoinPath(pSim->outputFolderPath, buffer);
    pfile_output = fopen(file_path.c_str(), "w");
    if (pfile_output == NULL)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to open the output file for OWC: " << idOWC << "\n";
        throw IOError(ss.str());
    }
    // Write header
    fprintf(pfile_output, "Time    Displacement    Velocity    Pressure\n");
    fprintf(pfile_output, "%f    %f    %f    %f\n", 0.0, displacement, velocity, air_pressure);
    // fflush(pfile_output);
}

void OWC::WriteOut(double time)
{
    // Write output to file
    fprintf(pfile_output, "%f    %f    %f    %f\n", time, displacement, velocity, air_pressure);
    // fflush(pfile_output);
}

void OWC::Finalize(void)
{
    // Close output file
    fclose(pfile_output);
}

void OWC::ComputeForces(double time)
{
    // Get the position and velocity of the OWC body
    // TODO: Review this. Ref frames?
    displacement = arma::as_scalar(pBodyOWC->pos(2, 0));
    velocity = arma::as_scalar(pBodyOWC->vel(2, 0));

    // Get the position and velocity of the floater body
    // TODO: Review this. Ref frames?
    floater_displacement = arma::as_scalar(pBodyFloater->pos(2, 0));
    floater_velocity = arma::as_scalar(pBodyFloater->vel(2, 0));

    // Compute the pneumatic forces from the pressure difference
    double force_pneumatic = 0.0;
    if (turbine_type >= 0)
    {
        ComputePressure(time);
        force_pneumatic = rel_pressure * waterplane_area;
    }

    // Save the forces in the bodies
    // For the OWC body, simply add the pneumatic force in the vertical direction
    pBodyOWC->owcForces(2, 0) = -force_pneumatic;
    // For the floater body
    // get the rotation matrix and add the pneumatic force in the local direction in global coordinates
    // and the moment in local coordinates
    arma::mat tmpForceLoc = arma::zeros(3, 1);
    tmpForceLoc(2, 0) = force_pneumatic;
    pBodyFloater->owcForces.rows(0, 2) = pBodyFloater->owcForces.rows(0, 2) + pBodyFloater->rotMat * tmpForceLoc;
    pBodyFloater->owcForces.rows(3, 5) = pBodyFloater->owcForces.rows(3, 5) + arma::cross(pos_local, tmpForceLoc);
}

void OWC::ComputePressure(double time)
{
    air_volume = reference_air_volume + waterplane_area * (floater_displacement - displacement);
    air_volume_dot = waterplane_area * (floater_velocity - velocity);
    air_pressure = pSim->airAtmPres * rel_pressure;
    air_density = ComputeAirDensity(air_pressure);
    air_mass = air_density * air_volume;

    if (turbine_type == 0)
    {
        double hole_pressure_diff = pSim->airAtmPres - air_pressure;
        double hole_stagnation_air_density = std::max(air_density, pSim->airAtmPresDensity);

        air_mass_dot = ComputeHoleMassFlowRate(hole_pressure_diff, hole_stagnation_air_density, hole_area, hole_discharge_coef);
    }
    else if (turbine_type > 0)
    {
        std::stringstream ss;
        ss << "Turbines are not implemented yet for OWCs! \n";
        throw ValueError(ss.str());
    }

    // Compute the relative pressure rate
    rel_pressure_dot = ComputeRelPressureRate(air_pressure, air_mass_dot, air_mass, air_volume_dot, air_volume);
}

double OWC::ComputeAirDensity(double pressure)
{
    if (pressure < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the air density in OWC: " << idOWC << "\n";
        throw ValueError(ss.str());
    }

    // Compute the air density using the ideal gas law
    double air_density = pSim->airAtmPresDensity * pow(pressure / pSim->airAtmPres, 1.0 / pSim->airAdiabaticDilation);
    return air_density;
}

double OWC::ComputeHoleMassFlowRate(double pressure_diff, double stagnation_density, double hole_area, double hole_discharge_coef)
{
    if (rel_pressure < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the relative pressure rate in OWC: " << idOWC << "\n";
        throw ValueError(ss.str());
    }
    if (stagnation_density < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the stagnation density in OWC: " << idOWC << "\n";
        throw ValueError(ss.str());
    }
    // Compute the hole mass flow rate using the orifice equation
    double hole_mass_flow_rate = hole_discharge_coef * hole_area * sqrt(2.0 * abs(pressure_diff) * stagnation_density) * arma::sign(pressure_diff);
    return hole_mass_flow_rate;
}

double OWC::ComputeRelPressureRate(double air_pressure, double mass_flow_rate, double air_mass, double air_volume_dot, double air_volume)
{
    if (rel_pressure < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the relative pressure rate in OWC: " << idOWC << "\n";
        throw ValueError(ss.str());
    }
    // Compute the relative pressure rate using the orifice equation
    double rel_pressure_rate = pSim->airAdiabaticDilation * (mass_flow_rate / air_mass - air_volume_dot / air_volume) * air_pressure / pSim->airAtmPres;
    return rel_pressure_rate;
}