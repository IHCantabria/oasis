
#include <armadillo>
#include <algorithm>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include "OWC.hpp"
#include "../MathTools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"

// OWC Turbine type methods ------------------------------------------------------

OWCTurbineType::OWCTurbineType(int n, Simulation *pSimInp)
{
    idOWCTurbineType = n + 1;
    pSim = pSimInp;
}

void OWCTurbineType::ReadPropertiesASCII(FILE *pFile)
{
    // Declare variables
    char buffer_line[1000];
    char cBuckCurvesFileName[1000];
    // Ignore the first three lines, where it says "OWC Turbine Type #"
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), pFile);
    }
    // Read and assign the properties of the turbine type
    fscanf(pFile, "%lf %[^\n]\n", &inertia, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &damping, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &rotor_diameter, buffer_line);
    fscanf(pFile, "%d %[^\n]\n", &valve_type, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &bypass_valve_area, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &bypass_valve_discharge_coef, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &throttle_valve_area, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &throttle_valve_discharge_coef, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &throttle_gap_air_volume, buffer_line);
    fscanf(pFile, "%d %[^\n]\n", &controller_type, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &generator_rated_power, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &generator_max_torque, buffer_line);
    fscanf(pFile, "%s %[^\n]\n", cBuckCurvesFileName, buffer_line);

    // Get the file name of the buckingham curves
    std::string filename = cBuckCurvesFileName;
    std::string filepath = JoinPath(this->pSim->inputFolderPath, filename);

    // Load the buckingham curves from the file, a CSV file with columns:
    // adim_pressure, adim_mass_flow_rate, efficiency and adim_power
    std::cout << "  --> Reading Buckingham curves from file: " << filepath << std::endl;

    arma::field<std::string> header(4);
    header(0) = "adim_pressure";
    header(1) = "adim_mass_flow_rate";
    header(2) = "efficiency";
    header(3) = "adim_power";
    arma::mat buck_curve_data;
    buck_curve_data.load(arma::csv_name(filepath, header));
    buck_curve_adim_pressure = buck_curve_data.col(0);
    buck_curve_adim_mass_flow_rate = buck_curve_data.col(1);
    buck_curve_efficiency = buck_curve_data.col(2);
    buck_curve_adim_power = buck_curve_data.col(3);

    std::cout << "  --> ...done! Number of points in curves: " << buck_curve_adim_pressure.n_elem << std::endl;
}

void OWCTurbineType::Initialize(FILE *pFile)
{
    // Read properties from the file
    ReadPropertiesASCII(pFile);
    // Get the maximum efficiency of the turbine
    arma::uword max_efficiency_index = buck_curve_efficiency.index_max();
    // Get the best efficiency point (BEP) of the turbine for the adimensional pressure and power
    adim_pressure_bep = buck_curve_adim_pressure(max_efficiency_index);
    adim_power_bep = buck_curve_adim_power(max_efficiency_index);
    // Get the maximum power of the turbine
    arma::uword max_power_index = buck_curve_adim_power.index_max();
    // Get the adimensional critical pressure and mass flow rate of the turbine
    adim_critical_pressure = buck_curve_adim_pressure(max_power_index);
    adim_critical_mass_flow_rate = buck_curve_adim_mass_flow_rate(max_power_index);
    // Set the algebraic controller parameters
    param_a = pSim->airAtmPres * pow(rotor_diameter, 5) * adim_power_bep;
    param_b = 2.0;
}

void OWCTurbineType::Finalize(void)
{
    // Free memory
    delete[] buck_curve_adim_pressure.memptr();
    delete[] buck_curve_adim_power.memptr();
    delete[] buck_curve_adim_mass_flow_rate.memptr();
    delete[] buck_curve_efficiency.memptr();
}

// OWC Turbine methods ------------------------------------------------------

OWCTurbine::OWCTurbine(Simulation *pSimInp, int id_OWCTurbineType)
{
    pSim = pSimInp;
    pOWCTurbineType = pSim->pOWCTurbines[id_OWCTurbineType - 1];
}

void OWCTurbine::Initialize(double initial_angular_velocity)
{
    mass_flow_rate = 0.0;
    angular_velocity = initial_angular_velocity;
    angular_acceleration = 0.0;
    air_torque = 0.0;
    generator_torque = 0.0;
    power = 0.0;
    critical_pressure = 0.0;
    critical_mass_flow_rate = 0.0;
    gap_rel_pressure = 1.0;
    gap_rel_pressure_dot = 0.0;
    gap_mass_flow_rate = 0.0;
    gap_air_pressure = pSim->airAtmPres;
    gap_air_density = pSim->airAtmPresDensity;
    gap_air_mass = pOWCTurbineType->throttle_gap_air_volume * gap_air_density;
}

bool OWCTurbine::IsAdimPressureInRange(double adim_pressure)
{
    if (adim_pressure < pOWCTurbineType->buck_curve_adim_pressure.min())
    {
        return false;
    }
    if (adim_pressure > pOWCTurbineType->buck_curve_adim_pressure.max())
    {
        return false;
    }
    return true;
}

void OWCTurbine::ComputeMassFlowRate(double pressure_diff, double stagnation_density)
{
    double adim_pressure = abs(pressure_diff) / (stagnation_density * pow(pOWCTurbineType->rotor_diameter, 2) * pow(angular_velocity, 2));
    if (!IsAdimPressureInRange(adim_pressure))
    {
        std::cout << "WARNING: The adimensional pressure is out of range. The turbine may not work properly." << std::endl;
        adim_pressure = std::clamp(adim_pressure, pOWCTurbineType->buck_curve_adim_pressure.min(), pOWCTurbineType->buck_curve_adim_pressure.max());
    }
    arma::vec adim_pressure_vec = adim_pressure * arma::ones<arma::vec>(1);
    arma::vec adim_mass_flow_rate_vec;
    arma::interp1(pOWCTurbineType->buck_curve_adim_pressure, pOWCTurbineType->buck_curve_adim_mass_flow_rate, adim_pressure_vec, adim_mass_flow_rate_vec);
    double adim_mass_flow_rate = arma::as_scalar(adim_mass_flow_rate_vec(0));
    mass_flow_rate = adim_mass_flow_rate * stagnation_density * pow(pOWCTurbineType->rotor_diameter, 3) * angular_velocity * arma::sign(pressure_diff);
}

void OWCTurbine::ComputeAirTorque(double pressure_diff, double stagnation_density)
{
    double adim_pressure = abs(pressure_diff) / (stagnation_density * pow(pOWCTurbineType->rotor_diameter, 2) * pow(angular_velocity, 2));
    if (!IsAdimPressureInRange(adim_pressure))
    {
        std::cout << "WARNING: The adimensional pressure is out of range. The turbine may not work properly." << std::endl;
        adim_pressure = std::clamp(adim_pressure, pOWCTurbineType->buck_curve_adim_pressure.min(), pOWCTurbineType->buck_curve_adim_pressure.max());
    }
    arma::vec adim_pressure_vec = adim_pressure * arma::ones<arma::vec>(1);
    arma::vec adim_power_vec;
    arma::interp1(pOWCTurbineType->buck_curve_adim_pressure, pOWCTurbineType->buck_curve_efficiency, adim_pressure_vec, adim_power_vec);
    double adim_power = arma::as_scalar(adim_power_vec(0));
    power = adim_power * stagnation_density * pow(pOWCTurbineType->rotor_diameter, 5) * pow(angular_velocity, 3);
    air_torque = power / angular_velocity;
}

void OWCTurbine::ComputeCriticalPressure(double stagnation_density)
{
    critical_pressure = pOWCTurbineType->adim_critical_pressure * stagnation_density * pow(pOWCTurbineType->rotor_diameter, 2) * pow(angular_velocity, 2);
}

void OWCTurbine::ComputeCriticalMassFlowRate(double stagnation_density)
{
    critical_mass_flow_rate = pOWCTurbineType->adim_critical_mass_flow_rate * stagnation_density * pow(pOWCTurbineType->rotor_diameter, 3) * angular_velocity;
}

void OWCTurbine::ComputeAngularAcceleration(double pressure_diff, double stagnation_density)
{
    if ((0.5 * angular_velocity * pOWCTurbineType->rotor_diameter) >= 180.0)
    {
        std::cout << "WARNING: The turbine is rotating too fast. Shock waves may occur." << std::endl;
    }
    ComputeAirTorque(pressure_diff, stagnation_density);
    angular_acceleration = (air_torque - generator_torque) / pOWCTurbineType->inertia;
}

void OWCTurbine::ComputeGenTorque(void)
{
    if (pOWCTurbineType->controller_type == 0)
    {
        generator_torque = 0.0;
    }
    else if (pOWCTurbineType->controller_type == 1)
    {
        generator_torque = pOWCTurbineType->param_a * pow(angular_velocity, pOWCTurbineType->param_b);
        generator_torque = std::min(generator_torque, pOWCTurbineType->generator_max_torque);
        generator_torque = std::min(generator_torque, pOWCTurbineType->generator_rated_power / angular_velocity);
    }
    else
    {
        // TODO: Implement a functional controller that avoids stall with negative generator torque
        // and uses the valves to control the pressure.
        std::stringstream ss;
        ss << "ERROR: The controller type " << pOWCTurbineType->controller_type << " is not implemented." << std::endl;
        throw NotImplementedError(ss.str());
    }
}

// OWC methods ------------------------------------------------------

OWC::OWC(int n, Simulation *pSimInp)
{
    idOWC = n + 1;
    pSim = pSimInp;
}

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
    fscanf(pFile, "%lf %[^\n]\n", &turb_initial_angular_velocity, buffer_line);
    double dtemp;
    for (int ii = 0; ii < 3; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the local position of the OWC: " << idOWC + 1 << "\n";
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
    rel_pressure = 1.0;
    rel_pressure_dot = 0.0;
    air_pressure = pSim->airAtmPres;
    air_volume = reference_air_volume;
    air_volume_dot = 0.0;
    air_density = pSim->airAtmPresDensity;
    air_mass = air_density * air_volume;
    air_mass_dot = 0.0;
    gap_pressure_diff = 0.0;
    gap_rel_pressure = 1.0;
    gap_rel_pressure_dot = 0.0;

    // Create the OWC turbine if it is defined
    if (turbine_type > 0)
    {
        if (pSim->numOWCTurbines < turbine_type)
        {
            std::stringstream ss;
            ss << "OWC turbine type: " << turbine_type << " in OWC " << idOWC + 1 << " was not declared in dataOWCTurbines.dat\n";
            throw ValueError(ss.str());
        }
        pOWCTurbine = new OWCTurbine(pSim, turbine_type);
        // TODO: set this in a init function
        pOWCTurbine->Initialize(turb_initial_angular_velocity);
    }
    else
    {
        pOWCTurbine = NULL;
    }

    // Open output file
    char buffer[50];
    int nn = sprintf(buffer, "OWC_%d.txt", idOWC);
    std::string file_path = JoinPath(pSim->outputFolderPath, buffer);
    pfile_output = fopen(file_path.c_str(), "w");
    if (pfile_output == NULL)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to open the output file for OWC: " << idOWC + 1 << "\n";
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
        force_pneumatic = pSim->airAtmPres * (rel_pressure - 1.0) * waterplane_area;
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
        int turb_valve_type = pOWCTurbine->pOWCTurbineType->valve_type;

        double turb_pressure_diff;
        double turb_stagnation_air_density;

        if (turb_valve_type == 0) // Case with no valves
        {
            turb_pressure_diff = pSim->airAtmPres - air_pressure;
            turb_stagnation_air_density = std::max(air_density, pSim->airAtmPresDensity);
            pOWCTurbine->ComputeMassFlowRate(turb_pressure_diff, turb_stagnation_air_density);
            air_mass_dot = pOWCTurbine->mass_flow_rate;
        }
        else if (turb_valve_type == 1) // Case with a bypass valve
        {
            turb_pressure_diff = pSim->airAtmPres - air_pressure;
            turb_stagnation_air_density = std::max(air_density, pSim->airAtmPresDensity);
            pOWCTurbine->ComputeMassFlowRate(turb_pressure_diff, turb_stagnation_air_density);
            double valve_mass_flow_rate = ComputeHoleMassFlowRate(turb_pressure_diff,
                                                                  turb_stagnation_air_density,
                                                                  pOWCTurbine->pOWCTurbineType->bypass_valve_area,
                                                                  pOWCTurbine->pOWCTurbineType->bypass_valve_discharge_coef);
            air_mass_dot = pOWCTurbine->mass_flow_rate + valve_mass_flow_rate;
        }
        else if (turb_valve_type == 2) // Case with a throttle valve
        {
            pOWCTurbine->gap_air_pressure = pSim->airAtmPres * pOWCTurbine->gap_rel_pressure;
            pOWCTurbine->gap_air_density = ComputeAirDensity(pOWCTurbine->gap_air_pressure);
            pOWCTurbine->gap_air_mass = pOWCTurbine->gap_air_density * pOWCTurbine->pOWCTurbineType->throttle_gap_air_volume;

            turb_pressure_diff = pSim->airAtmPres - pOWCTurbine->gap_air_pressure;
            turb_stagnation_air_density = std::max(pOWCTurbine->gap_air_density, pSim->airAtmPresDensity);
            pOWCTurbine->ComputeMassFlowRate(turb_pressure_diff, turb_stagnation_air_density);

            double valve_pressure_diff = pOWCTurbine->gap_air_pressure - air_pressure;
            double valve_stagnation_air_density = std::max(air_density, pOWCTurbine->gap_air_density);
            double valve_mass_flow_rate = ComputeHoleMassFlowRate(valve_pressure_diff,
                                                                  valve_stagnation_air_density,
                                                                  pOWCTurbine->pOWCTurbineType->throttle_valve_area,
                                                                  pOWCTurbine->pOWCTurbineType->throttle_valve_discharge_coef);

            air_mass_dot = valve_mass_flow_rate;
            pOWCTurbine->gap_mass_flow_rate = pOWCTurbine->mass_flow_rate - valve_mass_flow_rate;
        }
        else if (turb_valve_type == 3) // Case with both valves
        {
            pOWCTurbine->gap_air_pressure = pSim->airAtmPres * pOWCTurbine->gap_rel_pressure;
            pOWCTurbine->gap_air_density = ComputeAirDensity(pOWCTurbine->gap_air_pressure);
            pOWCTurbine->gap_air_mass = pOWCTurbine->gap_air_density * pOWCTurbine->pOWCTurbineType->throttle_gap_air_volume;

            turb_pressure_diff = pSim->airAtmPres - pOWCTurbine->gap_air_pressure;
            turb_stagnation_air_density = std::max(pOWCTurbine->gap_air_density, pSim->airAtmPresDensity);
            pOWCTurbine->ComputeMassFlowRate(turb_pressure_diff, turb_stagnation_air_density);

            double throttle_valve_pressure_diff = pOWCTurbine->gap_air_pressure - air_pressure;
            double throttle_valve_stagnation_air_density = std::max(air_density, pOWCTurbine->gap_air_density);
            double throttle_valve_mass_flow_rate = ComputeHoleMassFlowRate(throttle_valve_pressure_diff,
                                                                           throttle_valve_stagnation_air_density,
                                                                           pOWCTurbine->pOWCTurbineType->throttle_valve_area,
                                                                           pOWCTurbine->pOWCTurbineType->throttle_valve_discharge_coef);
            double bypass_valve_pressure_diff = pSim->airAtmPres - air_pressure;
            double bypass_valve_stagnation_air_density = std::max(air_density, pSim->airAtmPresDensity);
            double bypass_valve_mass_flow_rate = ComputeHoleMassFlowRate(bypass_valve_pressure_diff,
                                                                         bypass_valve_stagnation_air_density,
                                                                         pOWCTurbine->pOWCTurbineType->bypass_valve_area,
                                                                         pOWCTurbine->pOWCTurbineType->bypass_valve_discharge_coef);

            air_mass_dot = throttle_valve_mass_flow_rate + bypass_valve_mass_flow_rate;
            pOWCTurbine->gap_mass_flow_rate = pOWCTurbine->mass_flow_rate - throttle_valve_mass_flow_rate;
        }
        else
        {
            std::stringstream ss;
            ss << "In OWC: " << idOWC + 1 << ", valve type must be 0-3, but " << turb_valve_type << " was provided.\n";
            throw ValueError(ss.str());
        }

        // Compute the relative pressure rate in the gap
        if (turb_valve_type == 2 || turb_valve_type == 3)
        {
            pOWCTurbine->gap_rel_pressure_dot = ComputeRelPressureRate(pOWCTurbine->gap_air_pressure,
                                                                       pOWCTurbine->gap_mass_flow_rate,
                                                                       pOWCTurbine->gap_air_mass,
                                                                       0.0,
                                                                       pOWCTurbine->pOWCTurbineType->throttle_gap_air_volume);
        }

        // Compute the turbine angular acceleration
        if (pOWCTurbine->pOWCTurbineType->controller_type > 0)
        {
            pOWCTurbine->ComputeAngularAcceleration(turb_pressure_diff, turb_stagnation_air_density);
        }
        else
        {
            pOWCTurbine->angular_acceleration = 0.0;
        }
    }

    // Compute the relative pressure rate
    rel_pressure_dot = ComputeRelPressureRate(air_pressure, air_mass_dot, air_mass, air_volume_dot, air_volume);
}

double OWC::ComputeAirDensity(double pressure)
{
    if (pressure < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the air density in OWC: " << idOWC + 1 << "\n";
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
        ss << "An error occurred when trying to compute the relative pressure rate in OWC: " << idOWC + 1 << "\n";
        throw ValueError(ss.str());
    }
    if (stagnation_density < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the stagnation density in OWC: " << idOWC + 1 << "\n";
        throw ValueError(ss.str());
    }
    // Compute the hole mass flow rate using the orifice equation
    // Regularize sqrt(|Dp|)*sign(Dp) -> Dp/sqrt(|Dp|+eps) to avoid infinite derivative at Dp=0
    // which causes numerical Jacobian inaccuracy in the ESDIRK solver
    double eps_dp = 1.0; // [Pa] small regularization parameter
    double hole_mass_flow_rate = hole_discharge_coef * hole_area * sqrt(2.0 * stagnation_density) * pressure_diff / sqrt(abs(pressure_diff) + eps_dp);
    return hole_mass_flow_rate;
}

double OWC::ComputeRelPressureRate(double air_pressure, double mass_flow_rate, double air_mass, double air_volume_dot, double air_volume)
{
    if (rel_pressure < 0.0)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to compute the relative pressure rate in OWC: " << idOWC + 1 << "\n";
        throw ValueError(ss.str());
    }
    // Compute the relative pressure rate using the orifice equation
    return pSim->airAdiabaticDilation * (mass_flow_rate / air_mass - air_volume_dot / air_volume) * air_pressure / pSim->airAtmPres;
}