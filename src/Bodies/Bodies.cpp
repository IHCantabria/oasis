
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdio>
#include <cmath>
#include <armadillo>
#include "../Simulations/Simulation.hpp"
#include "Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../os_tools.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../Logger.hpp"

//
Body::Body(int n, Simulation* pIncSim)
{
    id = n;
    pSim = pIncSim;
}

// Computes the resultant force and moment at the COG from all BCP forces
void Body::ComputeBcpForces(void)
{

    // Local working variables
    arma::mat posG_temp;
    arma::mat ForceBCP_temp;
    arma::mat F_M;
    arma::mat M_F;
    double M_norm, det;
    double rx, ry, rz, Mx, My, Mz;

    for (int ii = 0; ii < numBcps; ii++)
    {

        posG_temp = pBodyBcps[ii]->posWrtCdgGlobal; // Store the global position of the BCP in a local variable
        ForceBCP_temp = pBodyBcps[ii]->forceBcp;    // Store the forces and moments at the BCP in a local variable

        M_F = arma::cross(posG_temp,
                          ForceBCP_temp.rows(0, 2)); // Moment about the COG caused by the BCP force, in global frame

        bcpForces.rows(0, 2) =
            bcpForces.rows(0, 2) + ForceBCP_temp.rows(0, 2); // Accumulate total force at the COG in global frame.
        bcpForces.rows(3, 5) =
            bcpForces.rows(3, 5) +
            rotMat.t() * (ForceBCP_temp.rows(3, 5) + M_F); // Accumulate total moment at the COG in local frame.
    }

    if (bcpForces.has_nan())
    {
        Logger::error("NaN Detected on body with id = " + std::to_string(id));
        throw std::exception();
    }
}

void Body::ComputeWindTurbForces(void)
{

    windTurbForces = arma::zeros(6, 1);
    for (int ii = 0; ii < numWindTurbs; ii++)
    {
        windTurbForces += pBodyWindTurbs[ii]->forceBodyCOG;
    }
}

int Body::GetId(void)
{
    return id;
}

void Body::LoadHydrodynamicDatabase(Body** hydroDatabaseBodies, int slotIndex)
{
    Logger::info("--> Reading Hydrodynamics Properties (HDF5 format)");

    // File path
    std::string file_path = JoinPath(this->pSim->inputFolderPath, this->hydroDatabaseName);

    // Load hydrodynamic database
    // slotIndex is the position of this body in the hydroDatabaseBodies array
    this->pHydro = new HydroDatabase(this->hydroDatabaseIndex, slotIndex, hydroDatabaseBodies, this->pSim);
    Logger::info("    --> Loading hydrodynamic database");
    this->pHydro->LoadHydrodynamicData(file_path);

    Logger::info("    --> Loading and checking the Hydrodynamic C.O.G position");
    // Load and check the Hydrodynamic C.O.G position
    if (this->takeCOGHydroDatabase == 1)
    {
        // Load COG equilibrium position
        this->pos_eq.rows(0, 2) = this->pHydro->GetCog().t();

        std::ostringstream oss_hdb;
        oss_hdb << "      --> Using HDB COG: " << this->pHydro->GetCog();
        Logger::debug(oss_hdb.str());

        // Add initial displacement to the COG equilibrium position
        this->pos = this->pos + this->pos_eq;
        pos_ini = pos;

        std::ostringstream oss_pos;
        oss_pos << "      --> Initial Position: " << this->pos.t();
        Logger::debug(oss_pos.str());
    }
    else if ((this->takeCOGHydroDatabase == 0) && (this->pHydro->GetNumBodies() > 1))
    {
        // Declare local variables
        bool xcond, ycond, zcond;

        // Check if the input C.O.G is in accordance with the hydrodynamic database
        double cog_tol = 1e-6;
        arma::mat cog = this->pHydro->GetCog();
        xcond = fabs(this->pos_eq(0, 0) - cog(0, 0)) > cog_tol;
        ycond = fabs(this->pos_eq(1, 0) - cog(0, 1)) > cog_tol;
        zcond = fabs(this->pos_eq(2, 0) - cog(0, 2)) > cog_tol;

        if (xcond || ycond || zcond)
        {
            std::stringstream ss;
            ss << "Specified Center of Gravity for body: " << this->GetId();
            ss << " mismatch with the C.O.G value in the Hydrodynamic database.\n";
            throw ValueError(ss.str());
        }
    }

    Logger::info("----> Hydrodynamic Properties Read");
}

// Read body properties from the input file
void Body::ReadPropertiesASCII(FILE* pFile)
{
    // Declare variables
    char buffer_line[1000];
    fpos_t carriage_init;
    char cHydroDatabaseName[1000];
    char cMovementsFileName[1000];
    char cHydrostaticMeshName[1000];
    double dtemp;
    int itemp;

    Logger::info("--> Reading Body: " + std::to_string(this->GetId() + 1));

    // Read flag to take COG from Hydrodynamic database
    if (fscanf(pFile, "%d %[^\n]\n", &takeCOGHydroDatabase, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag to take COG from Hydrodynamic database."
           << ".\n";
        throw ValueError(ss.str());
    }

    Logger::debug(" takeCOGHydroDatabase: " + std::to_string(takeCOGHydroDatabase));

    // Read dofs considered
    fgetpos(pFile, &carriage_init);
    while (fscanf(pFile, "%d", &itemp) == 1)
    {
        this->numDofs++;
    }
    fsetpos(pFile, &carriage_init);

    Logger::debug(" numDofs: " + std::to_string(numDofs));

    this->pDofs = new int[this->numDofs];
    for (int ii = 0; ii < this->numDofs; ii++)
    {
        if (fscanf(pFile, "%d", &itemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the degrees of freedom for body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        this->pDofs[ii] = itemp - 1;
        isDofActive(itemp - 1, 0) = 1.0;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    Logger::debug("    --> Degrees of Freedom: " + std::to_string(numDofs) + " dofs read");

    // Read the boundary condition points in the body
    fgetpos(pFile, &carriage_init);
    while (fscanf(pFile, "%d", &itemp) == 1)
    {
        this->numBcps++;
    }
    fsetpos(pFile, &carriage_init);

    this->pIndexBcps = new int[this->numBcps];
    fgetpos(pFile, &carriage_init);
    fgets(buffer_line, sizeof(buffer_line), pFile);
    fsetpos(pFile, &carriage_init);
    for (int ii = 0; ii < this->numBcps; ii++)
    {
        if (fscanf(pFile, "%d", &itemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the boundary condition points of the body: " << this->GetId()
               << " "
               << "\n";
            throw ValueError(ss.str());
        }
        this->pIndexBcps[ii] = itemp - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    Logger::debug("    --> Boundary Condition Points: " + std::to_string(numBcps) + " BCPs read");

    // Read the wind turbines in the body
    fgetpos(pFile, &carriage_init);
    while (fscanf(pFile, "%d", &itemp) == 1)
    {
        this->numWindTurbs++;
    }
    fsetpos(pFile, &carriage_init);

    this->pBodyWindTurbs = new WindTurbine*[this->numWindTurbs];
    this->pIndexWindTurbs = new int[this->numWindTurbs];
    fgetpos(pFile, &carriage_init);
    fgets(buffer_line, sizeof(buffer_line), pFile);
    fsetpos(pFile, &carriage_init);
    for (int ii = 0; ii < this->numWindTurbs; ii++)
    {
        if (fscanf(pFile, "%d", &itemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the wind turbines of the body: " << this->GetId() << " "
               << "\n";
            throw ValueError(ss.str());
        }
        this->pIndexWindTurbs[ii] = itemp - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    Logger::debug("    --> Wind Turbines: " + std::to_string(numWindTurbs) + " wind turbines read");

    // Read COG equilibrium position
    for (int ii = 0; ii < 6; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the initial position of the body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }

        if (this->takeCOGHydroDatabase == 0)
        {
            this->pos_eq(ii, 0) = dtemp;
            this->pos(ii, 0) = dtemp;
        }
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    if (this->takeCOGHydroDatabase == 0)
    {
        std::ostringstream oss_eq;
        oss_eq << "    --> Equilibrium Position: " << this->pos.t();
        Logger::debug(oss_eq.str());
    }

    // Read Initial displacement from reference position
    for (int ii = 0; ii < 6; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the initial displacement of the body: " << this->GetId()
               << "\n";
            throw ValueError(ss.str());
        }
        pos(ii, 0) += dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    pos_ini = pos;

    std::ostringstream oss_disp;
    oss_disp << "    --> Initial Disp: " << pos.t();
    Logger::debug(oss_disp.str());

    // Read hydrodynamic database filename
    if (fscanf(pFile, "%s %[^\n]\n", cHydroDatabaseName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the hydrodynamic database name of the body: " << this->GetId()
           << "\n";
        throw ValueError(ss.str());
    }
    this->hydroDatabaseName = cHydroDatabaseName;

    Logger::debug("    --> hydrodynamic Database: " + this->hydroDatabaseName);

    // Read body index in the associated database
    if (fscanf(pFile, "%d %[^\n]\n", &(this->hydroDatabaseIndex), buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the index of the body in the hydrodynamic database: "
           << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->hydroDatabaseIndex--;

    Logger::debug("    --> hydrodynamic Database Index: " + std::to_string(this->hydroDatabaseIndex + 1));

    // Read flag for blocking the body
    if (fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for blocking body."
           << ".\n";
        throw ValueError(ss.str());
    }
    if (itemp < 3)
    {
        flag_blocked = itemp;
    }
    else
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Flag for blocking body not available, must be 0, 1 or 2."
           << ".\n";
        throw ValueError(ss.str());
    }

    Logger::debug("    --> Flag Blocked: " + std::to_string(flag_blocked));

    // Read movements filename
    if (fscanf(pFile, "%s %[^\n]\n", cMovementsFileName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the movements file name of the body: " << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->movementsFileName = cMovementsFileName;

    Logger::debug("    --> Movements File Name: " + this->movementsFileName);

    // Read flag for hydrostatics of the body
    if (fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for hydrostatics of the body."
           << ".\n";
        throw ValueError(ss.str());
    }
    if (itemp == 0 || itemp == 1 || itemp == 2)
    {
        this->flag_hydrostatics = itemp;
    }
    else
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Flag for hydrostatics of the body not available, must be 0, 1 or 2."
           << ".\n";
        throw ValueError(ss.str());
    }

    Logger::debug("    --> Flag hydrostatics: " + std::to_string(flag_hydrostatics));

    // Read body mesh file name
    if (fscanf(pFile, "%s %[^\n]\n", cHydrostaticMeshName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the hydrostatics mesh name of the body: " << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->hydrostaticMeshName = cHydrostaticMeshName;

    Logger::debug("    --> Hydrostatic Mesh Name: " + this->hydrostaticMeshName);

    // Read flag for radiation force
    if (fscanf(pFile, "%d %[^\n]\n", &radiationFlag, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for the radiation force."
           << ".\n";
        throw ValueError(ss.str());
    }

    Logger::debug("    --> Radiation Force Flag: " + std::to_string(radiationFlag));

    // Read flag for first order excitation force
    if (fscanf(pFile, "%d %[^\n]\n", &firstOrderExcitationFlag, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for the first order excitation."
           << ".\n";
        throw ValueError(ss.str());
    }

    Logger::debug("    --> First Order Excitation Flag: " + std::to_string(firstOrderExcitationFlag));

    // Read flag for second order excitation force
    if (fscanf(pFile, "%d %[^\n]\n", &secondOrderExcitationFlag, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for the second order excitation."
           << ".\n";
        throw ValueError(ss.str());
    }

    Logger::debug("    --> Second Order Excitation Flag: " + std::to_string(secondOrderExcitationFlag));

    // Read viscous added mass coefficients
    for (int ii = 0; ii < 6; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the viscous added mass of body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        A_visc(ii, 0) += dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    std::ostringstream oss_av;
    oss_av << "    --> Viscous Added Mass: " << A_visc.t();
    Logger::debug(oss_av.str());

    // Read viscous damping coefficients
    for (int ii = 0; ii < 6; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the viscous damping of body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        B_visc(ii, 0) += dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    std::ostringstream oss_bv;
    oss_bv << "    --> Viscous Damping: " << B_visc.t();
    Logger::debug(oss_bv.str());

    // Read viscous damping coefficients
    for (int ii = 0; ii < 6; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the viscous damping 2 of body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        B_visc2(ii, 0) += dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    std::ostringstream oss_bv2;
    oss_bv2 << "    --> Viscous Damping 2: " << B_visc2.t();
    Logger::debug(oss_bv2.str());

    // Generate array of pointers in order to storage the BCPs pointers
    this->pBodyBcps = new BCP*[this->numBcps];

    if (flag_blocked == 2)
    {
        Logger::info("    --> Reading Body Imposed Movements...");
        ReadLockBodyMovements();
        Logger::info("    --> Body Imposed Movements Read");
    }

    // Generate object of hidrostatic mesh if needed
    if (flag_hydrostatics > 0)
    {
        std::string filename = JoinPath(pSim->inputFolderPath, hydrostaticMeshName);
        // HARCODED: mesh type, the user should be able to choose different meshes
        pNLHSMesh = new BodyTri2DMesh(this->id, filename, this);
        pNLHSMesh->ReadPropertiesASCII();
        pNLHSMesh->Preprocess();
    }
}

void Body::ReadPropertiesYAML(YAML::Node node)
{
    Logger::info("--> Reading Body: " + std::to_string(this->GetId() + 1));

    takeCOGHydroDatabase = node["take_cog_hydro_database"].as<int>();
    Logger::debug(" takeCOGHydroDatabase: " + std::to_string(takeCOGHydroDatabase));

    // Read dofs
    YAML::Node dofsNode = node["dofs"];
    numDofs = (int)dofsNode.size();
    Logger::debug(" numDofs: " + std::to_string(numDofs));
    pDofs = new int[numDofs];
    for (int ii = 0; ii < numDofs; ii++)
    {
        pDofs[ii] = dofsNode[ii].as<int>() - 1;
        isDofActive(pDofs[ii], 0) = 1.0;
    }

    // Read BCP indexes
    YAML::Node bcpsNode = node["bcps_indexes"];
    numBcps = (int)bcpsNode.size();
    pIndexBcps = new int[numBcps];
    for (int ii = 0; ii < numBcps; ii++)
    {
        pIndexBcps[ii] = bcpsNode[ii].as<int>() - 1;
    }

    // Read wind turbine indexes
    YAML::Node wtNode = node["wind_turbines_indexes"];
    numWindTurbs = (int)wtNode.size();
    pBodyWindTurbs = new WindTurbine*[numWindTurbs];
    pIndexWindTurbs = new int[numWindTurbs];
    for (int ii = 0; ii < numWindTurbs; ii++)
    {
        pIndexWindTurbs[ii] = wtNode[ii].as<int>() - 1;
    }

    // Read equilibrium position
    YAML::Node posNode = node["initial_position"];
    for (int ii = 0; ii < 6; ii++)
    {
        double dtemp = posNode[ii].as<double>();
        if (takeCOGHydroDatabase == 0)
        {
            pos_eq(ii, 0) = dtemp;
            pos(ii, 0) = dtemp;
        }
    }

    // Read initial displacement
    YAML::Node dispNode = node["initial_displacement"];
    for (int ii = 0; ii < 6; ii++)
    {
        pos(ii, 0) += dispNode[ii].as<double>();
    }
    pos_ini = pos;

    // Read hydro database
    hydroDatabaseName = node["hydro_database"].as<std::string>();
    hydroDatabaseIndex = node["hydro_database_index"].as<int>() - 1;

    // Read flags
    flag_blocked = node["freedom_flag"].as<int>();
    if (flag_blocked >= 3)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Flag for blocking body not available, must be 0, 1 or 2.\n";
        throw ValueError(ss.str());
    }

    movementsFileName = node["imposed_motion_file"].as<std::string>();
    flag_hydrostatics = node["hydrostatics_flag"].as<int>();
    hydrostaticMeshName = node["hydrostatics_mesh"].as<std::string>();
    radiationFlag = node["radiation_flag"].as<int>();
    firstOrderExcitationFlag = node["excitation_1st_flag"].as<int>();
    secondOrderExcitationFlag = node["excitation_2nd_flag"].as<int>();

    // Read viscous coefficients
    YAML::Node avNode = node["viscous_added_mass"];
    for (int ii = 0; ii < 6; ii++)
        A_visc(ii, 0) += avNode[ii].as<double>();

    YAML::Node bvNode = node["viscous_linear_damping"];
    for (int ii = 0; ii < 6; ii++)
        B_visc(ii, 0) += bvNode[ii].as<double>();

    YAML::Node bv2Node = node["viscous_quadratic_damping"];
    for (int ii = 0; ii < 6; ii++)
        B_visc2(ii, 0) += bv2Node[ii].as<double>();

    // Generate array of pointers for BCPs
    pBodyBcps = new BCP*[numBcps];

    if (flag_blocked == 2)
    {
        Logger::info("    --> Reading Body Imposed Movements...");
        ReadLockBodyMovements();
        Logger::info("    --> Body Imposed Movements Read");
    }

    // Generate hydrostatic mesh if needed
    if (flag_hydrostatics > 0)
    {
        std::string filename = JoinPath(pSim->inputFolderPath, hydrostaticMeshName);
        pNLHSMesh = new BodyTri2DMesh(this->id, filename, this);
        pNLHSMesh->ReadPropertiesASCII();
        pNLHSMesh->Preprocess();
    }
}

void Body::ReadLockBodyMovements(void)
{

    // Declare variables
    char buffer_line[1000];
    char cMovementsTimeSeriesFileName[1000];
    double dtemp;
    int itemp;

    std::string file_path = JoinPath(pSim->inputFolderPath, movementsFileName);
    Logger::info("             Reading from file: " + file_path);

    // Open file
    FILE* pFile = fopen(file_path.c_str(), "r");
    // Discard header lines and check for body movements type
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), pFile);
    }
    if (fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read type of movement"
           << ".\n";
        throw ValueError(ss.str());
    }
    movementTypeFlag = itemp;
    if (movementTypeFlag == 1)
    {
        Logger::info("             Reading Harmonic Analytic movements...");
        // Discard header lines and loop over dofs
        for (int ii = 0; ii < 3; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }
        for (int ii = 0; ii < 6; ii++)
        {
            Logger::debug("             Reading DOF " + std::to_string(ii));
            // Ignore dof name line
            fgets(buffer_line, sizeof(buffer_line), pFile);
            // Read offset
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement offset of body: " << this->GetId()
                   << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            offset(ii, 0) = dtemp;
            // Read amplitude
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement amplitude of body: " << this->GetId()
                   << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            amplitude(ii, 0) = dtemp;
            // Read period
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement period of body: " << this->GetId()
                   << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            period(ii, 0) = dtemp;
            // Read phase
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement phase of body: " << this->GetId()
                   << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            phase(ii, 0) = dtemp;
        }
        omega = 2.0 * arma::datum::pi / period;
        phase = phase * arma::datum::pi / 180.0;
    }
    else if (movementTypeFlag == 2)
    {
        // Discard lines from analytic solution
        for (int ii = 0; ii < 36; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }
        Logger::info("             Reading Time series data movements...");
        // Read movements filename
        if (fscanf(pFile, "%s %[^\n]\n", cMovementsTimeSeriesFileName, buffer_line) != 2)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the movements time series file name of the body: "
               << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        movementsTimeSeriesFileName = cMovementsTimeSeriesFileName;
    }
    else
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Type of movement can only be 1 or 2"
           << ".\n";
        throw ValueError(ss.str());
    }
    // Close file
    fclose(pFile);

    if (movementTypeFlag == 2)
    {

        // Inicializo la variable donde guardar el numero de pasos temporales
        int nt;
        // Abro el fichero
        file_path = JoinPath(pSim->inputFolderPath, movementsTimeSeriesFileName);
        Logger::info("             Reading from file: " + file_path);
        std::ifstream datosPosF(file_path);
        // Leo el numero de pasos temporales a leer
        datosPosF >> nt;
        datosPosF.ignore(std::numeric_limits<int>::max(), '\n');

        if (nt < 2)
        {
            std::stringstream ss;
            ss << "Body: " << this->GetId() << " - Movement time series does not have enough data. nt = " << nt
               << ".\n";
            throw ValueError(ss.str());
        }

        // Alocato la matriz que contiene la informacion
        timeFixed = arma::zeros(nt, 1);
        posFixed = arma::zeros(nt, 6);
        velFixed = arma::zeros(nt, 6);
        accFixed = arma::zeros(nt, 6);
        // Leo toda la info
        for (int i = 0; i < nt; i = i + 1)
        {
            datosPosF >> timeFixed(i, 0) >> posFixed(i, 0) >> posFixed(i, 1) >> posFixed(i, 2) >> posFixed(i, 3) >>
                posFixed(i, 4) >> posFixed(i, 5) >> velFixed(i, 0) >> velFixed(i, 1) >> velFixed(i, 2) >>
                velFixed(i, 3) >> velFixed(i, 4) >> velFixed(i, 5) >> accFixed(i, 0) >> accFixed(i, 1) >>
                accFixed(i, 2) >> accFixed(i, 3) >> accFixed(i, 4) >> accFixed(i, 5);
        }
        // Cierro el fichero
        datosPosF.close();

        if (timeFixed(0, 0) > 0.0)
        {
            std::stringstream ss;
            ss << "Body: " << this->GetId() << " - Movement time series does not start in zero"
               << ".\n";
            throw ValueError(ss.str());
        }

        if (timeFixed(nt - 1, 0) < pSim->simulationTime)
        {
            std::stringstream ss;
            ss << "Body: " << this->GetId() << " - Movement time series is not long enough for simulation time"
               << ".\n"
               << "pSim->simulationTime = " << pSim->simulationTime
               << " s;  timeFixed(nt-1,0) = " << timeFixed(nt - 1, 0) << " s; nt = " << nt << ".\n";
            throw ValueError(ss.str());
        }
    }

    UpdateLockBody(0.0);
    pos_ini = pos;
}

//
void Body::StoreVelocities(bool restoreMatrix)
{
    if (!restoreMatrix)
    {
        velBuffer.submat(0, pSim->timeBufferCount, 5, pSim->timeBufferCount) = vel;
    }
    else
    {
        int num_points_irf = this->pHydro->GetNumPointsIrf();
        arma::mat velBufferNew = arma::zeros(6, pSim->timeBufferSize);
        velBufferNew.cols(0, num_points_irf - 1) =
            velBuffer.cols(pSim->timeBufferSize - num_points_irf, pSim->timeBufferSize - 1);
        velBuffer = velBufferNew;
        velBufferCount = num_points_irf - 1;
    }
}

// Actualiza valores del BCP
void Body::UpdateBcps(void)
{
    // Datos necesarios para las matrices de rotacion
    double roll = arma::as_scalar(pos(3, 0));
    double pitch = arma::as_scalar(pos(4, 0));
    double yaw = arma::as_scalar(pos(5, 0));
    double roll_dot = arma::as_scalar(vel(3, 0));
    double pitch_dot = arma::as_scalar(vel(4, 0));
    double yaw_dot = arma::as_scalar(vel(5, 0));
    double roll_dot2 = arma::as_scalar(acc(3, 0));
    double pitch_dot2 = arma::as_scalar(acc(4, 0));
    double yaw_dot2 = arma::as_scalar(acc(5, 0));
    double cr = cos(roll);
    double sr = sin(roll);
    double cp = cos(pitch);
    double sp = sin(pitch);
    double cy = cos(yaw);
    double sy = sin(yaw);

    // Matriz de rotacion calculada como Mz*My*Mx
    arma::mat R1 = arma::eye(3, 3);
    R1(1, 1) = cr;
    R1(1, 2) = -sr;
    R1(2, 1) = sr;
    R1(2, 2) = cr;
    arma::mat R2 = arma::eye(3, 3);
    R2(2, 2) = cp;
    R2(2, 0) = -sp;
    R2(0, 2) = sp;
    R2(0, 0) = cp;
    arma::mat R3 = arma::eye(3, 3);
    R3(0, 0) = cy;
    R3(0, 1) = -sy;
    R3(1, 0) = sy;
    R3(1, 1) = cy;

    rotMat = R3 * R2 * R1;
    invRotMat = rotMat.t();

    arma::mat R1_dot = arma::zeros(3, 3);
    R1_dot(1, 1) = -sr;
    R1_dot(1, 2) = -cr;
    R1_dot(2, 1) = cr;
    R1_dot(2, 2) = -sr;
    R1_dot = roll_dot * R1_dot;
    arma::mat R2_dot = arma::zeros(3, 3);
    R2_dot(2, 2) = -sp;
    R2_dot(2, 0) = -cp;
    R2_dot(0, 2) = cp;
    R2_dot(0, 0) = -sp;
    R2_dot = pitch_dot * R2_dot;
    arma::mat R3_dot = arma::zeros(3, 3);
    R3_dot(0, 0) = -sy;
    R3_dot(0, 1) = -cy;
    R3_dot(1, 0) = cy;
    R3_dot(1, 1) = -sy;
    R3_dot = yaw_dot * R3_dot;

    rotMat_dot = R3_dot * R2 * R1 + R3 * R2_dot * R1 + R3 * R2 * R1_dot;

    arma::mat R1_dot2 = -R1;
    R1_dot2(0, 0) = 0;
    R1_dot2 = roll_dot2 * R1_dot + roll_dot * roll_dot * R1_dot2;
    arma::mat R2_dot2 = -R2;
    R2_dot2(1, 1) = 0;
    R2_dot2 = pitch_dot2 * R2_dot + pitch_dot * pitch_dot * R2_dot2;
    arma::mat R3_dot2 = -R3;
    R3_dot2(2, 2) = 0;
    R3_dot2 = yaw_dot2 * R3_dot + yaw_dot * yaw_dot * R3_dot2;

    rotMat_dot2 = R3_dot2 * R2 * R1 + R3_dot * R2_dot * R1 + R3_dot * R2 * R1_dot + R3_dot * R2_dot * R1 +
                  R3 * R2_dot2 * R1 + R3 * R2_dot * R1_dot + R3_dot * R2 * R1_dot + R3 * R2_dot * R1_dot +
                  R3 * R2 * R1_dot2;

    aMat = R3;
    aMat(0, 0) = cp * cy;
    aMat(1, 0) = cp * sy;
    aMat(2, 0) = -sp;

    aMat_dot = R3_dot;
    aMat_dot(0, 0) = -pitch_dot * sp * cy - yaw_dot * cp * sy;
    aMat_dot(1, 0) = -pitch_dot * sp * sy + yaw_dot * cp * cy;
    aMat_dot(2, 0) = -pitch_dot * cp;

    // Variable temporal
    arma::mat posG_temp, posL_temp;
    for (int ii = 0; ii < numBcps; ii++)
    {
        // Bucle sobre todos los BCPs
        posL_temp = pBodyBcps[ii]->posWrtCdgLocal;
        posG_temp = rotMat * posL_temp;             // Brazo COG-bcp en global
        pBodyBcps[ii]->posWrtCdgGlobal = posG_temp; // Brazo COG-bcp en global
        pBodyBcps[ii]->rotMat = rotMat;             // Matriz de rotacion
        pBodyBcps[ii]->rotMat_dot = rotMat_dot;     // Matriz de rotacion
        // Posicion del BCP en global, lo mismo para vel y acc.
        pBodyBcps[ii]->posG_BCP.rows(0, 2) = pos.rows(0, 2) + posG_temp; // Posicion del BCP en global
        pBodyBcps[ii]->velG_BCP.rows(0, 2) = vel.rows(0, 2) + rotMat_dot * posL_temp;
        pBodyBcps[ii]->accG_BCP.rows(0, 2) = acc.rows(0, 2) + rotMat_dot2 * posL_temp;
    }
}

// Resetea fuerzas BCPs
void Body::ResetBcps(void)
{
    for (int ii = 0; ii < numBcps; ii++)
    {
        // Reseteo a cero la fuerza sobre el BCP
        pBodyBcps[ii]->forceBcp = arma::zeros(6, 1);
        pBodyBcps[ii]->temp = arma::zeros(6, 1);
    }
    // Reseteo a cero la fuerza total de todos los BCPs
    bcpForces = arma::zeros(6, 1);
}

void Body::UpdateLockBody(double time)
{
    if (flag_blocked == 1)
    {
        pos = pos_ini;
        vel = arma::zeros(6, 1);
        acc = arma::zeros(6, 1);
    }
    else if (flag_blocked == 2)
    {

        if (movementTypeFlag == 1)
        {
            for (int ii = 0; ii < 6; ii++)
            {
                pos(ii, 0) = offset(ii, 0) + amplitude(ii, 0) * sin(time * omega(ii, 0) + phase(ii, 0));
                vel(ii, 0) = omega(ii, 0) * amplitude(ii, 0) * cos(time * omega(ii, 0) + phase(ii, 0));
                vel(ii, 0) = -omega(ii, 0) * omega(ii, 0) * amplitude(ii, 0) * sin(time * omega(ii, 0) + phase(ii, 0));
            }
        }
        else if (movementTypeFlag == 2)
        {
            arma::mat tmp;
            tmp = time * arma::ones(1, 1);
            pos = (interp1(timeFixed, posFixed, tmp)).t();
            vel = (interp1(timeFixed, velFixed, tmp)).t();
            acc = (interp1(timeFixed, accFixed, tmp)).t();
        }
    }
}

//
void Body::OpenOutputFilesASCII(std::string path)
{
    char buffer1[50], buffer2[50], buffer3[50], buffer4[50], buffer5[50], buffer6[50];
    char buffer7[50], buffer8[50], buffer9[50], buffer10[50], buffer11[50];

    int nn1 = sprintf(buffer1, "DOF_1_Body_%d.txt", GetId());
    int nn2 = sprintf(buffer2, "DOF_2_Body_%d.txt", GetId());
    int nn3 = sprintf(buffer3, "DOF_3_Body_%d.txt", GetId());
    int nn4 = sprintf(buffer4, "DOF_4_Body_%d.txt", GetId());
    int nn5 = sprintf(buffer5, "DOF_5_Body_%d.txt", GetId());
    int nn6 = sprintf(buffer6, "DOF_6_Body_%d.txt", GetId());
    int nn7 = sprintf(buffer7, "HydroStiffnessForce_Body_%d.txt", GetId());
    int nn8 = sprintf(buffer8, "WaveRadiationForce_Body_%d.txt", GetId());
    int nn9 = sprintf(buffer9, "BCPForce_Body_%d.txt", GetId());
    int nn10 = sprintf(buffer10, "WaveExcitationForce_Body_%d.txt", GetId());
    int nn11 = sprintf(buffer11, "WindTurbineForce_Body_%d.txt", GetId());

    std::string file_path1 = JoinPath(path, buffer1);
    std::string file_path2 = JoinPath(path, buffer2);
    std::string file_path3 = JoinPath(path, buffer3);
    std::string file_path4 = JoinPath(path, buffer4);
    std::string file_path5 = JoinPath(path, buffer5);
    std::string file_path6 = JoinPath(path, buffer6);
    std::string file_path7 = JoinPath(path, buffer7);
    std::string file_path8 = JoinPath(path, buffer8);
    std::string file_path9 = JoinPath(path, buffer9);
    std::string file_path10 = JoinPath(path, buffer10);
    std::string file_path11 = JoinPath(path, buffer11);

    pfile_DOF_1 = fopen(file_path1.c_str(), "w");
    if (pfile_DOF_1 == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn1 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_DOF_2 = fopen(file_path2.c_str(), "w");
    if (pfile_DOF_2 == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn2 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_DOF_3 = fopen(file_path3.c_str(), "w");
    if (pfile_DOF_3 == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn3 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_DOF_4 = fopen(file_path4.c_str(), "w");
    if (pfile_DOF_4 == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn4 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_DOF_5 = fopen(file_path5.c_str(), "w");
    if (pfile_DOF_5 == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn5 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_DOF_6 = fopen(file_path6.c_str(), "w");
    if (pfile_DOF_6 == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn6 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_HSF = fopen(file_path7.c_str(), "w");
    if (pfile_HSF == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn7 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_WRF = fopen(file_path8.c_str(), "w");
    if (pfile_WRF == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn8 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_BCPF = fopen(file_path9.c_str(), "w");
    if (pfile_BCPF == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn9 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_WEF = fopen(file_path10.c_str(), "w");
    if (pfile_WEF == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn10 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    pfile_WindTurb = fopen(file_path11.c_str(), "w");
    if (pfile_WindTurb == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << nn11 << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
}

//
void Body::CloseOutputFilesASCII(void)
{
    fclose(pfile_DOF_1);
    fclose(pfile_DOF_2);
    fclose(pfile_DOF_3);
    fclose(pfile_DOF_4);
    fclose(pfile_DOF_5);
    fclose(pfile_DOF_6);
    fclose(pfile_HSF);
    fclose(pfile_WRF);
    fclose(pfile_BCPF);
    fclose(pfile_WEF);
    fclose(pfile_WindTurb);
}

//
void Body::OpenOutputFilesCSV(std::string path)
{
    int bid = GetId();
    char buffer[100];

    // --- Motion CSV ---
    sprintf(buffer, "body_%d_motion.csv", bid);
    std::string motion_path = JoinPath(path, buffer);
    pfile_motion_csv = fopen(motion_path.c_str(), "w");
    if (pfile_motion_csv == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << buffer << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    // Write header: time, positions (x,y,z,rl,pt,yw), velocities, accelerations
    fprintf(pfile_motion_csv, "time");
    const char* dof_names[] = {"x", "y", "z", "rl", "pt", "yw"};
    for (int i = 0; i < 6; i++)
        fprintf(pfile_motion_csv, ",b%d_%s", bid, dof_names[i]);
    for (int i = 0; i < 6; i++)
        fprintf(pfile_motion_csv, ",b%d_v%s", bid, dof_names[i]);
    for (int i = 0; i < 6; i++)
        fprintf(pfile_motion_csv, ",b%d_a%s", bid, dof_names[i]);
    fprintf(pfile_motion_csv, "\n");

    // --- Forces CSV ---
    sprintf(buffer, "body_%d_forces.csv", bid);
    std::string forces_path = JoinPath(path, buffer);
    pfile_forces_csv = fopen(forces_path.c_str(), "w");
    if (pfile_forces_csv == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: " << buffer << "\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    // Write header
    const char* force_names[] = {"fx", "fy", "fz", "mx", "my", "mz"};
    fprintf(pfile_forces_csv, "time");
    // Hydrostatic forces
    for (int i = 0; i < 6; i++)
        fprintf(pfile_forces_csv, ",b%d_hs_%s", bid, force_names[i]);
    // Radiation forces
    for (int i = 0; i < 6; i++)
        fprintf(pfile_forces_csv, ",b%d_rad_%s", bid, force_names[i]);
    // Excitation forces (1st order)
    for (int i = 0; i < 6; i++)
        fprintf(pfile_forces_csv, ",b%d_exc1_%s", bid, force_names[i]);
    // Excitation forces (2nd order)
    for (int i = 0; i < 6; i++)
        fprintf(pfile_forces_csv, ",b%d_exc2_%s", bid, force_names[i]);
    // BCP forces (per BCP + total)
    for (int b = 0; b < numBcps; b++)
        for (int i = 0; i < 6; i++)
            fprintf(pfile_forces_csv, ",b%d_bcp%d_%s", bid, b, force_names[i]);
    for (int i = 0; i < 6; i++)
        fprintf(pfile_forces_csv, ",b%d_bcp_total_%s", bid, force_names[i]);
    // Wind turbine forces (per turbine + total)
    for (int w = 0; w < numWindTurbs; w++)
        for (int i = 0; i < 6; i++)
            fprintf(pfile_forces_csv, ",b%d_wt%d_%s", bid, w, force_names[i]);
    for (int i = 0; i < 6; i++)
        fprintf(pfile_forces_csv, ",b%d_wt_total_%s", bid, force_names[i]);
    fprintf(pfile_forces_csv, "\n");
}

void Body::CloseOutputFilesCSV(void)
{
    if (pfile_motion_csv)
        fclose(pfile_motion_csv);
    if (pfile_forces_csv)
        fclose(pfile_forces_csv);
}

void Body::OpenOutputFiles(std::string path)
{
    if (pSim->outputFormat == 1)
        OpenOutputFilesCSV(path);
    else
        OpenOutputFilesASCII(path);
}

void Body::CloseOutputFiles(void)
{
    if (pSim->outputFormat == 1)
        CloseOutputFilesCSV();
    else
        CloseOutputFilesASCII();
}

// Escibir datos a fichero
void Body::WriteOut(double t)
{
    if (pSim->outputFormat == 1)
    {
        // CSV format: motion file
        fprintf(pfile_motion_csv, "%f", t);
        for (int i = 0; i < 6; i++)
            fprintf(pfile_motion_csv, ",%f", this->pos(i, 0));
        for (int i = 0; i < 6; i++)
            fprintf(pfile_motion_csv, ",%f", this->vel(i, 0));
        for (int i = 0; i < 6; i++)
            fprintf(pfile_motion_csv, ",%f", this->acc(i, 0));
        fprintf(pfile_motion_csv, "\n");

        // CSV format: forces file
        fprintf(pfile_forces_csv, "%f", t);
        for (int i = 0; i < 6; i++)
            fprintf(pfile_forces_csv, ",%f", hydrostaticForces(i, 0));
        for (int i = 0; i < 6; i++)
            fprintf(pfile_forces_csv, ",%f", radiationForces(i, 0));
        for (int i = 0; i < 6; i++)
            fprintf(pfile_forces_csv, ",%f", excitationForces_1(i, 0));
        for (int i = 0; i < 6; i++)
            fprintf(pfile_forces_csv, ",%f", excitationForces_2(i, 0));
        for (int ii = 0; ii < numBcps; ii++)
            for (int jj = 0; jj < 6; jj++)
                fprintf(pfile_forces_csv, ",%f", pBodyBcps[ii]->temp(jj, 0));
        for (int jj = 0; jj < 6; jj++)
            fprintf(pfile_forces_csv, ",%f", bcpForces(jj, 0));
        for (int ii = 0; ii < numWindTurbs; ii++)
            for (int jj = 0; jj < 6; jj++)
                fprintf(pfile_forces_csv, ",%f", pBodyWindTurbs[ii]->forceBodyCOG(jj, 0));
        for (int jj = 0; jj < 6; jj++)
            fprintf(pfile_forces_csv, ",%f", windTurbForces(jj, 0));
        fprintf(pfile_forces_csv, "\n");
    }
    else
    {
        // ASCII format (original)
        fprintf(pfile_DOF_1, "%f    %f    %f    %f \n", t, this->pos(0, 0), this->vel(0, 0), this->acc(0, 0));
        fprintf(pfile_DOF_2, "%f    %f    %f    %f \n", t, this->pos(1, 0), this->vel(1, 0), this->acc(1, 0));
        fprintf(pfile_DOF_3, "%f    %f    %f    %f \n", t, this->pos(2, 0), this->vel(2, 0), this->acc(2, 0));
        fprintf(pfile_DOF_4, "%f    %f    %f    %f \n", t, this->pos(3, 0), this->vel(3, 0), this->acc(3, 0));
        fprintf(pfile_DOF_5, "%f    %f    %f    %f \n", t, this->pos(4, 0), this->vel(4, 0), this->acc(4, 0));
        fprintf(pfile_DOF_6, "%f    %f    %f    %f \n", t, this->pos(5, 0), this->vel(5, 0), this->acc(5, 0));
        fprintf(pfile_HSF, "%f    ", t);
        for (int ii = 0; ii < 6; ii = ii + 1)
            fprintf(pfile_HSF, "%f    ", hydrostaticForces(ii, 0));
        fprintf(pfile_HSF, "\n");

        fprintf(pfile_WRF, "%f    ", t);
        for (int ii = 0; ii < 6; ii = ii + 1)
            fprintf(pfile_WRF, "%f    ", radiationForces(ii, 0));
        fprintf(pfile_WRF, "\n");

        fprintf(pfile_WEF, "%f    ", t);
        for (int ii = 0; ii < 6; ii = ii + 1)
            fprintf(pfile_WEF, "%f    ", excitationForces_1(ii, 0));
        for (int ii = 0; ii < 6; ii = ii + 1)
            fprintf(pfile_WEF, "%f    ", excitationForces_2(ii, 0));
        fprintf(pfile_WEF, "\n");

        fprintf(pfile_BCPF, "%f    ", t);
        for (int ii = 0; ii < numBcps; ii = ii + 1)
        {
            for (int jj = 0; jj < 6; jj = jj + 1)
            {
                fprintf(pfile_BCPF, "%f    ", pBodyBcps[ii]->temp(jj, 0));
            }
        }
        for (int jj = 0; jj < 6; jj = jj + 1)
        {
            fprintf(pfile_BCPF, "%f    ", bcpForces(jj, 0));
        }
        fprintf(pfile_BCPF, "\n");

        fprintf(pfile_WindTurb, "%f    ", t);
        for (int ii = 0; ii < numWindTurbs; ii = ii + 1)
        {
            for (int jj = 0; jj < 6; jj = jj + 1)
            {
                fprintf(pfile_WindTurb, "%f    ", pBodyWindTurbs[ii]->forceBodyCOG(jj, 0));
            }
        }
        for (int jj = 0; jj < 6; jj = jj + 1)
        {
            fprintf(pfile_WindTurb, "%f    ", windTurbForces(jj, 0));
        }
        fprintf(pfile_WindTurb, "\n");
    }
}
