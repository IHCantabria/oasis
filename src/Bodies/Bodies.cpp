
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

//
Body::Body(int n, Simulation *pIncSim)
{
    id = n;
    pSim = pIncSim;
}

// Calcula el efecto de las fuerzas sobre los BCPs sobre su CDG
void Body::ComputeBcpForces(void)
{

    // Variables locales necesarias usadas a continuacion
    arma::mat posG_temp;
    arma::mat ForceBCP_temp;
    arma::mat F_M;
    arma::mat M_F;
    double M_norm, det;
    double rx, ry, rz, Mx, My, Mz;

    for (int ii = 0; ii < numBcps; ii++)
    {

        posG_temp = pBodyBcps[ii]->posWrtCdgGlobal; // Guardo en variable temporal la posicion global del BCP
        ForceBCP_temp = pBodyBcps[ii]->forceBcp;    // Guardo en variable temporal las fuerzas y momentos sobre el BCP

        M_F = arma::cross(posG_temp, ForceBCP_temp.rows(0, 2)); // Momento sobre el cdg causado por la fuerza en el bcp, en global

        bcpForces.rows(0, 2) = bcpForces.rows(0, 2) + ForceBCP_temp.rows(0, 2);                      // Acumulo la fuerza total sobre el cdg en global.
        bcpForces.rows(3, 5) = bcpForces.rows(3, 5) + rotMat.t() * (ForceBCP_temp.rows(3, 5) + M_F); // Acumulo el momento total sobre el cdg en local.
    }

    if (bcpForces.has_nan())
    {
        std::cout << std::endl
                  << "ERROR: NaN Detected on body with id = " << id << std::endl;
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

void Body::LoadHydrodynamicDatabase(Body **hydroDatabaseBodies)
{
    std::cout << "--> Reading Hydrodynamics Properties (HDF5 format)" << std::endl;

    // File path
    std::string file_path = JoinPath(this->pSim->inputFolderPath, this->hydroDatabaseName);

    // Load hydrodynamic database
    this->pHydro = new HydroDatabase(this->hydroDatabaseIndex, this->id, hydroDatabaseBodies, this->pSim);
    std::cout << "    --> Loading hydrodynamic database" << std::endl;
    this->pHydro->LoadHydrodynamicData(file_path);

    std::cout << "    --> Loading and checking the Hydrodynamic C.O.G position" << std::endl;
    // Load and check the Hydrodynamic C.O.G position
    if (this->takeCOGHydroDatabase == 1)
    {
        // Load COG equilibrium position
        this->pos_eq.rows(0, 2) = this->pHydro->GetCog().t();

        std::cout << "      --> Using HDB COG: " << this->pHydro->GetCog();

        // Add initial displacement to the COG equilibrium position
        this->pos = this->pos + this->pos_eq;
        pos_ini = pos;

        std::cout << "      --> Initial Position: " << this->pos.t();
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

    std::cout << "----> Hydrodynamic Properties Read" << std::endl;
}

// Read body properties from the input file
void Body::ReadPropertiesASCII(FILE *pFile)
{
    // Declare variables
    char buffer_line[1000];
    fpos_t carriage_init;
    char cHydroDatabaseName[1000];
    char cMovementsFileName[1000];
    char cHydrostaticMeshName[1000];
    double dtemp;
    int itemp;

    std::cout << "--> Reading Body: " << this->GetId() + 1 << std::endl;

    // Read flag to take COG from Hydrodynamic database
    if (fscanf(pFile, "%d %[^\n]\n", &takeCOGHydroDatabase, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag to take COG from Hydrodynamic database."
           << ".\n";
        throw ValueError(ss.str());
    }

    std::cout << " takeCOGHydroDatabase: " << takeCOGHydroDatabase << std::endl;

    // Read dofs considered
    fgetpos(pFile, &carriage_init);
    while (fscanf(pFile, "%d", &itemp) == 1)
    {
        this->numDofs++;
    }
    fsetpos(pFile, &carriage_init);

    std::cout << " numDofs: " << numDofs << std::endl;

    this->pDofs = new int[this->numDofs];
    for (int ii = 0; ii < this->numDofs; ii++)
    {
        if (fscanf(pFile, "%d", &itemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the deegres of freedom for body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        this->pDofs[ii] = itemp - 1;
        isDofActive(itemp - 1, 0) = 1.0;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    std::cout << "    --> Degrees of Freedom: " << this->pDofs << std::endl;

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
            ss << "An error occurred when trying to read the boundary condition points of the body: " << this->GetId() << " "
               << "\n";
            throw ValueError(ss.str());
        }
        this->pIndexBcps[ii] = itemp - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    std::cout << "    --> Boundary Condition Points: " << this->pIndexBcps << std::endl;

    // Read the wind turbines in the body
    fgetpos(pFile, &carriage_init);
    while (fscanf(pFile, "%d", &itemp) == 1)
    {
        this->numWindTurbs++;
    }
    fsetpos(pFile, &carriage_init);

    this->pBodyWindTurbs = new WindTurbine *[this->numWindTurbs];
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

    std::cout << "    --> Wind Turbines: " << this->pIndexWindTurbs << std::endl;

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
        std::cout << "    --> Equilibrium Position: " << this->pos.t() << std::endl;
    }

    // Read Initial displacement from reference position
    for (int ii = 0; ii < 6; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the initial displacement of the body: " << this->GetId() << "\n";
            throw ValueError(ss.str());
        }
        pos(ii, 0) += dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    pos_ini = pos;

    std::cout << "    --> Initial Disp: " << pos.t() << std::endl;

    // Read hydrodynamic database filename
    if (fscanf(pFile, "%s %[^\n]\n", cHydroDatabaseName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the hydrodynamic database name of the body: " << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->hydroDatabaseName = cHydroDatabaseName;

    std::cout << "    --> hydrodynamic Database: " << this->hydroDatabaseName << std::endl;

    // Read body index in the associated database
    if (fscanf(pFile, "%d %[^\n]\n", &(this->hydroDatabaseIndex), buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the index of the body in the hydrodynamic database: " << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->hydroDatabaseIndex--;

    std::cout << "    --> hydrodynamic Database Index: " << this->hydroDatabaseIndex + 1 << std::endl;

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

    std::cout << "    --> Flag Blocked: " << flag_blocked << std::endl;

    // Read movements filename
    if (fscanf(pFile, "%s %[^\n]\n", cMovementsFileName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the movements file name of the body: " << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->movementsFileName = cMovementsFileName;

    std::cout << "    --> Movements File Name: " << this->movementsFileName << std::endl;

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

    std::cout << "    --> Flag hydrostatics: " << flag_hydrostatics << std::endl;

    // Read body mesh file name
    if (fscanf(pFile, "%s %[^\n]\n", cHydrostaticMeshName, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the hydrostatics mesh name of the body: " << this->GetId() << "\n";
        throw ValueError(ss.str());
    }
    this->hydrostaticMeshName = cHydrostaticMeshName;

    std::cout << "    --> Hydrostatic Mesh Name: " << this->hydrostaticMeshName << std::endl;

    // Read flag for radiation force
    if (fscanf(pFile, "%d %[^\n]\n", &radiationFlag, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for the radiation force."
           << ".\n";
        throw ValueError(ss.str());
    }

    std::cout << "    --> Radiation Force Flag: " << firstOrderExcitationFlag << std::endl;

    // Read flag for first order excitation force
    if (fscanf(pFile, "%d %[^\n]\n", &firstOrderExcitationFlag, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for the first order excitation."
           << ".\n";
        throw ValueError(ss.str());
    }

    std::cout << "    --> First Order Excitation Flag: " << firstOrderExcitationFlag << std::endl;

    // Read flag for second order excitation force
    if (fscanf(pFile, "%d %[^\n]\n", &secondOrderExcitationFlag, buffer_line) != 2)
    {
        std::stringstream ss;
        ss << "Body: " << this->GetId() << " - Not possible to read flag for the second order excitation."
           << ".\n";
        throw ValueError(ss.str());
    }

    std::cout << "    --> Second Order Excitation Flag: " << secondOrderExcitationFlag << std::endl;

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

    std::cout << "    --> Viscous Added Mass: " << A_visc.t() << std::endl;

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

    std::cout << "    --> Viscous Damping: " << B_visc.t() << std::endl;

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

    std::cout << "    --> Viscous Damping 2: " << B_visc2.t() << std::endl;

    // Generate array of pointers in order to storage the BCPs pointers
    this->pBodyBcps = new BCP *[this->numBcps];

    if (flag_blocked == 2)
    {
        std::cout << "    --> Reading Body Imposed Movements..." << std::endl;
        ReadLockBodyMovements();
        std::cout << "    --> Body Imposed Movements Read" << std::endl;
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

void Body::ReadLockBodyMovements(void)
{

    // Declare variables
    char buffer_line[1000];
    char cMovementsTimeSeriesFileName[1000];
    double dtemp;
    int itemp;

    std::string file_path = JoinPath(pSim->inputFolderPath, movementsFileName);
    std::cout << "             Reading from file: " << file_path << std::endl;

    // Open file
    FILE *pFile = fopen(file_path.c_str(), "r");
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
        std::cout << "             Reading Harmonic Analytic movements..." << std::endl;
        // Discard header lines and loop over dofs
        for (int ii = 0; ii < 3; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }
        for (int ii = 0; ii < 6; ii++)
        {
            std::cout << "             Reading DOF " << ii << std::endl;
            // Ignore dof name line
            fgets(buffer_line, sizeof(buffer_line), pFile);
            // Read offset
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement offset of body: " << this->GetId() << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            offset(ii, 0) = dtemp;
            // Read amplitude
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement amplitude of body: " << this->GetId() << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            amplitude(ii, 0) = dtemp;
            // Read period
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement period of body: " << this->GetId() << " in DOF " << ii << "\n";
                throw ValueError(ss.str());
            }
            period(ii, 0) = dtemp;
            // Read phase
            if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read analytic movement phase of body: " << this->GetId() << " in DOF " << ii << "\n";
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
        std::cout << "             Reading Time series data movements..." << std::endl;
        // Read movements filename
        if (fscanf(pFile, "%s %[^\n]\n", cMovementsTimeSeriesFileName, buffer_line) != 2)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the movements time series file name of the body: " << this->GetId() << "\n";
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
        std::cout << "             Reading from file: " << file_path << std::endl;
        std::ifstream datosPosF(file_path);
        // Leo el numero de pasos temporales a leer
        datosPosF >> nt;
        datosPosF.ignore(std::numeric_limits<int>::max(), '\n');

        if (nt < 2)
        {
            std::stringstream ss;
            ss << "Body: " << this->GetId() << " - Movement time series does not have enough data. nt = " << nt << ".\n";
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
            datosPosF >> timeFixed(i, 0) >> posFixed(i, 0) >> posFixed(i, 1) >> posFixed(i, 2) >> posFixed(i, 3) >> posFixed(i, 4) >> posFixed(i, 5) >> velFixed(i, 0) >> velFixed(i, 1) >> velFixed(i, 2) >> velFixed(i, 3) >> velFixed(i, 4) >> velFixed(i, 5) >> accFixed(i, 0) >> accFixed(i, 1) >> accFixed(i, 2) >> accFixed(i, 3) >> accFixed(i, 4) >> accFixed(i, 5);
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
               << "pSim->simulationTime = " << pSim->simulationTime << " s;  timeFixed(nt-1,0) = " << timeFixed(nt - 1, 0) << " s; nt = " << nt << ".\n";
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
        velBufferNew.cols(0, num_points_irf - 1) = velBuffer.cols(pSim->timeBufferSize - num_points_irf, pSim->timeBufferSize - 1);
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

    rotMat_dot2 = R3_dot2 * R2 * R1 + R3_dot * R2_dot * R1 + R3_dot * R2 * R1_dot +
                  R3_dot * R2_dot * R1 + R3 * R2_dot2 * R1 + R3 * R2_dot * R1_dot +
                  R3_dot * R2 * R1_dot + R3 * R2_dot * R1_dot + R3 * R2 * R1_dot2;

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
        posG_temp = rotMat * posL_temp;             // Brazo cdg-bcp en global
        pBodyBcps[ii]->posWrtCdgGlobal = posG_temp; // Brazo cdg-bcp en global
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

// Escibir datos a fichero
void Body::WriteOut(double t)
{
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
            // fprintf(pfile_BCPF, "%f    ", pBodyBcps[ii]->forceBcp(jj, 0));
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
