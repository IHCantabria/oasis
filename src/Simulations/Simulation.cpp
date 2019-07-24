
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <ctime>
#include "Simulation.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"


arma::mat Simulation::CalculateSystemDynamics(double time, arma::mat y)
{
    numCallsSysFun++;

    //std::cout << "Main::fun - At first" << std::endl; 
	arma::mat yprime = arma::zeros(size(y));
	int i0;
	//std::cout << "Main::fun - Before copy y to objects" << std::endl;
	// Copy info from y to the objects.
	int ini = 0;
	for(int ii=0; ii<numBodies; ii++)
	{
		pBodies[ii]->pos = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0; ii<numLines;ii++)
	{
		for(int jj=0; jj<pLines[ii]->N; jj++)
        {
			pLines[ii]->pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	ini = ini + numWinches;
	for(int ii=0; ii<numBodies; ii++)
	{
		pBodies[ii]->vel = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0; ii<numLines; ii++)
	{
		for(int jj=0; jj<pLines[ii]->N; jj++)
		{
			pLines[ii]->vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0; ii<numWinches; ii++)
    {
		pWinches[ii]->theta = arma::as_scalar(y.row(numSystem2-(ii+1)));
		pWinches[ii]->omega = arma::as_scalar(y.row(numSystem-(ii+1)));
		pWinches[ii]->omega = arma::as_scalar(y.row(numSystem-(ii+1)));
	}
	
	// Update BodyBCP positions and velocities
	//std::cout << "Main::fun - Update BCP positions and velocities" << std::endl;
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii]->UpdateBcps();
	}
	// Set boundary conditions on pos and vel of Lines if the BCP is not a joint
	//std::cout << "Main::fun - Set Boundary conditios" << std::endl;
	if (numLines >= 1){
		if (pLines[0]->pLineBcps[0]->tBCP != time){
			for(int ii=0; ii<numLines; ii++)
            {		
				if(pLines[ii]->pLineBcps[0]->GetType() != 3){
					pLines[ii]->pLineBcps[0]->GetValues(time);
				}
				if(pLines[ii]->pLineBcps[1]->GetType() != 3){
					pLines[ii]->pLineBcps[1]->GetValues(time);
				}
			}
		}
	}
	for(int ii=0; ii<numLines; ii++)
    {
		if(pLines[ii]->pLineBcps[0]->GetType() != 3){
			pLines[ii]->pos.row(0)                = pLines[ii]->pLineBcps[0]->pos.t();
			pLines[ii]->vel.row(0)                = pLines[ii]->pLineBcps[0]->vel.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() != 3){
			pLines[ii]->pos.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->pos.t();
			pLines[ii]->vel.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->vel.t();
		}
	}
	// Set boundary conditions on pos and vel of Lines if the BCP is a joint
	//std::cout << "Main::fun - Set Boundary conditios if BCP is a Joint" << std::endl;
	for(int ii=0; ii<numLines; ii++)
    {
		if(pLines[ii]->pLineBcps[0]->GetType() == 3){
			pLines[ii]->pLineBcps[0]->posLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->pos.row(0);
			pLines[ii]->pLineBcps[0]->velLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->vel.row(0);
			pLines[ii]->pLineBcps[0]->iLJ = pLines[ii]->pLineBcps[0]->iLJ + 1;
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->pLineBcps[1]->posLines.row(pLines[ii]->pLineBcps[1]->iLJ) = pLines[ii]->pos.row(pLines[ii]->N-1);
			pLines[ii]->pLineBcps[1]->velLines.row(pLines[ii]->pLineBcps[1]->iLJ) = pLines[ii]->vel.row(pLines[ii]->N-1);
			pLines[ii]->pLineBcps[1]->iLJ = pLines[ii]->pLineBcps[1]->iLJ + 1;
		}
	}
	for(int ii=0; ii<numLines; ii++)
    {
		if(pLines[ii]->pLineBcps[0]->GetType() == 3){
			pLines[ii]->pLineBcps[0]->GetValues(time);
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->pLineBcps[1]->GetValues(time);
		}
	}
	for(int ii=0; ii<numLines; ii++)
    {
		if(pLines[ii]->pLineBcps[0]->GetType() == 3){
			pLines[ii]->pos.row(0)                = pLines[ii]->pLineBcps[0]->pos.t();
			pLines[ii]->vel.row(0)                = pLines[ii]->pLineBcps[0]->vel.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->pos.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->pos.t();
			pLines[ii]->vel.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->vel.t();
		}
	}
	// Compute forces vector for the different Lines
	//std::cout << "Main::fun - Compute forces vector for different lines" << std::endl;
	for(int ii=0; ii<numLines; ii++)
    {
		pLines[ii]->SEM_computeF();
	}
	// Compute forces of Springs
	//std::cout << "Main::fun - Compute spring" << std::endl;
	for(int ii=0; ii<numSprings; ii++)
    {
		pSprings[ii]->computeSpringForces();
	}
	// Compute hydrostatic and hidrodynamic forces
    //std::cout << "Main::fun - Compute hydrodynamic and hydrostatic forces" << std::endl;
    arma::mat Fb = arma::zeros(6*numBodies, 1);
	arma::mat Fh;
	arma::mat Fr;
	arma::mat Fe;
    for (int ii=0; ii<numBodies; ii++)
    {
        Fh = pBodies[ii]->hydro->ComputeHydrostaticForces();
        Fr = pBodies[ii]->hydro->ComputeRadiationForces();
        Fe = pBodies[ii]->hydro->ComputeFirstWaveExcForce();
        Fb(arma::span(6*ii,6*(ii+1)-1), 0) =  Fe - (Fr + Fh);
    }
	// Compute forces on BCPs
	//std::cout << "Main::fun - Compute forces on BCPs" << std::endl;
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii]->ComputeBcpForces();
		Fb(arma::span(6*ii,6*(ii+1)-1), 0) =  Fb(arma::span(6*ii,6*(ii+1)-1), 0) + pBodies[ii]->bcpForces;
	}
	//WriteASCII("output/bcpForces.dat", pBodies[0]->bcpForces, true);
	// Compute body acceleration
	//(**SD.Water->pAddedMassHf).print();
	//arma::mat accB = (*SD.Water->pStructuralMass+**SD.Water->pAddedMassHf)*Fb;
	//std::cout << "Main::fun - Compute accelerations" << std::endl;
	arma::mat accB;
	arma::mat total_mass;
	//WriteASCII("output/accB.dat", accB, true);
	for(int ii=0; ii<numBodies; ii++)
    {
		total_mass = *pBodies[ii]->hydro->pStructuralMass+*pBodies[ii]->hydro->pAddedMassHf[ii];
		//std::cout << "Total mass matrix(2,2): " << mmm(2,2) << std::endl;
		accB = arma::solve(total_mass,Fb);
		//accB = arma::solve(*SD.Water->pStructuralMass,Fb);
		pBodies[ii]->acc = accB(arma::span(6*ii,6*(ii+1)-1), 0);
	}

	// Update BodyBCP accelerations
	//std::cout << "Main::fun - Compute BCP accelerations" << std::endl;
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii]->UpdateBcps();
	}
	//WriteASCII("output/bcpAcceleration.dat", SD.Bodies[0]->pBodyBcps[0]->accG_BCP, true);
	// Obtain Lines accelerations, imposing boundary conditions if the BCP is not a joint
	//std::cout << "Main::fun - Compute lines accelerations" << std::endl;
	for(int ii=0; ii<numLines; ii++)
    {
		if(pLines[ii]->pLineBcps[0]->GetType() != 3){
			pLines[ii]->F.row(0)                = pLines[ii]->pLineBcps[0]->acc.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() != 3){
			pLines[ii]->F.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->acc.t();
		}
		if((pLines[ii]->pLineBcps[0]->GetType() != 3)&&(pLines[ii]->pLineBcps[1]->GetType() != 3)){
			pLines[ii]->acc = pLines[ii]->inv_MM_1N * pLines[ii]->F / pLines[ii]->dL;
		} else if ((pLines[ii]->pLineBcps[0]->GetType() == 3)&&(pLines[ii]->pLineBcps[1]->GetType() != 3)){
			pLines[ii]->acc = pLines[ii]->inv_MM_N * pLines[ii]->F / pLines[ii]->dL;
		} else if ((pLines[ii]->pLineBcps[0]->GetType() != 3)&&(pLines[ii]->pLineBcps[1]->GetType() == 3)){
			pLines[ii]->acc = pLines[ii]->inv_MM_1 * pLines[ii]->F / pLines[ii]->dL;
		} else if ((pLines[ii]->pLineBcps[0]->GetType() == 3)&&(pLines[ii]->pLineBcps[1]->GetType() == 3)){
			pLines[ii]->acc = pLines[ii]->inv_MM * pLines[ii]->F / pLines[ii]->dL;
		}
	}
	// Obtain Lines accelerations, imposing boundary conditions if the BCP is a joint
	//std::cout << "Main::fun - Lines accelerations if the Line is a Joint" << std::endl;
	for(int ii=0;ii<numLines;ii=ii+1){
		if(pLines[ii]->pLineBcps[0]->GetType() == 3){
			pLines[ii]->pLineBcps[0]->accLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->acc.row(0);
			pLines[ii]->pLineBcps[0]->iLJ = pLines[ii]->pLineBcps[0]->iLJ + 1;
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->pLineBcps[1]->accLines.row(pLines[ii]->pLineBcps[1]->iLJ) = pLines[ii]->acc.row(pLines[ii]->N-1);
			pLines[ii]->pLineBcps[1]->iLJ = pLines[ii]->pLineBcps[1]->iLJ + 1;
		}
	}
	for(int ii=0;ii<numLines;ii=ii+1){
		if(pLines[ii]->pLineBcps[0]->GetType() == 3){
			pLines[ii]->pLineBcps[0]->GetValues(time);
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->pLineBcps[1]->GetValues(time);
		}
	}
	for(int ii=0;ii<numLines;ii=ii+1){
		if(pLines[ii]->pLineBcps[0]->GetType() == 3){
			pLines[ii]->acc.row(0)                = pLines[ii]->pLineBcps[0]->acc.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->acc.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->acc.t();
		}
	}

	// Compute Winchies
	//std::cout << "Main::fun - Compute Winchies" << std::endl;
	for(int ii=0;ii<numWinches;ii=ii+1){
		pWinches[ii]->computeWinchie();
	}

	// Copy info from the objects to yprime
	//std::cout << "Main::fun - Copy info to yprime" << std::endl;
	yprime.rows(0,numSystem2-1) = y.rows(numSystem2, numSystem-1);
	yprime.rows(numSystem2,numSystem2+6*numBodies-1) = accB;
	ini = numSystem2+6*numBodies;
	for(int ii=0;ii<numLines;ii=ii+1){
		for(int jj=0;jj<pLines[ii]->N;jj=jj+1){
			yprime.rows(ini,ini+2) = pLines[ii]->acc.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<numWinches;ii=ii+1){
		yprime.row(numSystem2-(ii+1)) = pWinches[ii]->omega;
		yprime.row(numSystem-(ii+1)) = pWinches[ii]->alpha;
	}
	//std::cout << "Main::fun - Check if yprime has a NaN" << std::endl;
	//WriteASCII("output/yprime.dat", yprime, true);
	if (yprime.has_nan()){
		std::cout << std::endl << "ERROR: NaN Detected!" << std::endl;
		throw std::exception();
	}
	//std::cout << "Main::fun - End of fcn" << std::endl;
	return yprime;
}


void Simulation::Initialize()
{
    /**
    timeBufferCount = 3;
    timeBuffer(0, 0) = 0.0;
    timeBuffer(0, 1) = 50.0;
    timeBuffer(0, 2) = 100.0;
    double start_time = 117.0;
    **/
    double start_time = 0.0;

    // Initialize system vector
    std::cout << "Num DOFs Total: " << this->numDofTotal << std::endl;
    std::cout << "Num NumBodies: " << this->numBodies << std::endl;
    std::cout << "Num NumWinchies: " << this->numWinches << std::endl;
    numSystem2 = 3*this->numDofTotal + 6*this->numBodies + this->numWinches;
    numSystem = 2*numSystem2;

    printf("Sistema size: %d\n", numSystem);
    printf("Sistema2 size: %d\n", numSystem2);

    arma::mat y = arma::zeros(numSystem,1);
    arma::mat yprime = arma::zeros(numSystem,1);

    int ini=0;
    std::string file_path;
    if(this->readEquilibrium==0){
        for(int ii=0; ii<this->numBodies; ii=ii+1){
                y.rows(ini,ini+5) = this->pBodies[ii]->pos;
                ini=ini+6;
        }
        for(int ii=0; ii<this->numLines; ii=ii+1){
            for(int jj=0; jj<this->pLines[ii]->N; jj=jj+1){
                    y.rows(ini,ini+2) = this->pLines[ii]->pos.row(jj).t();
                    ini=ini+3;
            }
        }
    }else if (this->readEquilibrium){
        file_path = JoinPath(this->inputFolderPath, "Equilibrio.dat");
        std::ifstream equi(file_path);
            for(int ii=6*this->numBodies;ii<this->numSystem2;ii=ii+1) {
                equi >> y(ii,0);    equi.ignore(std::numeric_limits<int>::max(), '\n');
            }
        equi.close();

        ini=6*this->numBodies;
        for(int ii=0; ii<this->numLines; ii=ii+1){
            for(int jj=0; jj<this->pLines[ii]->N; jj=jj+1){
                this->pLines[ii]->pos.row(jj) = y.rows(ini,ini+2).t();
                ini=ini+3;
            }
        }
    }

     // Initialize Temporal Solver
    if (this->timeIntMethod == 1)
	{
        std::cout << "Initializing temporal solver..." << std::endl;
        pTimeSolver = new BDF(start_time, this->simulationTime, this->maxTimeStep, y, this);
        pTimeSolver->Initialize();
        pTimeSolver->atol = this->timeIntAbsTol;
        pTimeSolver->rtol = this->timeIntRelTol;
        pTimeSolver->nIterMax = this->maxIterStep;
    }
    std::cout << "Simulation::Initialize - Time: " << this->pTimeSolver->t << std::endl;

    // Write initial condition to files
    for(int ii=0; ii<this->numLines; ii=ii+1) this->pLines[ii]->WriteOut(start_time);
    for(int ii=0; ii<this->numBodies; ii=ii+1) this->pBodies[ii]->WriteOut(start_time);

    // Save first data
    std::cout << "Antes de update system" << std::endl;
    this->UpdateSystem();
    
}


void Simulation::LoadCase()
{
    // Read Simulation Properties
    this->ReadProperties();

    // Read Components Data
    this->ReadBodies();
    this->ReadHydrodynamicsHDF5();
    this->ReadLines();
    this->ReadBcps();
    this->ReadWinches();
    this->ReadSprings();

    // Setup case
    this->SetupCase();
}


void Simulation::ReadBcps()
{
    (this->*pReadBcps)();
}


void Simulation::ReadBcpsASCII()
{
    std::cout << "--> Reading BCPs Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosBCPs.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosBCPs.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
    
    // Read data
    fscanf(file_pointer, "%d %[^\n]\n", &numFairBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numAnchorBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numJointBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numBodyBcps, bufferLine);
	numBcps = numFairBcps + numAnchorBcps + numJointBcps + numBodyBcps;

    printf("Number of fairleads: %d\n", numFairBcps);
    printf("Number of anchor: %d\n", numAnchorBcps);
    printf("Number of joint: %d\n", numJointBcps);
    printf("Number of body: %d\n", numBodyBcps);

	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA BCP
	pBcps = new BCP* [numBcps];
	int bcp_count = 0;

    // Se leen los BCPs
	pFairleadBcps = new FairleadBCP* [numFairBcps];
	for(int ii=0; ii<numFairBcps; ii++)
    {
		pFairleadBcps[ii] = new FairleadBCP(bcp_count);
		pFairleadBcps[ii]->ReadPropertiesASCII(file_pointer, inputFolderPath);
		dynamic_cast<FairleadBCP*>(pFairleadBcps[ii])->Initialize(inputFolderPath);
		pBcps[bcp_count] = pFairleadBcps[ii];
		bcp_count++;
	}
	pAnchorBcps = new AnchorBCP* [numAnchorBcps];
	for(int ii=0; ii<numAnchorBcps; ii++)
    {
		pAnchorBcps[ii] = new AnchorBCP(bcp_count);
		pAnchorBcps[ii]->ReadPropertiesASCII(file_pointer);
		pBcps[bcp_count] = pAnchorBcps[ii];
		bcp_count++;
	}
	pJointBcps = new JointBCP* [numJointBcps];
	for(int ii=0; ii<numJointBcps; ii++)
    {
		pJointBcps[ii] = new JointBCP(bcp_count);
		pJointBcps[ii]->ReadPropertiesASCII(file_pointer);
		pBcps[bcp_count] = pJointBcps[ii];
		bcp_count++;
	}
	pBodyBcps = new BodyBCP* [numBodyBcps];
	for(int ii=0; ii<numBodyBcps; ii++)
    {
		pBodyBcps[ii] = new BodyBCP(bcp_count);
		pBodyBcps[ii]->ReadPropertiesASCII(file_pointer);
		pBcps[bcp_count] = pBodyBcps[ii];
		bcp_count++;
	}

    // Loop over BCPs to check if it is necessary to read the winches file
    for (int ii=0; ii<numBcps; ii++)
    {
        if (pBcps[ii]->winchId !=0)
        {
            useWinches = true;
            break;
        }
    }

    // Close file
    fclose(file_pointer);
    std::cout << "----> BCPs Properties Read" << std::endl;
}


void Simulation::ReadBcpsHDF5()
{
    std::cout << "--> Reading BCPs Properties (ASCII format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBcpsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> BCPs Properties Read" << std::endl;
}


void Simulation::ReadBodies()
{
    (this->*pReadBodies)();
}


void Simulation::ReadBodiesASCII()
{
    std::cout << "--> Reading Bodies Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosBodies.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosBodies.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
    
    // Read total number of bodies to read
    fscanf(file_pointer, "%d %[^\n]\n", &numBodies, bufferLine);

	//Read all bodies
    pBodies = new Body* [numBodies];
    printf("Total number of bodies: %d\n", numBodies);
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii] = new Body(ii, this);
		pBodies[ii]->ReadPropertiesASCII(file_pointer);

	}
    /**
	for(int ii=0; ii<nBCPs; ii=ii+1)
    {
		BCPs[ii]->getValues(0.0);
	}
    **/

    fclose(file_pointer);
    std::cout << "----> Bodies Properties Read" << std::endl;
}


void Simulation::ReadBodiesHDF5()
{
    std::cout << "--> Reading Bodies Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBodiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Bodies Properties Read" << std::endl;
}


void Simulation::ReadHydrodynamicsHDF5()
{
    std::cout << "--> Reading Hydrodynamics Properties (HDF5 format)" << std::endl;
    // Read associated hydrodynamics
    std::string file_path = JoinPath(inputFolderPath, "cajon1_LC1.ehydb");
    for (int ii=0; ii<numBodies; ii++)
    {
        pBodies[ii]->hydro = new HydroDatabase(ii, pBodies, this);
        pBodies[ii]->hydro->ReadHydroMechanicsHDF5(file_path);
        pBodies[ii]->hydro->ComputeIRF();
        std::cout << "Structural Mass (2,2): " << (*pBodies[ii]->hydro->pStructuralMass)(2,2) << std::endl;
        std::cout << "Added Mass Hf (2,2): " << (*pBodies[ii]->hydro->pAddedMassHf[ii])(2,2) << std::endl;
        std::cout << "Stiffness (2,2): " << (*pBodies[ii]->hydro->pHydrostaticStiffness)(2,2) << std::endl;
        std::cout << "This is the time vector..." << std::endl;
    }

    // Check time buffere w.r.t IRF size
    if (this->timeBufferSize < 10*pBodies[0]->hydro->numPointsIRF)
    {
        this->timeBufferSize = 10*pBodies[0]->hydro->numPointsIRF;
        this->timeBuffer = arma::zeros(1, this->timeBufferSize);

        for (int ii=0; ii<numBodies; ii++)
        {
            pBodies[ii]->velBufferSize = 10*pBodies[ii]->hydro->numPointsIRF;
            pBodies[ii]->velBuffer = arma::zeros(6, pBodies[ii]->velBufferSize);
        }
    }
    
    std::cout << "----> Hydrodynamic Properties Read" << std::endl;
}


void Simulation::ReadLines()
{
    (this->*pReadLines)();
}


void Simulation::ReadLinesASCII()
{
    std::cout << "--> Reading Lines Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
	std::string file_path = JoinPath(inputFolderPath, "datosLines.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosLines.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

   // Read total number of springs to read 
    fscanf(file_pointer, "%d %[^\n]\n", &numLines, bufferLine);

	pLines = new Line*[numLines];
	//INICIO LAS LINEAS
	for(int ii=0; ii<numLines; ii++)
    {
		pLines[ii] = new Line(ii, gravity, waterDensity, waterDepth);
		try
		{
            pLines[ii]->ReadPropertiesASCII(file_pointer);
            //pLines[ii]->print_out();
            /**
            pLines[ii]->LineBCP[0] = BCPs[pLines[ii]->BCP_1-1];
            pLines[ii]->LineBCP[1] = BCPs[pLines[ii]->BCP_N-1];
            pLines[ii]->pos_1 = BCPs[pLines[ii].BCP_1-1]->pos;
            pLines[ii]->pos_N = BCPs[pLines[ii].BCP_N-1]->pos;
            nNodosTotal=nNodosTotal+pLines[ii].N;
            if(flag_read_eq==0) pLines[ii]->initLine();
            Lines[ii].SEM_getBaseFunctions();
            **/
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is under the floor level." << std::endl << std::endl;
			if (e==1) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor and it shouldn't. " << std::endl << std::endl;
			if (e==2) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl << std::endl;
			if (e==3) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and vertical. It should be pretensed. " << std::endl << std::endl;
			if (e==4) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and it should. " << std::endl << std::endl;
			if (e==5) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " initial shape can't be computed with QS method. " << std::endl;
			if (e==6) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor althoug none of its ends are there. " << std::endl << std::endl;
		}
        
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Lines Properties Read" << std::endl;
}


void Simulation::ReadLinesHDF5()
{
    std::cout << "--> Reading Lines Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadLinesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Lines Properties Read" << std::endl;
}


void Simulation::ReadSprings()
{
    (this->*pReadSprings)();
}


void Simulation::ReadSpringsASCII()
{
    std::cout << "--> Reading Spring Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
	std::string file_path = JoinPath(inputFolderPath, "datosSprings.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosSprings.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

    // Read file contents
    fscanf(file_pointer, "%d %[^\n]", &numSprings, bufferLine);

	//ALOCATO UN VECTOR DE POINTERS A OBJETOS, UNO PARA CADA MUELLE
	pSprings = new Spring*[numSprings];

	//INICIO LLOS MUELLES
	for(int ii=0; ii<numSprings; ii++)
    {
		pSprings[ii] = new Spring(ii);
		pSprings[ii]->ReadPropertiesASCII(file_path);
		//pSprings[ii].SpringBCP[0] = BCPs[Springs[ii].BCP_1-1];
		//pSprings[ii].SpringBCP[1] = BCPs[Springs[ii].BCP_2-1];
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Spring Properties Read" << std::endl;
}


void Simulation::ReadSpringsHDF5()
{
    std::cout << "--> Reading Spring Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadSpringsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Spring Properties Read" << std::endl;
}


void Simulation::ReadProperties()
{
    (this->*pReadProperties)();
}


void Simulation::ReadPropertiesASCII()
{
    std::cout << "--> Reading Simulation Properties (ASCII format)" << std::endl;
    std::string file_path = JoinPath(inputFolderPath, "datosProblema.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosProblema.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
	
    char bufferLine [1000];
    int dummyBool;
    fscanf(file_pointer, "%lf %[^\n]\n", &gravity, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &waterDensity, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &waterDepth, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &maxTimeStep, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &simulationTime, bufferLine);
	fscanf(file_pointer, "%d %[^\n]\n", &timeIntMethod, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &timeIntAbsTol, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &timeIntRelTol, bufferLine);
	fscanf(file_pointer, "%d %[^\n]\n", &maxIterStep, bufferLine);
	fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    readEquilibrium = dummyBool;
	fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    writeEquilibrium = dummyBool;

    // Close file
    fclose(file_pointer);


    // Show inputs
    if (true)
    {
        std::cout << "Gravity: " << gravity << std::endl;
        std::cout << "Water Density: " << waterDensity << std::endl;
        std::cout << "Water Depth: " << waterDepth << std::endl;
        std::cout << "Max Time Step: " << maxTimeStep << std::endl;
        std::cout << "Simulation Time: " << simulationTime << std::endl;
        std::cout << "Time Integration Method: " << timeIntMethod << std::endl;
        std::cout << "Time Integration Absolute Tolerace: " << timeIntAbsTol << std::endl;
        std::cout << "Time Integration Relative Tolerace: " << timeIntRelTol << std::endl;
        std::cout << "Max Iterations per Step: " << maxIterStep << std::endl;
        std::cout << "Read Equilibrium: " << readEquilibrium << std::endl;
        std::cout << "Write Equilibrium: " << writeEquilibrium << std::endl;
    }

    std::cout << "----> Simulation Properties Read" << std::endl;
}


void Simulation::ReadPropertiesHDF5()
{
    std::cout << "--> Reading Simulation Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadPropertiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Simulation Properties Read" << std::endl;
}


void Simulation::ReadWinches()
{
    (this->*pReadWinches)();
}


void Simulation::ReadWinchesASCII()
{
    std::cout << "--> Reading Winches Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
	std::string file_path = JoinPath(inputFolderPath, "datosWinchies.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosWinchies.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

    // Read number of winches defined in the file
	fscanf(file_pointer, "%d %[^\n]\n", &numWinches, bufferLine);

    if ((numWinches ==0) && useWinches)
    {
        throw ValueError("Use of winches is requested when loading BCPs but there is no winches specified in datosWinches.dat\n");
    }

	// Allocate a vector of pointers to Winch class objects
	pWinches = new Winchie*[numWinches];

	// Read Winches
	for(int ii=0; ii<numWinches; ii++)
    {
		pWinches[ii] = new Winchie(ii);
		pWinches[ii]->ReadPropertiesASCII(file_pointer);
		pWinches[ii]->LineW = pLines[pWinches[ii]->nLine - 1];
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Winches Properties Read" << std::endl;
}


void Simulation::ReadWinchesHDF5()
{
    std::cout << "--> Reading Winches Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWinchesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Winches Properties Read" << std::endl;
}


void Simulation::Run()
{
    time_t tstart, tend;
	double wallTime = 0.0;
    tstart = time(0);
    std::cout<< "    t = " << wallTime << " s" << std::endl;
    do
    {
        pTimeSolver->step();
        
        UpdateSystem();
        
        // Print out time if any
        std::cout<< "    t = " << pTimeSolver->t << " s"  << std::endl;
        for(int ii=0; ii<numLines; ii=ii+1) pLines[ii]->WriteOut(pTimeSolver->t);
        for(int ii=0; ii<numBodies; ii=ii+1) pBodies[ii]->WriteOut(pTimeSolver->t);
        /**
        if (pTimeSolver->t >= wallTime + maxTimeStep)
        {
            //wallTime = wallTime + maxTimeStep;
            std::cout<< "    t = " << wallTime << " s"  << std::endl;
            //if (nWinchies>0) CW.controlWinchies();
            for(int ii=0; ii<numLines; ii=ii+1) pLines[ii]->WriteOut(wallTime);
            for(int ii=0; ii<numBodies; ii=ii+1) pBodies[ii]->WriteOut(wallTime);
        }
        **/
    } while (pTimeSolver->t <= simulationTime);
    
    tend = time(0); 
    std::cout << std::endl << "    Computational time  : " << difftime(tend, tstart) << " seconds" << std::endl;
    std::cout << "    Total function calls: " << numCallsSysFun << std::endl;
    std::cout << "    Total jac calls: " << pTimeSolver->iJ << std::endl << std::endl;

    /**
    if (writeEquilibrium == 1) {

        std::cout << "  Writting data to Equilibrio.dat ..." << std::endl << std::endl; //////////////////////////////////////////

        std::ofstream equi("output/Equilibrio.dat");
            for(int ii=6*nBodies;ii<nSistema2;ii=ii+1) equi << y(ii,0) << std::endl;
        equi.close();
    }
    **/
}

void Simulation::SetupCase()
{
    std::cout << "--> Setting up the case configuration..." << std::endl;
    // Count the number of BCP in each body and create pointer array
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<pBodies[ii]->numBcps; jj++)
        {
            if (pBodies[ii]->pIndexBcps[jj]+1 > numBcps)
            {
                std::stringstream ss;
                ss << "BCP index: " << pBodies[ii]->pIndexBcps[jj] << " in Body: " << pBodies[ii]->GetId() \
                    << " is out of range when compare with the Number of BCPs(" << numBcps << ") defined in" \
                    << " datosBCPs.dat";
                throw ValueError(ss.str());
            }
            pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp++;
        }
    }
    bool defined_body_bcps [numBcps] = {0}; 
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<pBodies[ii]->numBcps; jj++)
        {
            if (!defined_body_bcps[pBodies[ii]->pIndexBcps[jj]])
            {
                pBcps[pBodies[ii]->pIndexBcps[jj]]->pBodies = new Body* [pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp];
                defined_body_bcps[pBodies[ii]->pIndexBcps[jj]] = true;
            }
        }
    }

    // Assing to each BCP the corresponding Body pointer
    for (int ii=0; ii<numBodies; ii++)
    {
        for(int jj=0; jj<pBodies[ii]->numBcps; jj++)
        {
            pBodies[ii]->pBodyBcps[jj] = pBcps[pBodies[ii]->pIndexBcps[jj]];
            pBodies[ii]->pBodyBcps[jj]->pBodies[pBodies[ii]->pBodyBcps[jj]->countBody] = pBodies[ii];
            pBodies[ii]->pBodyBcps[jj]->countBody++;
        }
        pBodies[ii]->UpdateBcps();
    }

    // Check every thing is correct
    /**
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<pBodies[ii]->numBcps; jj++)
        {
            std::cout << "Body ID: " << pBodies[ii]->GetId() << " - BCP ID: " << pBodies[ii]->pBodyBcps[jj]->GetId() << std::endl;
        }
    }
    
    for (int ii=0; ii<numBcps; ii++)
    {
        for (int jj=0; jj<pBcps[ii]->numBodiesBcp; jj++)
        {
            std::cout << "BCP ID: " << pBcps[ii]->GetId() << " - Body ID: " <<  pBcps[ii]->pBodies[jj]->GetId() << std::endl;
            std::cout << "BCP ID: " << pBcps[ii]->GetId() << " - X pos: " <<  pBcps[ii]->pos[0] << std::endl;
            std::cout << "BCP ID: " << pBcps[ii]->GetId() << " - Y pos: " <<  pBcps[ii]->pos[1] << std::endl;
            std::cout << "BCP ID: " << pBcps[ii]->GetId() << " - Z pos: " <<  pBcps[ii]->pos[2] << std::endl;
        }
        std::cout << "heree2" << std::endl;
    }
    **/


    // Count the number of Lines in each body and create pointer array
    for (int ii=0; ii<numLines; ii++)
    {
        for (int jj=0; jj<pLines[ii]->numBcps; jj++)
        {
            if (pLines[ii]->indexBcps[jj]+1 > numBcps)
            {
                std::stringstream ss;
                ss << "BCP index: " << pLines[ii]->indexBcps[jj] << " in Line: " << pLines[ii]->GetId() \
                    << " is out of range when compare with the Number of BCPs(" << numBcps << ") defined in" \
                    << " datosBCPs.dat";
                throw ValueError(ss.str());
            }
            pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp++;
        }
    }
    bool defined_lines_bcps [numBcps] = {0}; 
    for (int ii=0; ii<numLines; ii++)
    {
        for (int jj=0; jj<pLines[ii]->numBcps; jj++)
        {
            if (!defined_lines_bcps[pLines[ii]->indexBcps[jj]])
            {
                pBcps[pLines[ii]->indexBcps[jj]]->pLines = new Line* [pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                defined_lines_bcps[pLines[ii]->indexBcps[jj]] = true;
            }
        }
        
    }
    // Assing to each Line the corresponding Body pointer
    for (int ii=0; ii<numLines; ii++)
    {
        for(int jj=0; jj<pLines[ii]->numBcps; jj++)
        {
            pLines[ii]->pLineBcps[jj] = pBcps[pLines[ii]->indexBcps[jj]];
            pLines[ii]->pLineBcps[jj]->pLines[pLines[ii]->pLineBcps[jj]->countLine] = pLines[ii];
            pLines[ii]->pLineBcps[jj]->countLine++;
        }
        numDofTotal += pLines[ii]->N;
        try
        {
            if(!readEquilibrium) pLines[ii]->initLine();
            pLines[ii]->print_out();
            pLines[ii]->SEM_getBaseFunctions();
        }
        catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is under the floor level." << std::endl << std::endl;
			if (e==1) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor and it shouldn't. " << std::endl << std::endl;
			if (e==2) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl << std::endl;
			if (e==3) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and vertical. It should be pretensed. " << std::endl << std::endl;
			if (e==4) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " is not tense and it should. " << std::endl << std::endl;
			if (e==5) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " initial shape can't be computed with QS method. " << std::endl;
			if (e==6) std::cout<< "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor althoug none of its ends are there. " << std::endl << std::endl;
		}
    }
    // Check every thing is correct
    /**
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<pLines[ii]->numBcps; jj++)
        {
            std::cout << "Line ID: " << pLines[ii]->GetId() << " - BCP ID: " << pLines[ii]->pLineBcps[jj]->GetId() << std::endl;
        }
    }
    for (int ii=0; ii<numBcps; ii++)
    {
        for (int jj=0; jj<pBcps[ii]->numLinesBcp; jj++)
        {
            std::cout << "BCP ID: " << pBcps[ii]->GetId() << " - Line ID: " <<  pBcps[ii]->pLines[jj]->GetId() << std::endl;
        }
    }
    **/

    // Setup Springs
	for(int ii=0; ii<numSprings; ii++)
    {
		pSprings[ii]->SpringBCP[0] = pBcps[pSprings[ii]->BCP_1];
		pSprings[ii]->SpringBCP[1] = pBcps[pSprings[ii]->BCP_2];
	}
    std::cout << "----> Case configuration done" << std::endl;
}


Simulation::Simulation(std::string incProjectPath, std::string incDataFormat)
{
    // Assign direct variables
    projectFolderPath = incProjectPath;

    // Process indirect variables
    if (!incDataFormat.compare("ASCII"))
    {
        dataFormat = 0;
        dataFormatStr = "ASCII";
        inputFolderPath = JoinPath(incProjectPath, "input");
        outputFolderPath = JoinPath(incProjectPath, "output");
        pReadProperties = &(Simulation::ReadPropertiesASCII);
        pReadBcps = &(Simulation::ReadBcpsASCII);
        pReadBodies = &(Simulation::ReadBodiesASCII);
        pReadLines = &(Simulation::ReadLinesASCII);
        pReadSprings = &(Simulation::ReadSpringsASCII);
        pReadWinches = &(Simulation::ReadWinchesASCII);
    }
    else if (!incDataFormat.compare("HDF5"))
    {
        dataFormat = 1;
        dataFormatStr = "HDF5";
        inputFolderPath = incProjectPath;
        outputFolderPath = incProjectPath;
        pReadProperties = &(Simulation::ReadPropertiesHDF5);
        pReadBcps = &(Simulation::ReadBcpsHDF5);
        pReadBodies = &(Simulation::ReadBodiesHDF5);
        pReadLines = &(Simulation::ReadLinesHDF5);
        pReadSprings = &(Simulation::ReadSpringsHDF5);
        pReadWinches = &(Simulation::ReadWinchesHDF5);
    }
    else
    {
        std::stringstream ss;
        ss << "Simulation data format --> " << incDataFormat << " is not available.\n    Available formats: ASCII | HDF5.";
        throw ValueError(ss.str());
    }

    // 

}


void Simulation::UpdateSystem()
{
    // Update time vector if any
    timeBufferCount++;
    if (timeBufferCount < timeBufferSize)
    {
        timeBuffer(0, timeBufferCount) = pTimeSolver->t;
    }
    else
    {
        arma::mat timeBufferNew = arma::zeros(1, timeBufferSize);
        timeBufferNew.cols(0, pBodies[0]->hydro->numPointsIRF-1) = timeBuffer.cols(timeBufferSize-pBodies[0]->hydro->numPointsIRF, timeBufferSize-1);
        timeBuffer = timeBufferNew;
        timeBufferCount = pBodies[0]->hydro->numPointsIRF-1;
    }

    // Update bodies velocity
    int ini = 0;
    for(int ii=0; ii<numBodies; ii++)
	{
		pBodies[ii]->StoreVelocities();
		ini = ini + 6;
	}

    // Update hydrodynamic properties
    if (timeBufferCount > 0)
    {
        for(int ii=0; ii<numBodies; ii++)
        {
            pBodies[ii]->UpdateHydrostaticForces();
            pBodies[ii]->UpdateRadiationForces();
        }
    }
    
}