
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <ctime>
#include "Simulation.hpp"
#include "../CommonTools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../BCPs/Winchies.hpp"
#include "../BCPs/WinchiesController.hpp"
#include "../Waves/Wave.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"
#include "../WindTurbine/WindTurbine.hpp"

#ifndef __has_include
  static_assert(false, "__has_include not supported");
#else
#  if __cplusplus >= 201703L && __has_include(<filesystem>)
#    include <filesystem>
     namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
     namespace fs = std::experimental::filesystem;
#  elif __has_include(<boost/filesystem.hpp>)
#    include <boost/filesystem.hpp>
     namespace fs = boost::filesystem;
#  endif
#endif


arma::mat Simulation::CalculateSystemDynamics(double time, arma::mat y)
{
    // std::cout << "Time: " << time << " s\n";
    numCallsSysFun++;
    // std::cout << "Simulation::CalculateSystemDynamics - At first" << std::endl; 
	arma::mat yprime = arma::zeros(size(y));
	int i0;
	// std::cout << "Simulation::CalculateSystemDynamics - Before copy y to objects" << std::endl;
	// Copy info from y to the objects.
	int ini = 0;
	for(int ii=0; ii<numBodiesFree; ii++)
	{
		pBodiesFree[ii]->pos = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0; ii<numLines;ii++)
	{
		for(int jj=pLines[ii]->first_node; jj<pLines[ii]->last_node; jj=jj+1)
        {
			pLines[ii]->pos.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0; ii<numWinches; ii++)
    {
		pWinches[ii]->theta = arma::as_scalar(y.row(ini));
	    ini = ini + 1;
	}
    for(int ii=0; ii<numWindTurbines; ii++)
    {
		pWindTurbines[ii]->rotPos = arma::as_scalar(y.row(ini));
	    ini = ini + 1;
	}
	for(int ii=0; ii<numBodiesFree; ii++)
	{
		pBodiesFree[ii]->vel = y.rows(ini,ini+5);
		ini = ini + 6;
	}
	for(int ii=0; ii<numLines; ii++)
	{
		for(int jj=pLines[ii]->first_node; jj<pLines[ii]->last_node; jj=jj+1)
		{
			pLines[ii]->vel.row(jj) = y.rows(ini,ini+2).t();
			ini = ini + 3;
		}
	}
	for(int ii=0; ii<numWinches; ii++)
    {
		pWinches[ii]->omega = arma::as_scalar(y.row(ini));
	    ini = ini + 1;
	}
    for(int ii=0; ii<numWindTurbines; ii++)
    {
		pWindTurbines[ii]->rotSpeed = arma::as_scalar(y.row(ini));
	    ini = ini + 1;
	}

    for(int ii=0; ii<numBodiesLock; ii++)
	{
		pBodiesLock[ii]->UpdateLockBody(time);
	}
	
	// Update BodyBCP positions and velocities
	// std::cout << "Simulation::CalculateSystemDynamics - Update BCP positions and velocities" << std::endl;
	for(int ii=0; ii<numBodies; ii++)
    {
		pBodies[ii]->UpdateBcps();
		pBodies[ii]->ResetBcps();
	}
	// Set boundary conditions on pos and vel of Lines if the BCP is not a joint
	// std::cout << "Simulation::CalculateSystemDynamics - Set Boundary conditios" << std::endl;
	for(int ii=0; ii<numLines; ii++)
    {
		if(pLines[ii]->pLineBcps[0]->GetType() != 3){
			pLines[ii]->pLineBcps[0]->GetValues(time);
			pLines[ii]->pos.row(0)                = pLines[ii]->pLineBcps[0]->pos.t();
			pLines[ii]->vel.row(0)                = pLines[ii]->pLineBcps[0]->vel.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() != 3){
			pLines[ii]->pLineBcps[1]->GetValues(time);
			pLines[ii]->pos.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->pos.t();
			pLines[ii]->vel.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->vel.t();
		}
	}	
	// Set boundary conditions on pos and vel of Lines if the BCP is a joint
	// std::cout << "Simulation::CalculateSystemDynamics - Set Boundary conditios if BCP is a Joint" << std::endl;
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
			pLines[ii]->pos.row(0)                = pLines[ii]->pLineBcps[0]->pos.t();
			pLines[ii]->vel.row(0)                = pLines[ii]->pLineBcps[0]->vel.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 3){
			pLines[ii]->pLineBcps[1]->GetValues(time);
			pLines[ii]->pos.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->pos.t();
			pLines[ii]->vel.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->vel.t();
		}
	}

	// Compute forces vector for the different Lines
	// std::cout << "Simulation::CalculateSystemDynamics - Compute forces vector for different lines" << std::endl;
	for(int ii=0; ii<numLines; ii++)
    {
		pLines[ii]->SEM_computeF();
	}

	// Compute forces of Springs
	// std::cout << "Simulation::CalculateSystemDynamics - Compute spring" << std::endl;
	for(int ii=0; ii<numSprings; ii++)
    {
		pSprings[ii]->computeSpringForces();
	}

	// Update hydrostatic parameters if there is sinking
    // std::cout << "Simulation::CalculateSystemDynamics - Update Sinking Hydrostatics" << std::endl;
	for(int ii=0; ii<numSinking; ii=ii+1) {
    	pSinking[ii]->UpdateSinkingHydrostatics(time);
    }

	// Compute hydrostatic and hidrodynamic forces
    // std::cout << "Simulation::CalculateSystemDynamics - Compute hydrodynamic and hydrostatic forces" << std::endl;
    arma::mat Fb = arma::zeros(6*numBodiesFree, 1);
    for (int ii=0; ii<numBodiesFree; ii++)
    {    
		Fb(arma::span(6*ii,6*(ii+1)-1), 0) =  pBodiesFree[ii]->Fb + pBodiesFree[ii]->pHydro->CalculateHydrostaticForces();
    }

	// Compute forces on BCPs
	// std::cout << "Simulation::CalculateSystemDynamics - Compute forces on BCPs" << std::endl;
	for(int ii=0; ii<numBodiesFree; ii++)
    {
		pBodiesFree[ii]->ComputeBcpForces();
		Fb(arma::span(6*ii,6*(ii+1)-1), 0) =  Fb(arma::span(6*ii,6*(ii+1)-1), 0) + pBodiesFree[ii]->bcpForces;
	}

    // Add wind turbine forces
	// std::cout << "Simulation::CalculateSystemDynamics - Add wind turbine forces" << std::endl;
	for(int ii=0; ii<numBodiesFree; ii++)
    {
        pBodiesFree[ii]->ComputeWindTurbForces();
		Fb(arma::span(6*ii,6*(ii+1)-1), 0) =  Fb(arma::span(6*ii,6*(ii+1)-1), 0) + pBodiesFree[ii]->windTurbForces;
	}

    // Compute also everything for locked bodies so it can be displayed on the output files
    arma::mat dummy;
    for(int ii=0; ii<numBodiesLock; ii++)
	{
        dummy = pBodiesLock[ii]->pHydro->CalculateHydrostaticForces();
		pBodiesLock[ii]->ComputeWindTurbForces();
        pBodiesLock[ii]->ComputeBcpForces();
	}

	// Compute body acceleration
    // std::cout << "Simulation::CalculateSystemDynamics - Compute Bodies accelerations" << std::endl;
    arma::mat accB;
    if (numBodiesFree>0){
        accB = (*pSystemMatrixInv) * Fb;
    }
    for (int ii=0; ii<numBodiesFree; ii++)
    {
        pBodiesFree[ii]->acc = accB(arma::span(6*ii,6*(ii+1)-1), 0)%pBodiesFree[ii]->isDofActive; //////////////////////////////////////////////////////////////////////////////////////////  HARCODEO
    }

	// Update BodyBCP accelerations
	// std::cout << "Simulation::CalculateSystemDynamics - Compute BCP accelerations" << std::endl;
	for(int ii=0; ii<numBodiesFree; ii++)
    {
		pBodiesFree[ii]->UpdateBcps();
	}

	// Obtain Lines accelerations, imposing boundary conditions if the BCP is not a joint
	// std::cout << "Simulation::CalculateSystemDynamics - Compute lines accelerations" << std::endl;
	for(int ii=0; ii<numLines; ii++)
    {
    	// For body BCPs, the accelertion has changed, so GetValues() routine is called again
    	if(pLines[ii]->pLineBcps[0]->GetType() == 4){
					pLines[ii]->pLineBcps[0]->GetValues(time);
		}
		if(pLines[ii]->pLineBcps[1]->GetType() == 4){
					pLines[ii]->pLineBcps[1]->GetValues(time);
		}
		// Impose BCP accelerations on lines forces vectors
		if(pLines[ii]->pLineBcps[0]->GetType() != 3){
			pLines[ii]->F.row(0)                = pLines[ii]->pLineBcps[0]->acc.t();
		}
		if(pLines[ii]->pLineBcps[1]->GetType() != 3){
			pLines[ii]->F.row(pLines[ii]->N-1) = pLines[ii]->pLineBcps[1]->acc.t();
		}
	}

    // std::cout << "Simulation::CalculateSystemDynamics - Assemble lines forces vectors" << std::endl;
    arma::mat LinesCouplingVector = arma::zeros(numAllLinesNodes,3);
    for(int ii=0; ii<numLines; ii++)
    {
        LinesCouplingVector.rows(pLines[ii]->ind4CouplingMat) += pLines[ii]->F;
    }
    for(int ii=0; ii<numBcps; ii++){
        if((pBcps[ii]->GetType() == 3) && (pBcps[ii]->flag_assigned == 1)){
            LinesCouplingVector.row(pBcps[ii]->couplingMatIndex) += pBcps[ii]->JointForce;
        }
    }

    // std::cout << "Simulation::CalculateSystemDynamics - Solve lines accelerations" << std::endl;
    arma::mat LinesAccelerations;
    if(numAllLinesNodes>=100){
        arma::superlu_opts opts; opts.allow_ugly  = false;
        arma::spsolve(LinesAccelerations,*pLinesCouplingMatrix_sp,LinesCouplingVector,"superlu",opts);
    } else {
        if (useWinches){
            LinesAccelerations = arma::solve(*pLinesCouplingMatrix,LinesCouplingVector);
        } else {
            LinesAccelerations = (*pLinesCouplingMatrixInv) * LinesCouplingVector;
        }
    }
    for(int ii=0; ii<numLines; ii++)
    {
        pLines[ii]->acc = LinesAccelerations.rows(pLines[ii]->ind4CouplingMat);
    }

	// Compute Winchies
	// std::cout << "Simulation::CalculateSystemDynamics - Compute Winchies" << std::endl;
	for(int ii=0;ii<numWinches;ii=ii+1){
		pWinches[ii]->computeWinchie();
	}
    // Recompute Lines coupling matrix for new dL values
    if(useWinches && numJointBcps>0){
        ComputeLinesCouplingMatrix();
    }

    // Compute Wind Turbines rotor acceleration
	// std::cout << "Simulation::CalculateSystemDynamics - Compute Wind Turbines" << std::endl;
	for(int ii=0;ii<numWindTurbines;ii=ii+1){
		pWindTurbines[ii]->ComputeRotorAcc();
	}

	// Copy info from the objects to yprime
	// std::cout << "Simulation::CalculateSystemDynamics - Copy info to yprime" << std::endl;
	ini = 0;
	for(int ii=0;ii<numBodiesFree;ii=ii+1){
		yprime.rows(ini,ini+5) = pBodiesFree[ii]->vel;
		ini = ini + 6;
	}
	for(int ii=0;ii<numLines;ii=ii+1){
		for(int jj=pLines[ii]->first_node; jj<pLines[ii]->last_node; jj=jj+1){
			yprime.rows(ini,ini+2) = pLines[ii]->vel.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<numWinches;ii=ii+1){
		yprime.row(ini) = pWinches[ii]->omega;
        ini = ini + 1;
	}
    for(int ii=0;ii<numWindTurbines;ii=ii+1){
		yprime.row(ini) = pWindTurbines[ii]->rotSpeed;
        ini = ini + 1;
	}
	for(int ii=0;ii<numBodiesFree;ii=ii+1){
		yprime.rows(ini,ini+5) = pBodiesFree[ii]->acc;
		ini = ini + 6;
	}
	for(int ii=0;ii<numLines;ii=ii+1){
		for(int jj=pLines[ii]->first_node; jj<pLines[ii]->last_node; jj=jj+1){
			yprime.rows(ini,ini+2) = pLines[ii]->acc.row(jj).t();
			ini = ini + 3;
		}
	}
	for(int ii=0;ii<numWinches;ii=ii+1){
		yprime.row(ini) = pWinches[ii]->alpha;
        ini = ini + 1;
	}
    for(int ii=0;ii<numWindTurbines;ii=ii+1){
		yprime.row(ini) = pWindTurbines[ii]->rotAcc;
        ini = ini + 1;
	}

	// std::cout << "Simulation::CalculateSystemDynamics - Check if yprime has a NaN" << std::endl;
	if (yprime.has_nan() | yprime.has_inf()){
		std::cout << std::endl << "ERROR: NaN or Inf Detected! yprime = " << std::endl;
		std::cout << yprime << std::endl;
		throw std::exception();
	}
	
	// std::cout << "Simulation::CalculateSystemDynamics - End of fcn" << std::endl;
	return yprime;
}


void Simulation::CloseCase()
{

    for (int ii=0; ii<numLines; ii++)
    {
    	pLines[ii]->CloseOutputFilesASCII();
    }

    for (int ii=0; ii<numBodies; ii++)
    {
    	pBodies[ii]->CloseOutputFilesASCII();
    }

    for (int ii=0; ii<numSinking; ii++)
    {
    	pSinking[ii]->CloseOutputFilesASCII();
    }

    if (useWinches) {
    	WinchesController.CloseOutputFilesASCII();
	}

	for (int ii=0; ii<numWindTurbines; ii++)
    {
    	pWindTurbines[ii]->Finalize();
    }

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
    std::cout << "Num. Bodies Free: " << this->numBodiesFree << std::endl;
    std::cout << "Num. Mooring DOFs Total: " << this->numDofTotal << std::endl;
    std::cout << "Num. Winchies: " << this->numWinches << std::endl;
    std::cout << "Num. Wind Turbines: " << this->numWindTurbines << std::endl;
    numSystem2 = 3*this->numDofTotal + 6*this->numBodiesFree + this->numWinches + this->numWindTurbines;
    numSystem = 2*numSystem2;

    printf("Sistem size: %d\n", numSystem);

    arma::mat y = arma::zeros(numSystem,1);
    arma::mat yprime = arma::zeros(numSystem,1);

    int ini=0;
    std::string file_path;
    std::cout << "Initiallizing system vector..." << std::endl;
    if(this->readEquilibrium==0){
        for(int ii=0; ii<this->numBodiesFree; ii=ii+1){
                std::cout << "  ... Including Body: " << ii+1 << std::endl;
                y.rows(ini,ini+5) = this->pBodiesFree[ii]->pos;
                ini=ini+6;
        }
        for(int ii=0; ii<this->numLines; ii=ii+1){
            std::cout << "  ... Including Line: " << ii+1 << std::endl;
            for(int jj=pLines[ii]->first_node; jj<pLines[ii]->last_node; jj=jj+1){
                    y.rows(ini,ini+2) = this->pLines[ii]->pos.row(jj).t();
                    ini=ini+3;
            }
        }
        for(int ii=0; ii<this->numWinches; ii=ii+1){
            std::cout << "  ... Including Winch: " << ii+1 << std::endl;
            ini=ini+1;
        }
        for(int ii=0; ii<this->numWindTurbines; ii=ii+1){
            std::cout << "  ... Including Wind Turbine: " << ii+1 << std::endl;
            y(ini) = this->pWindTurbines[ii]->rotPos;
            y(numSystem2+ini) = this->pWindTurbines[ii]->rotSpeed;
            ini=ini+1;
        }
    } else {
        std::cout << "Reading Equilibrium.dat..." << std::endl;
        file_path = JoinPath(inputFolderPath, "Equilibrio.dat");
        y.load(file_path,arma::arma_ascii);

        ini=6*this->numBodiesFree;
        for(int ii=0; ii<this->numLines; ii=ii+1){
            for(int jj=this->pLines[ii]->first_node; jj<this->pLines[ii]->last_node; jj=jj+1){
                this->pLines[ii]->pos.row(jj) = y.rows(ini,ini+2).t();
                ini=ini+3;
            }
            if(this->pLines[ii]->first_node==0){
                this->pLines[ii]->pLineBcps[0]->pos = this->pLines[ii]->pos.row(0).t();
            } else {
                this->pLines[ii]->pos.row(0) = this->pLines[ii]->pLineBcps[0]->pos.t();
            }
            if(this->pLines[ii]->last_node==this->pLines[ii]->N){
                this->pLines[ii]->pLineBcps[1]->pos = this->pLines[ii]->pos.row(this->pLines[ii]->N-1).t();
            } else {
                this->pLines[ii]->pos.row(this->pLines[ii]->N-1) = this->pLines[ii]->pLineBcps[1]->pos.t();
            }
        }
    }

     // Initialize Temporal Solver
    if (this->timeIntMethod == 1)
	{
        std::cout << "Initializing temporal solver..." << std::endl;
        pTimeSolver = new BDF(start_time, this->simulationTime, this->maxTimeStep, y, this);
        std::cout << "  Temporal solver constructor done!" << std::endl;
        pTimeSolver->Initialize();
        pTimeSolver->atol = this->timeIntAbsTol;
        pTimeSolver->rtol = this->timeIntRelTol;
        pTimeSolver->nIterMax = this->maxIterStep;
        std::cout << "  Temporal solver initiallized!" << std::endl;
    }

    // Write initial condition to files
    for(int ii=0; ii<this->numLines; ii=ii+1) this->pLines[ii]->WriteOut(start_time);
    for(int ii=0; ii<this->numBodies; ii=ii+1) this->pBodies[ii]->WriteOut(start_time);
    std::cout << "Update system" << std::endl;
    // Save first data
    this->UpdateSystem();
    std::cout << "System updated!" << std::endl;
}


void Simulation::LoadCase()
{

    fs::remove_all(outputFolderPath);
    fs::create_directory(outputFolderPath);
    
    // Read Simulation Properties
    this->ReadProperties();

    // Read Components Data
    this->ReadBodies();
    this->ReadWaves();
    this->ReadLines();
    this->ReadBcps();
    if (useWinches)
    {
        this->ReadWinches();
    }
    this->ReadSprings();
    this->ReadSinking();
    this->ReadWindTurbines();

    // Setup case
    this->SetupCase();

    this->PrintSetup();
}


void Simulation::PrintSetup(void)
{

    // Print BCP properties
    for (int i=0; i<this->numBcps; i++)
    {
        pBcps[i]->Print();
    }
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
        dynamic_cast<JointBCP*>(pJointBcps[ii])->Initialize(this->gravity,this->waterDensity,this->waterDepth);
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
    int body_count = 0;
    char buffer_line [1000];
    int diff_count = 0;
    int hydro_database_count=0;
    std::string hydro_databases_name [300];
    int max_num_bodies_database=100;
    int pos_body = 0;
    int pos_database=0;

    // Parse file in order to guess the number of bodies
    std::string file_path = JoinPath(inputFolderPath, "datosBodies.dat");
    this->numBodies = parse_file(file_path);

    if (numBodies>0){

    // Open file
    FILE* pFile = fopen(file_path.c_str(), "r");
	
	if (pFile == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosBodies.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

	//Read all bodies
    pBodies = new Body* [numBodies];
	for(int ii=0; ii<numBodies; ii++)
    {
        // Discard header lines and check for body type
        for(int ii=0; ii<3; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }

        // Get body type line
        fgets(buffer_line, sizeof(buffer_line), pFile);

        // Read body
        if (strncmp(buffer_line, "RAD_DIFF", 8) == 0)
        {
            pBodies[ii] = new Body(ii, this);
            pBodies[ii]->ReadPropertiesASCII(pFile);
            pBodies[ii]->OpenOutputFilesASCII(outputFolderPath);
        }
        else
        {
            std::stringstream ss;
            ss << "Error while parsing file: datosBodies.dat\n --> Expected body: " << ii <<" type definition\n";
            throw ValueError(ss.str());
        }
	}
    std::cout << "All bodies read" << std::endl;
    // Loop over bodies in order to get the number of hydrodynamic databases
    hydro_databases_name[hydro_database_count] = pBodies[0]->hydroDatabaseName;
    hydro_database_count++;
    for (int ii=1; ii<this->numBodies; ii++)
    {
        diff_count = 0;
        for (int jj=0; jj<hydro_database_count; jj++)
        {
            if (pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[jj]) != 0)
            {
                diff_count++;
            }
        }

        if (diff_count == hydro_database_count)
        {
            hydro_databases_name[hydro_database_count] = pBodies[ii]->hydroDatabaseName;
            hydro_database_count++;
        }
    }

    for (int ii=0; ii<hydro_database_count+1; ii++)
    {
        std::cout << hydro_databases_name[ii].c_str() << std::endl;
    }

    // Arrange all the bodies by database
    Body** pBodiesSort = new Body* [numBodies];
    int* pBody_found = new int [numBodies];
    for (int ii=0; ii<numBodies; ii++)
    {
        pBody_found[ii] = 0;
    }

    for (int ii=0; ii<hydro_database_count; ii++)
    {
        for (int jj=0; jj<numBodies; jj++)
        {
            if ((hydro_databases_name[ii].compare(pBodies[jj]->hydroDatabaseName) == 0) && (pBody_found[jj] == 0))
            {
                pBody_found[jj] = 1;
                pBodiesSort[body_count] = pBodies[jj];
                body_count++;
            }
        }
    }
    delete [] pBody_found;

    for (int ii=0; ii<numBodies; ii++)
    {
        pBodies[ii] = pBodiesSort[ii];
    }

    delete [] pBodiesSort;

    // Create an array in order to store the indexes of the bodies in each database
    int **check_hydro_bodies_id = new int* [hydro_database_count];
    Body*** check_hydro_bodies = new Body** [hydro_database_count];
    for (int ii=0; ii<hydro_database_count; ii++)
    {
        check_hydro_bodies_id[ii] = new int [max_num_bodies_database+1];
        check_hydro_bodies[ii] = new Body* [max_num_bodies_database];
    }
    for (int ii=0; ii<hydro_database_count; ii++)
    {
        for (int jj=0; jj<max_num_bodies_database+1; jj++)
        {
            check_hydro_bodies_id[ii][jj] = 0;
        }
    }

    // Check if there is some repeated body definition in each database
    std::cout << "Looking for body repetition..." << std::endl;
    for (int ii=0; ii<this->numBodies; ii++)
    {
        std::cout << "Looking for body repetition..." << std::endl;
        // Find database position inside the array of names generated previously
        pos_database = 0;
        while (true)
        {
            if (pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[pos_database]) == 0)
            {
                break;
            }
            pos_database++;
            if (pos_database >= hydro_database_count)
            {
                std::stringstream ss;
                ss << "It is not possible to find the name of the hydro database: " << pBodies[ii]->hydroDatabaseName;
                ss << " in the list of the hydrodatabase names done with bodies definition.";
                throw ValueError(ss.str());
            }
        }

        // Check if the body id already exist
        for (int jj=1; jj<=check_hydro_bodies_id[pos_database][0]; jj++)
        {
            if (check_hydro_bodies_id[pos_database][jj] == pBodies[ii]->hydroDatabaseIndex)
            {
                std::stringstream ss;
                ss << "Repeated Hydrodynamic Bodoy Index(" << pBodies[ii]->hydroDatabaseIndex <<") definition for Body: ";
                ss << ii << " and hydrodynamic database: " << hydro_databases_name[pos_database];
                throw ValueError(ss.str());
            }
        }
        std::cout << check_hydro_bodies_id[pos_database][0] << std::endl;
        check_hydro_bodies_id[pos_database][0]++;
        check_hydro_bodies_id[pos_database][check_hydro_bodies_id[pos_database][0]] = pBodies[ii]->hydroDatabaseIndex;
        check_hydro_bodies[pos_database][check_hydro_bodies_id[pos_database][0]-1] = pBodies[ii];
    }

    // Read hydrodynamic databases
    /**
    HydroDatabase** hydro_databases = new HydroDatabase* [hydro_database_count];
    for (int ii=0; ii<hydro_database_count; ii++)
    {
        
        
    }
    **/

    // Set hydrodynamic database to each body
    std::string hydro_file_path;
    for (int ii=0; ii<this->numBodies; ii++)
    {
        // Look for position of the database
        pos_database = 0;
        while (true)
        {
            if (this->pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[pos_database])==0)
            {
                break;
            }
            pos_database++;
        }

        // Set database to the target Body object
        hydro_file_path = JoinPath(this->inputFolderPath, hydro_databases_name[pos_database]);
        this->pBodies[ii]->LoadHydrodynamicDatabase(check_hydro_bodies[pos_database]);
    }

    // Fill System Matrix
    std::cout << "Fill system matrix...\n";
    arma::span a1;
    arma::span a2;
    double db_shift;
    double body_shift;
    this->pSystemMatrix = new arma::mat(6*this->numBodies, 6*this->numBodies, arma::fill::zeros);
    this->pSystemMatrixInv = new arma::mat(6*this->numBodies, 6*this->numBodies, arma::fill::zeros);
    for (int ii=0; ii<this->numBodies; ii++)
    {
        // Look for position of the database
        pos_database = 0;
        while (true)
        {
            if (this->pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[pos_database])==0)
            {
                break;
            }
            pos_database++;
        }

        db_shift = 0;
        for (int jj=0; jj<pos_database; jj++)
        {
            db_shift += 6*check_hydro_bodies_id[pos_database][0];
        }

        // Look for position of the body
        body_shift = 6*pBodies[ii]->hydroDatabaseIndex;

        // Fill system matrix
        a1 = arma::span(db_shift+body_shift,db_shift+body_shift+5);
        //std::cout << "a1: " << db_shift+body_shift << " - " << db_shift+body_shift+5 << "\n";
        a2 = arma::span(db_shift,db_shift+6*check_hydro_bodies_id[pos_database][0]-1);
        //std::cout << "a2: " << db_shift <<  " - " << db_shift+6*check_hydro_bodies_id[pos_database][0]-1 << "\n";
        (*pSystemMatrix)(a1, a2) += pBodies[ii]->pHydro->GetTotalMass();
        pBodies[ii]->sysMatSpan1 = a1;
        pBodies[ii]->sysMatSpan2 = a2;
    }
    std::cout << "Inverting system matrix...\n";
    *pSystemMatrixInv = arma::solve(*pSystemMatrix,eye(size(*pSystemMatrix)));
    //std::string filename = JoinPath(outputFolderPath, "SysyemMatrix.txt");  
    //pSystemMatrix->save(filename,arma::raw_ascii);
    std::cout << "System matrix inverted...\n";
    // Check the simulation time
    int time_buffer_size = this->timeBufferSize;
    for (int ii=0; ii<this->numBodies; ii++)
    {
        if (time_buffer_size < 10*this->pBodies[ii]->pHydro->GetNumPointsIrf())
        {
            time_buffer_size = 10*this->pBodies[ii]->pHydro->GetNumPointsIrf();
        }
    }

    for (int ii=0; ii<this->numBodies; ii++)
    {
        this->pBodies[ii]->velBufferSize = time_buffer_size;
        this->pBodies[ii]->velBuffer = arma::zeros(6, time_buffer_size);
    }
    this->timeBufferSize = time_buffer_size;
    this->timeBuffer = arma::zeros(1, time_buffer_size);

    // Free memory
    //delete[] hydro_databases;
    delete[] check_hydro_bodies_id;
    delete[] check_hydro_bodies;

    // Close file
    fclose(pFile);

    std::cout << "----> Bodies Properties Read" << std::endl;
    }
}


void Simulation::ReadBodiesHDF5()
{
    std::cout << "--> Reading Bodies Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBodiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Bodies Properties Read" << std::endl;
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
            pLines[ii]->OpenOutputFilesASCII(outputFolderPath);
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


void Simulation::ReadSinking()
{
    (this->*pReadSinking)();
}


void Simulation::ReadSinkingASCII()
{
    std::cout << "--> Reading Sinking Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine [1000];
    int sinkingBodyIndex;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataSinking.dat");
    FILE* pFile = fopen(file_path.c_str(), "r");
	
	if (pFile == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: dataSinking.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

	// Read total number of sinking bodies to read 
    fscanf(pFile, "%d %[^\n]\n", &numSinking, bufferLine);

    if (numSinking>0) {

        if (numSinking>1) {
            std::stringstream ss;
            ss << "Multiple bodies sinking is not implemented yet." << std::endl;
            throw IOError(ss.str());
        }

	    //Read all sinking bodies
	    pSinking = new Sinking* [numSinking];
	    for(int ii=0; ii<numSinking; ii++)
	    {
	    	// Discard header lines
	        for(int ii=0; ii<3; ii++)
	        {
	            fgets(bufferLine, sizeof(bufferLine), pFile);
	        }

	        // Get sinking body id
	        fscanf(pFile, "%d %[^\n]\n", &sinkingBodyIndex, bufferLine);

	        // Initiallice Sinking object
	    	pSinking[ii] = new Sinking(sinkingBodyIndex, this);

	    	// Read sinking body properties
	    	pSinking[ii]->ReadPropertiesASCII(pFile,inputFolderPath);

            pSinking[ii]->OpenOutputFilesASCII(outputFolderPath);
	    }

	}

    // Close file
    fclose(pFile);

    std::cout << "----> Sinking Properties Read" << std::endl;
}


void Simulation::ReadSinkingHDF5()
{
    std::cout << "--> Reading Sinking Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadSinkingHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Sinking Properties Read" << std::endl;
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
	fscanf(file_pointer, "%lf %[^\n]\n", &hydroTimeStep, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &fastTimeStep, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &fastControllerTimeStep, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIRF, bufferLine);
	fscanf(file_pointer, "%lf %[^\n]\n", &sinkingTimeStep, bufferLine);
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


void Simulation::ReadWaves()
{
    (this->*pReadWaves)();
}


void Simulation::ReadWavesASCII()
{
	std::cout << "--> Reading Waves (ASCII format)" << std::endl;
    std::string file_path = JoinPath(inputFolderPath, "dataWaves.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: dataWaves.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
	
    char bufferLine [1000];

    //Ignoro las tres primeras lineas
	for(int ii=0; ii<3; ii++)
	{
		fgets(bufferLine, sizeof(bufferLine), file_pointer);
	}

    // Get wave type line
    char wave_type [1000];
    fgets(wave_type, sizeof(wave_type), file_pointer);
    double H; double T; double D;
    fscanf(file_pointer, "%lf %[^\n]\n", &H, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &T, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &D, bufferLine);

    // Read body
    if (strncmp(wave_type, "REG", 3) == 0)
    {
        pWave = new RegularWave(H,T,D);
    }
    else
    {
        if (strncmp(wave_type, "IRR", 3) == 0)
		{
		    pWave = new IrregularWave(H,T,D);
		    for(int ii=0; ii<3; ii++)
			{
				fgets(bufferLine, sizeof(bufferLine), file_pointer);
			}
			fscanf(file_pointer, "%d %[^\n]\n", &pWave->specType_flag, bufferLine);
			fgets(bufferLine, sizeof(bufferLine), file_pointer);
			fscanf(file_pointer, "%lf %[^\n]\n", &pWave->gamma, bufferLine);
			fscanf(file_pointer, "%lf %[^\n]\n", &pWave->s, bufferLine);
			fscanf(file_pointer, "%lf %[^\n]\n", &pWave->dtheta, bufferLine);
			fscanf(file_pointer, "%lf %[^\n]\n", &pWave->rel_tol, bufferLine);
			fscanf(file_pointer, "%lf %[^\n]\n", &pWave->dt, bufferLine);
			fscanf(file_pointer, "%lf %[^\n]\n", &pWave->factor, bufferLine);
			fgets(bufferLine, sizeof(bufferLine), file_pointer);
			char cWaveDatabaseName [1000];
			fscanf(file_pointer, "%s %[^\n]\n", cWaveDatabaseName, bufferLine);
			pWave->waveDatabaseName = cWaveDatabaseName;
			pWave->file_path = JoinPath(inputFolderPath, pWave->waveDatabaseName);
		}
		else
		{
		    std::stringstream ss;
		    ss << "Error while parsing file: dataWaves.dat; Unexpected wave type. \n";
		    throw ValueError(ss.str());
		}
    }

    // Close file
    fclose(file_pointer);

    // Show inputs
    if (true)
    {
    	if (strncmp(wave_type, "REG", 3) == 0)
    	{
    		std::cout << "Wave type: Regular" << std::endl;
    		std::cout << "Wave height: " << H << std::endl;
    		std::cout << "Wave period: " << T << std::endl;
    		std::cout << "Wave heading: " << D << std::endl;
    	}
    	else
    	{
    		std::cout << "Wave type: Irregular" << std::endl;
    		std::cout << "Wave significant height: " << H << std::endl;
    		std::cout << "Wave peak period: " << T << std::endl;
    		std::cout << "Wave heading: " << D << std::endl;
    		if (pWave->specType_flag==1)
    		{
    			std::cout << "Wave peak enhacement factor: " << pWave->gamma << std::endl;
    			std::cout << "Wave directional spreading: " << pWave->s << std::endl;
    			std::cout << "Wave directional step: " << pWave->dtheta << std::endl;
    			std::cout << "Relative tolerance for wave check: " << pWave->rel_tol << std::endl;
    			std::cout << "Time step for wave check: " << pWave->dt << std::endl;
    		}
    		else
    		{
    			std::cout << "Wave base data file name: " << pWave->waveDatabaseName << std::endl;
    		}
    	}
    }

    std::cout << "----> Waves Read" << std::endl;

    pWave->simulationTime = simulationTime;
    pWave->gravity = gravity;
    pWave->waterDepth = abs(waterDepth);
    pWave->CheckBreakingWave();
    pWave->GetWaveSpectrum();
    if (strncmp(wave_type, "IRR", 3) == 0)
    {
        pWave->WriteOut(outputFolderPath);
    }

    //std::string filename;
    //filename = JoinPath(outputFolderPath, "freqs.txt");  
    //pWave->freqs.save(filename,arma::raw_ascii);
    //filename = JoinPath(outputFolderPath, "headings.txt");  
    //pWave->headings.save(filename,arma::raw_ascii);
    //filename = JoinPath(outputFolderPath, "k.txt");  
    //pWave->k.save(filename,arma::raw_ascii);
    //filename = JoinPath(outputFolderPath, "kx.txt");  
    //pWave->kx.save(filename,arma::raw_ascii);
    //filename = JoinPath(outputFolderPath, "ky.txt");  
    //pWave->ky.save(filename,arma::raw_ascii);
}


void Simulation::ReadWavesHDF5()
{
	std::cout << "--> Reading Waves (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWavesHDF5 in class Simulation not implemented yet.";
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
    printf("NumWinches: %d - UseWinches: %d\n", numWinches, useWinches);
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

    std::cout << "--> Reading Winches Controller Properties (ASCII format)" << std::endl;
    file_path = JoinPath(inputFolderPath, "datosWinchiesController.dat");
    file_pointer = fopen(file_path.c_str(), "r");
    WinchesController = WinchieController(numWinches,pWinches,this);
    WinchesController.ReadPropertiesASCII(file_pointer);
    fclose(file_pointer);
    WinchesController.OpenOutputFilesASCII(outputFolderPath);
    std::cout << "----> Winches Controller Properties Read" << std::endl;
}


void Simulation::ReadWinchesHDF5()
{
    std::cout << "--> Reading Winches Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWinchesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Winches Properties Read" << std::endl;
}


void Simulation::ReadWindTurbines(void)
{
     (this->*pReadWindTurbines)();
}


void Simulation::ReadWindTurbinesASCII(void)
{
    std::cout << "--> Reading Wind Turbines Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine [1000];

    // Open file
	std::string file_path = JoinPath(inputFolderPath, "datosWindTurbines.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: datosWindTurbines.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
	}

    // Read number of Wind Turbines defined in the file
	fscanf(file_pointer, "%d %[^\n]\n", &numWindTurbines, bufferLine);
    // Allocate a vector of pointers to WindTurbine class objects
    pWindTurbines = new WindTurbine* [numWindTurbines];
    for(int ii=0; ii<numWindTurbines; ii++)
    {
		pWindTurbines[ii] = new WindTurbine(ii,this);
		pWindTurbines[ii]->ReadPropertiesASCII(file_pointer);
	}

    // Close the file
    fclose(file_pointer);
    std::cout << "----> Wind Turbines Properties Read" << std::endl;

}


void Simulation::ReadWindTurbinesHDF5(void)
{
    std::cout << "--> Reading Wind Turbines (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWindTurbinesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Wind Turbines Properties Read" << std::endl;
}


void Simulation::Run()
{
    time_t tstart, tend;
	double wallTime = 0.0;
	double wallTimeHydro = 0.0;
    double wallTimeFAST = 0.0;
    double wallTimeControllerFAST = 0.0;
	double wallTimeSinking = 0.0;
    tstart = time(0);
    std::cout<< "    t = " << wallTime << " s" << std::endl;
    //pTimeSolver->dt_max = 0.1;
    // std::cout<< "In Simulation::Run --> Starting temporal integration loop "<< std::endl;
    do
    {
        pTimeSolver->step();
        // std::cout<< "In Simulation::Run --> Call to step() was succesfull "<< std::endl;
        // Print out time if any
        //std::cout<< "    t = " << pTimeSolver->t << " s"  << std::endl;
        /**
        for(int ii=0; ii<numLines; ii=ii+1) pLines[ii]->WriteOut(pTimeSolver->t);
        for(int ii=0; ii<numBodies; ii=ii+1) pBodies[ii]->WriteOut(pTimeSolver->t);
        **/
        if (pTimeSolver->t >= wallTime + maxTimeStep)
        {
            wallTime = wallTime + maxTimeStep;
            std::cout<< "    t = " << wallTime << " s"  << std::endl;
            if (numWinches>0) {
            	WinchesController.controlWinchies(wallTime);
            	WinchesController.WriteOut(wallTime);
            }
            for(int ii=0; ii<numLines; ii=ii+1) pLines[ii]->WriteOut(wallTime);
            for(int ii=0; ii<numBodies; ii=ii+1) pBodies[ii]->WriteOut(wallTime);
        }

    	if (pTimeSolver->t >= wallTimeHydro + hydroTimeStep)
    	{
            wallTimeHydro += hydroTimeStep;
    		if (numBodies>0){   
		        UpdateSystem();
		    }
            for(int ii=0; ii<numBodies; ii=ii+1) 
            {
            	pBodies[ii]->Fb = pBodies[ii]->pHydro->CalculateHydrodynamicForces(wallTimeHydro);
            }
    	}

        if (numWindTurbines>0) {
            if (pTimeSolver->t >= wallTimeFAST + fastTimeStep)
            {
                wallTimeFAST += fastTimeStep;
                for(int ii=0; ii<numWindTurbines; ii=ii+1){
                    pWindTurbines[ii]->SetInputsFAST();
                    pWindTurbines[ii]->ComputeForces(wallTimeFAST);
                }
            }
            if (pTimeSolver->t >= wallTimeControllerFAST + fastControllerTimeStep)
            {
                wallTimeControllerFAST += fastControllerTimeStep;
                for(int ii=0; ii<numWindTurbines; ii=ii+1){
                    pWindTurbines[ii]->SetInputsFAST();
                    pWindTurbines[ii]->ComputeControler(wallTimeControllerFAST);
                }
            }
        }

    	
    	if (numSinking>0) {   
	        if (pTimeSolver->t >= wallTimeSinking + sinkingTimeStep) {
	        	wallTimeSinking += sinkingTimeStep;
	        	for(int ii=0; ii<numSinking; ii=ii+1) {
	            	pSinking[ii]->UpdateSinkingHydrodynamics(wallTime);
	            }
	            UpdateSystemMatrix();
                for(int ii=0; ii<numSinking; ii=ii+1) pSinking[ii]->WriteOut(wallTime);
	        }
	    }
	    
        
    } while (pTimeSolver->t <= simulationTime);
    
    tend = time(0); 
    std::cout << std::endl << "    Computational time  : " << difftime(tend, tstart) << " seconds" << std::endl;
    std::cout << "    Total function calls: " << numCallsSysFun << std::endl;
    std::cout << "    Total jac calls: " << pTimeSolver->iJ << std::endl << std::endl;

    if (writeEquilibrium == 1) {
        std::cout << "  Writting data to Equilibrio.dat ..." << std::endl << std::endl; //////////////////////////////////////////
        std::string filename = JoinPath(outputFolderPath, "Equilibrio.dat");
        pTimeSolver->y.save(filename,arma::arma_ascii);
    }
}


void Simulation::SetupCase()
{
    std::cout << "----> Setting up the case configuration ..." << std::endl;

    // Checking free bodies
    std::cout << "        Checking for free bodies ..." << std::endl;
    for (int ii=0; ii<numBodies; ii++)
    {
        if (pBodies[ii]->flag_blocked>0){
            numBodiesLock++;
        } else {
            numBodiesFree++;
        }
    }
    pBodiesLock = new Body* [numBodiesLock];
    pBodiesFree = new Body* [numBodiesFree];
    int indL = 0; int indF = 0;
    for (int ii=0; ii<numBodies; ii++)
    {
        if (pBodies[ii]->flag_blocked>0){
            pBodiesLock[indL] = pBodies[ii];
            indL++;
        } else {
            pBodiesFree[indF] = pBodies[ii];
            indF++;
        }
    }

    // Count the number of BCP in each body and create pointer array
    std::cout << "        Counting the number of BCPs in each body ..." << std::endl;
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
    bool* pDefined_body_bcps = new bool[numBcps];
    for (int ii=0; ii<numBcps; ii++)
    {
        pDefined_body_bcps[ii] = 0;
    }
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<pBodies[ii]->numBcps; jj++)
        {
            if (!pDefined_body_bcps[pBodies[ii]->pIndexBcps[jj]])
            {
                pBcps[pBodies[ii]->pIndexBcps[jj]]->pBodies = new Body* [pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp];
                pDefined_body_bcps[pBodies[ii]->pIndexBcps[jj]] = true;
            }
        }
    }
    delete [] pDefined_body_bcps;    
    std::cout << "        ... done!" << std::endl;

    // Assing to each BCP the corresponding Body pointer
    std::cout << "        Assigning to each BCP the corresponding body  ..." << std::endl;
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
    std::cout << "        ... done!" << std::endl;

    // Count the number of Wind Turbines in each body and create pointer array
    std::cout << "        Checking the number of Wind Turbines in each body ..." << std::endl;
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<pBodies[ii]->numWindTurbs; jj++)
        {
            if (pBodies[ii]->pIndexWindTurbs[jj]+1 > numWindTurbines)
            {
                std::stringstream ss;
                ss << "BCP index: " << pBodies[ii]->pIndexWindTurbs[jj] << " in Body: " << pBodies[ii]->GetId() \
                    << " is out of range when compare with the Number of BCPs(" << numBcps << ") defined in" \
                    << " datosBCPs.dat";
                throw ValueError(ss.str());
            }
        }
    }
    for (int ii=0; ii<numBodies; ii++)
    {
        for (int jj=0; jj<numBodies; jj++)
        {
            if (ii!=jj && pBodies[ii]->numWindTurbs>0 && pBodies[jj]->numWindTurbs>0){
                for (int kii=0; kii<pBodies[ii]->numWindTurbs; kii++)
                {
                    int indWTloc = pBodies[ii]->pIndexWindTurbs[kii];
                    for (int kjj=0; kjj<numBodies; kjj++)
                    {
                        if(indWTloc==pBodies[jj]->pIndexWindTurbs[kjj])
                        {
                            std::stringstream ss;
                            ss << "Bodies " << ii+1 << " and " << jj+1 << " share wind turbine" << indWTloc << "!";
                            throw ValueError(ss.str());
                        }
                    }
                }
            }
        }
    }
    std::cout << "        ... done!" << std::endl;

    
    // Assing to each Body the corresponding Wind Turbine pointers
    std::cout << "       Assingning to each Body the corresponding Wind Turbine pointers  ..." << std::endl;
    for (int ii=0; ii<numBodies; ii++)
    {
        for(int jj=0; jj<pBodies[ii]->numWindTurbs; jj++)
        {
            pBodies[ii]->pBodyWindTurbs[jj] = pWindTurbines[pBodies[ii]->pIndexWindTurbs[jj]];
        }
    }    
    std::cout << "        ... done!" << std::endl;

    // Count the number of lines in each joint BCP
    std::cout << "        Counting the number of lines in each joint BCP ..." << std::endl;
    for (int ii=0; ii<numBcps; ii++)
    {
        if (pBcps[ii]->GetType()==3)
        {
            int temp_nL = 0;
            int temp_BCP_Id = pBcps[ii]->GetId();
            for (int jj=0; jj<numLines; jj++)
            {
                if (pLines[jj]->indexBcps[0]==temp_BCP_Id)
                {
                    temp_nL++;
                }
                if (pLines[jj]->indexBcps[1]==temp_BCP_Id)
                {
                    temp_nL++;
                }
            }
            pBcps[ii]->posLines = arma::zeros(temp_nL,3);
            pBcps[ii]->velLines = arma::zeros(temp_nL,3);
        }
    }
    std::cout << "        ... done!" << std::endl;

    // Count the number of Lines in each body and create pointer array
    std::cout << "        Counting the number of lines in each body ..." << std::endl;
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
    bool* pDefined_lines_bcps = new bool[numBcps]; 
    for (int ii=0; ii<numBcps; ii++)
    {
        pDefined_lines_bcps[ii] = 0;
    }
    for (int ii=0; ii<numLines; ii++)
    {
        for (int jj=0; jj<pLines[ii]->numBcps; jj++)
        {
            if (!pDefined_lines_bcps[pLines[ii]->indexBcps[jj]])
            {
                pBcps[pLines[ii]->indexBcps[jj]]->pLines = new Line* [pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pDefined_lines_bcps[pLines[ii]->indexBcps[jj]] = true;
            }
        }
        
    }
    delete [] pDefined_lines_bcps;    
    std::cout << "        ... done!" << std::endl;

    // Assing to each Line the corresponding BCP pointer
    std::cout << "        Assigning to each Line the corresponding BCP ..." << std::endl;
    for (int ii=0; ii<numLines; ii++)
    {
        for(int jj=0; jj<pLines[ii]->numBcps; jj++)
        {
            pLines[ii]->pLineBcps[jj] = pBcps[pLines[ii]->indexBcps[jj]];
            pLines[ii]->pLineBcps[jj]->pLines[pLines[ii]->pLineBcps[jj]->countLine] = pLines[ii];
            pLines[ii]->pLineBcps[jj]->countLine++;
        }

        if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() != 3){
        	numDofTotal += (pLines[ii]->N-2);
        	pLines[ii]->first_node = 1;
        	pLines[ii]->last_node = pLines[ii]->N-1;
        } else if (pLines[ii]->pLineBcps[0]->GetType() == 3 && pLines[ii]->pLineBcps[1]->GetType() == 3){
        	numDofTotal += pLines[ii]->N;
        	pLines[ii]->first_node = 0;
        	pLines[ii]->last_node = pLines[ii]->N;
        } else if (pLines[ii]->pLineBcps[0]->GetType() == 3 && pLines[ii]->pLineBcps[1]->GetType() != 3){
        	numDofTotal += (pLines[ii]->N - 1);
        	pLines[ii]->first_node = 0;
        	pLines[ii]->last_node = pLines[ii]->N-1;
        } else if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() == 3){
        	numDofTotal += (pLines[ii]->N-1);
        	pLines[ii]->first_node = 1;
        	pLines[ii]->last_node = pLines[ii]->N;
        }

        try
        {
            if(!readEquilibrium) pLines[ii]->initLine();
            pLines[ii]->print_out();
            pLines[ii]->SEM_getBaseFunctions();
            std::string filename = JoinPath(outputFolderPath,"DerivativeMatrix_");
            filename = filename + "_Line_" + std::to_string(ii+1) + ".dat";  
            pLines[ii]->D.save(filename,arma::arma_ascii);
            std::string filename2 = JoinPath(outputFolderPath,"ArcLengthPoints_");
            filename2 = filename2 + "_Line_" + std::to_string(ii+1) + ".dat";  
            pLines[ii]->s.save(filename2,arma::arma_ascii);
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
    std::cout << "        ... done!" << std::endl;

    std::cout << "        Counting the number of line nodes without repetition of joint nodes ..." << std::endl;
    // Count the number of line nodes without repetition of joint nodes
    int numUsedJointBCPs = 0;
    for (int jj=0; jj<numLines; jj++){
        numAllLinesNodes += pLines[jj]->N;
        for (int kk=0; kk<2; kk++){
            if((pLines[jj]->pLineBcps[kk]->GetType() == 3)&&(pLines[jj]->pLineBcps[kk]->flag_counted == 0)){
                    pLines[jj]->pLineBcps[kk]->flag_counted = 1;
                    numUsedJointBCPs++;
            }
        }
    }
    numAllLinesNodes -= numUsedJointBCPs;

    // Build the lines coupling sparse matrix and store the lines index vectors
    int indFirstNodeAvail = 0;    
    for (int jj=0; jj<numLines; jj++){

        int numLineNodes_tmp = pLines[jj]->N;

        pLines[jj]->ind4CouplingMat = arma::zeros<arma::uvec>(numLineNodes_tmp);

        // First node
        if (pLines[jj]->pLineBcps[0]->flag_assigned == 0){
            pLines[jj]->pLineBcps[0]->flag_assigned = 1;
            pLines[jj]->pLineBcps[0]->couplingMatIndex = indFirstNodeAvail;
            indFirstNodeAvail++;
        }
        pLines[jj]->ind4CouplingMat(0) = pLines[jj]->pLineBcps[0]->couplingMatIndex;

        // Intermediate nodes
        for (int kk=1; kk<numLineNodes_tmp-1; kk++){
            pLines[jj]->ind4CouplingMat(kk) = indFirstNodeAvail;
            indFirstNodeAvail++;
        }
        
        // Last node
        if (pLines[jj]->pLineBcps[1]->flag_assigned == 0){
            pLines[jj]->pLineBcps[1]->flag_assigned = 1;
            pLines[jj]->pLineBcps[1]->couplingMatIndex = indFirstNodeAvail;
            indFirstNodeAvail++;
        }
        pLines[jj]->ind4CouplingMat(numLineNodes_tmp-1) = pLines[jj]->pLineBcps[1]->couplingMatIndex;

        std::cout << "        Line " << jj+1 << " indices: " << std::endl << pLines[jj]->ind4CouplingMat << std::endl;

    }

    std::cout << "        Computing Lines Coupling Matrix ..." << std::endl;
    pLinesCouplingMatrix = new arma::mat(numAllLinesNodes,numAllLinesNodes,arma::fill::zeros);
    pLinesCouplingMatrixInv = new arma::mat(numAllLinesNodes,numAllLinesNodes);
    pLinesCouplingMatrix_sp = new arma::sp_mat(numAllLinesNodes, numAllLinesNodes);
    ComputeLinesCouplingMatrix();
    if(numAllLinesNodes<100){
        *pLinesCouplingMatrixInv = arma::solve(*pLinesCouplingMatrix,eye(size(*pLinesCouplingMatrix)));
        std::string filename = JoinPath(outputFolderPath,"LinesCouplingMatrix.dat");
        (*pLinesCouplingMatrix).save(filename,arma::arma_ascii);
    }    
    std::cout << "        ... done!" << std::endl;

    // Setup Springs
    std::cout << "        Setting up springs ..." << std::endl;
	for(int ii=0; ii<numSprings; ii++)
    {
        std::cout << "            Spring: " << ii << "\n";
        std::cout << "            Spring:->BCP_1 " << pSprings[ii]->BCP_1 << "\n";
        std::cout << "            Spring:->BCP_2 " << pSprings[ii]->BCP_2 << "\n";
		pSprings[ii]->SpringBCP[0] = pBcps[pSprings[ii]->BCP_1];
		pSprings[ii]->SpringBCP[1] = pBcps[pSprings[ii]->BCP_2];
	}    
    std::cout << "        ... done!" << std::endl;

	// Setup hidro data bases
	std::cout << "        Setting up hydro data bases ..." << std::endl;
	for (int ii=0; ii<numBodies; ii++)
    {
        std::cout << "            Body: " << ii << "\n";
    	pBodies[ii]->pHydro->SetUp();
    }    
    std::cout << "        ... done!" << std::endl;

    // Setup winchies controller
    if (useWinches) {
		std::cout << "        Setting up winchies controller ..." << std::endl;
    	WinchesController.SetUpWinchiesController();
    std::cout << "        ... done!" << std::endl;
	}

    if (numWindTurbines>0)
    {
        // Setup wind turbines
        std::cout << "        Setting up wind turbines ..." << std::endl;
        for (int ii=0; ii<numWindTurbines; ii++)
        {
            std::cout << "            Turbine: " << ii << "\n";
            pWindTurbines[ii]->Initialize();
        }    
        std::cout << "        ... done!" << std::endl;


        // Computing bodies structural mass considering wind turbines
        std::cout << "        Computing bodies structural mass considering wind turbines ..." << std::endl;
        for (int ii=0; ii<numBodies; ii++)
        {
            if (pBodies[ii]->numWindTurbs>0) 
            {
                std::cout << "            Body: " << ii << "\n";
                // Aqui seria conveniente comprobar que las distintas turbinas 
                // la inercia de la HDB son al menos muy similares.
                pBodies[ii]->inertia = pBodies[ii]->pBodyWindTurbs[0]->bodyInerMat;
                for (int jj=0; jj<pBodies[ii]->numWindTurbs; jj++)
                {
                    
                    pBodies[ii]->inertia += pBodies[ii]->pBodyWindTurbs[0]->towrInerMat;
                    pBodies[ii]->inertia += pBodies[ii]->pBodyWindTurbs[0]->turbInerMat;
                }
            }
        }    
        std::cout << "        ... done!" << std::endl;
    }

    std::cout << "----> Case configuration done" << std::endl;
}


void Simulation::ComputeLinesCouplingMatrix(void){

    for (int jj=0; jj<numLines; jj++){

        int numLineNodes_tmp = pLines[jj]->N;

        arma::mat LineMassMat_tmp = pLines[jj]->MM * pLines[jj]->dL;
        if(pLines[jj]->pLineBcps[0]->GetType() != 3){
            LineMassMat_tmp.row(0) = arma::zeros(1,numLineNodes_tmp); 
            LineMassMat_tmp(0,0) = 1.0;
        }
        if(pLines[jj]->pLineBcps[1]->GetType() != 3){
            LineMassMat_tmp.row(numLineNodes_tmp-1) = arma::zeros(1,numLineNodes_tmp); 
            LineMassMat_tmp(numLineNodes_tmp-1,numLineNodes_tmp-1) = 1.0;
        }

        if(numAllLinesNodes>=100){
            for (int irow=0; irow<numLineNodes_tmp; irow++){
                for (int icol=0; icol<numLineNodes_tmp; icol++){
                    (*pLinesCouplingMatrix_sp)(pLines[jj]->ind4CouplingMat(irow),pLines[jj]->ind4CouplingMat(icol)) += LineMassMat_tmp(irow,icol);
                }
            }
        } else {
            (*pLinesCouplingMatrix).submat(pLines[jj]->ind4CouplingMat,pLines[jj]->ind4CouplingMat) += LineMassMat_tmp;
        }

    }


    for(int jj=0; jj<numBcps; jj++){
        if((pBcps[jj]->GetType() == 3) && (pBcps[jj]->flag_assigned == 1)){
            if(numAllLinesNodes>=100){
                (*pLinesCouplingMatrix_sp)(pBcps[jj]->couplingMatIndex,pBcps[jj]->couplingMatIndex) += pBcps[jj]->mass_Joint;
            } else {
                (*pLinesCouplingMatrix)(pBcps[jj]->couplingMatIndex,pBcps[jj]->couplingMatIndex) += pBcps[jj]->mass_Joint;
            }
        }
    }

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
        pReadProperties = &Simulation::ReadPropertiesASCII;
        pReadWaves = &Simulation::ReadWavesASCII;
        pReadBcps = &Simulation::ReadBcpsASCII;
        pReadBodies = &Simulation::ReadBodiesASCII;
        pReadSinking = &Simulation::ReadSinkingASCII;
        pReadLines = &Simulation::ReadLinesASCII;
        pReadSprings = &Simulation::ReadSpringsASCII;
        pReadWinches = &Simulation::ReadWinchesASCII;
        pReadWindTurbines = &Simulation::ReadWindTurbinesASCII;
    }
    else if (!incDataFormat.compare("HDF5"))
    {
        dataFormat = 1;
        dataFormatStr = "HDF5";
        inputFolderPath = incProjectPath;
        outputFolderPath = incProjectPath;
        pReadProperties = &Simulation::ReadPropertiesHDF5;
        pReadWaves = &Simulation::ReadWavesHDF5;
        pReadBcps = &Simulation::ReadBcpsHDF5;
        pReadBodies = &Simulation::ReadBodiesHDF5;
        pReadSinking = &Simulation::ReadSinkingHDF5;
        pReadLines = &Simulation::ReadLinesHDF5;
        pReadSprings = &Simulation::ReadSpringsHDF5;
        pReadWinches = &Simulation::ReadWinchesHDF5;
        pReadWindTurbines = &Simulation::ReadWindTurbinesHDF5;
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
    if (numBodies>0){
        // Update time vector if any
        bool restoreMatrix = false;
        timeBufferCount++;
        if (timeBufferCount < timeBufferSize)
        {
            timeBuffer(0, timeBufferCount) = pTimeSolver->t;
        }
        else
        {
            arma::mat timeBufferNew = arma::zeros(1, timeBufferSize);
            timeBufferNew.cols(0, pBodies[0]->pHydro->GetNumPointsIrf()-1) = timeBuffer.cols(timeBufferSize-pBodies[0]->pHydro->GetNumPointsIrf(), timeBufferSize-1);
            timeBuffer = timeBufferNew;
            timeBufferCount = pBodies[0]->pHydro->GetNumPointsIrf()-1;
            restoreMatrix = true;
        }

        // Update bodies velocity
        int ini = 0;
        for(int ii=0; ii<numBodies; ii++)
        {
            pBodies[ii]->StoreVelocities(restoreMatrix);
            ini = ini + 6;
        }
    }
}


void Simulation::UpdateSystemMatrix()
{
    *pSystemMatrix = arma::zeros(6*numBodies,6*numBodies);
    for (int ii=0; ii<this->numBodies; ii++)
    {
        (*pSystemMatrix)(pBodies[ii]->sysMatSpan1, pBodies[ii]->sysMatSpan2) += pBodies[ii]->pHydro->GetTotalMass();
    }
    *pSystemMatrixInv = arma::solve(*pSystemMatrix,eye(size(*pSystemMatrix)));
}
