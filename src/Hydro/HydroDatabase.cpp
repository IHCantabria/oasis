
#include <armadillo>
#include <cstdio>
#include <string>
#include <sstream>
#include <chrono>
#include "HydroDatabase.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"
#include "../Simulations/Simulation.hpp"
#include "../Bodies/Bodies.hpp"
#include "../MathTools.hpp"


arma::mat HydroDatabase::ComputeHydrostaticForces()
{	
	arma::mat hydrostatic_force = -(*pHydrostaticStiffness)*((*pBodies)[id].pos);
	return hydrostatic_force;
}


arma::mat HydroDatabase::ComputeRadiationForces(double t, solver_data SD)
{	
	// Allocate radiation force solution vector
	arma::mat radiation_force = arma::zeros(6, 1);

	// Declare local auxiliary variables
	arma::mat time_local;
	arma::mat irf_local;
	arma::mat irf_local_interp;
	arma::mat vel_local;
	arma::mat vel_local_interp;
	arma::mat vel_local_filter;
	int idx_begin, idx_end;
	double dt = IRFTime(0, 1) - IRFTime(0, 0);
	int aux;

	// Compute radiation forces for the 6DOFs
	for(int ib=0; ib<numBodies; ib++)
	{
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++)
			{
				if ((*SD.timeBuffer)(0, *SD.timeBufferCount) > IRFTime(0, (*pIRFPoints[ib])(i, j)))
				{
					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);
					
					// Calculate begin and end indexes for matrix slicing
					idx_begin = (*pBodies[ib]).velBufferCount-(*pIRFPoints[ib])(i, j)+1;
					idx_end = (*pBodies[ib]).velBufferCount-1;

					// Get velocity chunck from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, idx_begin, j, idx_end);

					// Interpolate local velocity vector in order to fit with IRF resolution
					vel_local_interp = interp1((*SD.timeBuffer).cols(idx_begin, idx_end)-(*SD.timeBuffer)(0, idx_end), vel_local, IRFTime.cols(0, (*pIRF[ib])(0, i, j)));

					// Calulate Duhamel integral
					vel_local_filter = arma::fliplr(irf_local)%(vel_local_interp.t());
					radiation_force(i) += trapz(vel_local_filter, dt);
				}
				else if ((*pBodies[ib]).velBufferCount > 0)
				{
					// Calculate begin and end indexes for matrix slicing
					idx_end = (*pBodies[ib]).velBufferCount-1;

					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);

					// Get velocity chunk from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, 0, j, idx_end);

					// Interpolate local velocity vector in order to fit with IRF resolution
					time_local = IRFTime.cols(0, numPointsIRF-1) - IRFTime(0, numPointsIRF-1) + (*SD.timeBuffer)(0, idx_end);
					irf_local_interp = interp1(time_local, arma::fliplr(irf_local), (*SD.timeBuffer).cols(0, idx_end));

					// Calulate Duhamel integral
					vel_local_filter = irf_local_interp%vel_local;
					radiation_force(i) += trapzi((*SD.timeBuffer).cols(0, idx_end), vel_local_filter);
				}
			}
		}
	}
	
	/**
	if((*pBodies)[id].velBufferCount > 0)
	{
		for(int ii=0; ii<6; ii++)
		{
			for(int jj=0; jj<6; jj++)
			{
				if(mySim.time > )
			}
		}
	}
	**/
	return radiation_force;
}


void HydroDatabase::ComputeIRF(void)
{
	// Declare local variables
	int count_max;
	int count_zero;
	int max_consec;
	arma::mat dummy_mat;
	arma::mat max_position;
	arma::mat zero_cross;
	arma::mat dampingFreq;
	
	// Check input arguments
	if ((*pFrequencies)(1)<(*pFrequencies)(0))
	{
		perror("Frequencies does not increase monotonically...\n");
	}
	
	// Calculate maximum time allowed
	double dt=0.05;
	double df = (*pFrequencies)(1)-(*pFrequencies)(0);
	double tmax = 1/df/2.0;
	IRFTime = arange(0, tmax, dt);
	
	// Allocate IRF matrix
	numPointsIRF = IRFTime.n_cols;
	pIRF = new arma::cube* [numBodies];
	pIRFPoints = new arma::mat* [numBodies];
	

	// Loop to find the IRF value for each body influence and DOF
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	for (int ib=0; ib<numBodies; ib++)
	{
		pIRFPoints[ib] = new arma::mat(6, 6, arma::fill::zeros);
		pIRF[ib] = new arma::cube(IRFTime.n_cols, 6, 6, arma::fill::zeros);
		for (int i=0; i<6; i++)
		{
			for (int j=0; j<6; j++)
			{
				// Clear previous results
				count_max = 1;
				count_zero = -1;
				max_consec = 0;
				max_position = arma::zeros(1, IRFTime.n_cols);
				zero_cross = arma::zeros(1, IRFTime.n_cols);
				
				// Start new Dof data
				
				dampingFreq = (*pDampingRadiation[ib]).subcube(i,j,0,i,j,numFrequencies-1);
				dummy_mat = dampingFreq%cos(2*M_PI*(*pFrequencies)*IRFTime(0, 0));
				(*pIRFPoints[ib])(i, j) = IRFTime.n_cols - 1;
				(*pIRF[ib])(0, i, j) = 2*trapz(dummy_mat, df);
				for (int k=1; k<IRFTime.n_cols; k++)
				{
					// Calculate new value of IRF
					dummy_mat = dampingFreq%cos(2*M_PI*(*pFrequencies)*IRFTime(0, k));
					(*pIRF[ib])(k, i, j) = 2*trapz(dummy_mat, df);
					
				}
			}
		}
	}
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time elapsed ComputeIRF: " << elapsed << std::endl;
}


int HydroDatabase::GetId(void)
{
	return id;
}


HydroDatabase::HydroDatabase(int incId, Body** incBody)
{
	id = incId;
	pBodies = incBody;
}


void HydroDatabase::Print()
{
	std::cout << "Number of bodies associated: " << numBodies << std::endl;
	std::cout << "Number of frequencies: " << numFrequencies << std::endl;
	std::cout << "Number of headings: " << numHeadings << std::endl;
}


void HydroDatabase::ReadHydroMechanicsHDF5(std::string filePath)
{
	// Read number of bodies
	arma::mat num_bodies_mat;
	num_bodies_mat.load(arma::hdf5_name(filePath, "num_bodies"));
	numBodies = num_bodies_mat(0);
	
	// Read frequencies
	std::stringstream frequencies_fn;
	pFrequencies = new arma::mat;
	frequencies_fn << "body_" << this->GetId() << "/frequencies";
	pFrequencies->load(arma::hdf5_name(filePath, frequencies_fn.str()));
	numFrequencies = pFrequencies->n_cols;
	
	// Read headings (radians)
	std::stringstream headings_fn;
	pHeadings = new arma::mat;
	headings_fn << "body_" << this->GetId() << "/headings";
	pHeadings->load(arma::hdf5_name(filePath, headings_fn.str()));
	numHeadings = pHeadings->n_cols;
	
	// Read hydrostatic stiffness
	std::stringstream hydrostatic_stiffness_fn;
	pHydrostaticStiffness = new arma::mat;
	hydrostatic_stiffness_fn << "body_" << this->GetId() << "/hydstiffness";
	pHydrostaticStiffness->load(arma::hdf5_name(filePath, hydrostatic_stiffness_fn.str(), arma::hdf5_opts::trans));
	
	// Read structural mass properties
	std::stringstream structural_mass_fn;
	pStructuralMass = new arma::mat;
	structural_mass_fn << "body_" << this->GetId() << "/mass";
	pStructuralMass->load(arma::hdf5_name(filePath, structural_mass_fn.str(), arma::hdf5_opts::trans));
	
	// Read Added Mass
	std::stringstream added_mass_fn;
	pAddedMass = new arma::cube* [numBodies];
	for (int ii=0; ii<numBodies; ii++)
	{
		pAddedMass[ii] = new arma::cube;
		added_mass_fn.str("");
		added_mass_fn << "body_" << this->GetId() << "/added_mass/body_" << ii;
		pAddedMass[ii]->load(arma::hdf5_name(filePath, added_mass_fn.str()));
	}
	
	// Read High frequency asymptotic added mass
	std::stringstream added_mass_hf_fn;
	pAddedMassHf = new arma::mat* [numBodies];
	for (int ii=0; ii<numBodies; ii++)
	{
		pAddedMassHf[ii] = new arma::mat;
		added_mass_hf_fn.str("");
		added_mass_hf_fn << "body_" << this->GetId() << "/added_mass_hf/body_" << ii;
		pAddedMassHf[ii]->load(arma::hdf5_name(filePath, added_mass_hf_fn.str()));
	}
	
	// Read Low frequency asymptotic added mass
	std::stringstream added_mass_lf_fn;
	pAddedMassLf = new arma::mat* [numBodies];
	for (int ii=0; ii<numBodies; ii++)
	{
		pAddedMassLf[ii] = new arma::mat;
		added_mass_lf_fn.str("");
		added_mass_lf_fn << "body_" << this->GetId() << "/added_mass_lf/body_" << ii;
		pAddedMassLf[ii]->load(arma::hdf5_name(filePath, added_mass_lf_fn.str()));
	}
	
	// Read wave radiation damping coeffients
	std::stringstream damping_radiation_fn;
	pDampingRadiation = new arma::cube* [numBodies];
	for (int ii=0; ii<numBodies; ii++)
	{
		pDampingRadiation[ii] = new arma::cube;
		damping_radiation_fn.str("");
		damping_radiation_fn << "body_" << this->GetId() << "/damping_radiation/body_" << ii;
		pDampingRadiation[ii]->load(arma::hdf5_name(filePath, damping_radiation_fn.str()));
	}
	
	// Read low frequency asymptotic wave radiation damping
	std::stringstream damping_radiation_lf_fn;
	pDampingRadiationLf = new arma::mat* [numBodies];
	for (int ii=0; ii<numBodies; ii++)
	{
		pDampingRadiationLf[ii] = new arma::mat;
		damping_radiation_lf_fn.str("");
		damping_radiation_lf_fn << "body_" << this->GetId() << "/damping_radiation_lf/body_" << ii;
		pDampingRadiationLf[ii]->load(arma::hdf5_name(filePath, damping_radiation_lf_fn.str()));
	}
	
	// Read mean drift coefficients
	std::stringstream mean_drift_fn;
	pMeanDrift = new arma::cube;
	mean_drift_fn << "body_" << this->GetId() << "/mean_drift";
	std::chrono::steady_clock::time_point begin_load = std::chrono::steady_clock::now();
	pMeanDrift->load(arma::hdf5_name(filePath, mean_drift_fn.str()));
	std::chrono::steady_clock::time_point end_load = std::chrono::steady_clock::now();
	
	// Read Wave exciting data
	std::stringstream wave_exciting_mag_fn;
	wave_exciting_mag_fn << "body_" << this->GetId() << "/wave_exciting_mag";
	pWaveExcitingMag = new arma::cube;
	pWaveExcitingMag->load(arma::hdf5_name(filePath, wave_exciting_mag_fn.str()));
	
	std::stringstream wave_exciting_pha_fn;
	pWaveExcitingPha = new arma::cube;
	wave_exciting_pha_fn << "body_" << this->GetId() << "/wave_exciting_pha";
	pWaveExcitingPha->load(arma::hdf5_name(filePath, wave_exciting_pha_fn.str()));
	
	// Read QTF data
	std::stringstream qtf_diff_fn;
	std::stringstream qtf_sum_fn;
	pQtfDiff = new arma::cube** [2];
	pQtfSum = new arma::cube** [2];
	for (int ii=0; ii<2; ii++)
	{
		pQtfDiff[ii] = new arma::cube* [activeDofs];
		pQtfSum[ii] = new arma::cube* [activeDofs];
		for (int jj=0; jj<activeDofs; jj++)
		{
			pQtfDiff[ii][jj] = new arma::cube;
			pQtfSum[ii][jj] = new arma::cube;
			
			qtf_diff_fn.str("");
			qtf_diff_fn << "body_" << this->GetId() << "/qtf_diff" << "/part_" << ii << "/dof_" << jj;
			pQtfDiff[ii][jj]->load(arma::hdf5_name(filePath, qtf_diff_fn.str()));
			
			qtf_sum_fn.str("");
			qtf_sum_fn << "body_" << this->GetId() << "/qtf_sum" << "/part_" << ii << "/dof_" << jj;
			pQtfSum[ii][jj]->load(arma::hdf5_name(filePath, qtf_sum_fn.str()));
		}
	}
	
}


/**
// Calcula la impulse response function
void Hydro::computeIRF(void){
	IRF = arma::zeros(6*nBodies,6*nBodies,nt_IRF);
}

// Calcula el espectro
void Hydro::computeWaveSpectrum(void){

	nComp = 10;
	wave_periods = arma::zeros(nComp,1); 
	wave_frequencies = arma::zeros(nComp,1); 
	wave_heights = arma::zeros(nComp,1); 
	wave_phases = arma::zeros(nComp,1); 
}

// Calcula la serie temporal de fuerzas de excitación
void Hydro::computeFe(void){

	Fe = arma::zeros(6*nBodies,nt_Fe);

}

// Obten las fuerzas hidroestaticas e hidrodinamicas en el tiempo deseado
void Hydro::computeHydroForces(double t){
	HydroForces = arma::zeros(6*nBodies,1);
	
	arma::mat positions = arma::zeros(6*nBodies,1);
	for(int ii=0;ii<nBodies;ii=ii+1){
		positions(arma::span(6*ii,6*(ii+1)-1),arma::span(0)) = Bodies[ii]->pos;
	}	

	HydroForces = HydroForces - hydro*positions;
	

}
**/