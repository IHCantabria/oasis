
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
#include "../os_tools.hpp"


arma::mat HydroDatabase::ComputeHydrostaticForces()
{	
	arma::mat hydrostatic_force = (*pHydrostaticStiffness)*((*pBodies)[id].pos - (*pBodies)[id].pos_init);
	return hydrostatic_force;
}


arma::mat HydroDatabase::ComputeRadiationForces()
{	
	// Allocate radiation force solution vector
	arma::mat radiation_force = arma::zeros(6, 1);

	// Declare local auxiliary variables
	arma::mat time_local;
	arma::mat time_irf;
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
				if (pSim->timeBuffer(0, pSim->timeBufferCount) > IRFTime(0, (*pIRFPoints[ib])(i, j)))
				{
					//std::cout << "HydroDatabase::ComputeRadiationForces - TimeBuffer exceed" << std::endl;
					//std::cout << "HydroDatabase::ComputeRadiationForces - timeBuffer size: " << arma::size(pSim->timeBuffer) << std::endl;
					//std::cout << "HydroDatabase::ComputeRadiationForces - IRFTime: " << IRFTime(0, (*pIRFPoints[ib])(i, j)) << std::endl;
					
					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);
					
					// Calculate begin and end indexes for matrix slicing
					for (int k=0; k<pSim->timeBufferCount; k++)
					{
						//std::cout << "here" << std::endl;
						if (pSim->timeBuffer(0, k) > (pSim->timeBuffer(0, pSim->timeBufferCount)-IRFTime(0, (*pIRFPoints[ib])(i, j))))
						{
							//std::cout << "here" << std::endl;
							(*pTimeStartPos)(ib, i, j) = k-1;
							//std::cout << "here" << std::endl;
							break;
						}
					}
					idx_begin = (*pTimeStartPos)(ib, i, j);
					idx_end = pSim->timeBufferCount;
					//std::cout << "HydroDatabase::ComputeRadiationForces - idx_begin: " << idx_begin << std::endl;
					//std::cout << "HydroDatabase::ComputeRadiationForces - idx_end: " << idx_end << std::endl;

					// Get velocity chunck from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, idx_begin, j, idx_end);

					// Interpolate local velocity vector in order to fit with IRF resolution
					//pSim->timeBuffer.cols(idx_begin, idx_end).print();
					//IRFTime.cols(0, (*pIRFPoints[ib])(i, j)).print();
					time_local = pSim->timeBuffer.cols(idx_begin, idx_end)-pSim->timeBuffer(0, idx_begin);
					time_irf = IRFTime.cols(0, (*pIRFPoints[ib])(i, j));
					vel_local_interp = interp1(time_local, vel_local, time_irf);

					// Calulate Duhamel integral
					vel_local_filter = arma::flipud(irf_local)%(vel_local_interp.t());
					vel_local_filter = vel_local_filter.t();
					radiation_force(i) += trapz(vel_local_filter, dt);
					/**
					if ((i==2) && (j==2))
					{
						std::cout << "dt: " << dt << std::endl;
						std::cout << "Radidation force (2,2): " << trapz(vel_local_filter, dt) << std::endl;
						std::cout << "Radidation force Unit time (2,2): " << trapz(vel_local_filter, 1.0) << std::endl;
						std::cout << "Radiation force: " << radiation_force(i) << std::endl;
						time_local.save(arma::hdf5_name("time_local_x.h5","irf"));
						time_irf.save(arma::hdf5_name("time_irf.h5","time_irf"));
						irf_local.save(arma::hdf5_name("irf_x.h5","irf"));
						vel_local.save(arma::hdf5_name("vel_x.h5","vel"));
						vel_local_interp.save(arma::hdf5_name("vel_local_interp_x.h5", "vel_local_interp"));
						vel_local_filter.save(arma::hdf5_name("vel_local_filter_x.h5","vel"));
						//time_local = (*SD.timeBuffer.cols(0, idx_end);
						//time_local.save(arma::hdf5_name("time_x.h5","time"));
					}
					**/
				}
				else if (pSim->timeBufferCount > 0)
				{
					// Calculate begin and end indexes for matrix slicing
					idx_end = pSim->timeBufferCount;

					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);

					// Get velocity chunk from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, 0, j, idx_end);

					// Interpolate local velocity vector in order to fit with IRF resolution
					time_local = IRFTime.cols(0, numPointsIRF-1) - IRFTime(0, numPointsIRF-1) + pSim->timeBuffer(0, idx_end);
					irf_local_interp = interp1(time_local, arma::flipud(irf_local), pSim->timeBuffer.cols(0, idx_end));

					// Calulate Duhamel integral
					vel_local_filter = irf_local_interp%vel_local;
					radiation_force(i) += trapzi(pSim->timeBuffer.cols(0, idx_end), vel_local_filter);

					/**
					if ((idx_end == 500) && (i==2) && (j==2))
					{
						time_local.save(arma::hdf5_name("time_local_2500.h5","irf"));
						irf_local.save(arma::hdf5_name("irf_2500.h5","irf"));
						irf_local_interp.save(arma::hdf5_name("irf_interp_2500.h5","irf_interp"));
						vel_local.save(arma::hdf5_name("vel_2500.h5","vel"));
						time_local = (*SD.timeBuffer).cols(0, idx_end);
						time_local.save(arma::hdf5_name("time_2500.h5","time"));
					}
					**/
				}
			}
		}
	}
	
	return radiation_force;
}


arma::mat HydroDatabase::ComputeFirstWaveExcForce()
{
	// Create local variables
	double time = pSim->pTimeSolver->t;
	// Get First Order Wave exciting data from storage
	arma::mat wave_exc_mag = pWaveExcitingMag->subcube(0, numPeriodExc, numHeadingExc, 5, numPeriodExc, numHeadingExc);
	arma::mat wave_exc_pha = pWaveExcitingPha->subcube(0, numPeriodExc, numHeadingExc, 5, numPeriodExc, numHeadingExc);
	// Calculate wave force
	double time_slope = 1/(*pFrequencies)(0, numPeriodExc);
	arma::mat wave_force = arma::zeros(6, 1);
	double angular_freq = 2*M_PI*(*pFrequencies)(0, numPeriodExc);
	for (int i=0; i<6; i++)
	{
		wave_force(i, 0) = wave_exc_mag(i, 0)*waveAmplitude*cos(angular_freq*time + wave_exc_pha(i, 0));

		if (time < time_slope)
		{
			wave_force(i, 0) = wave_force(i, 0)*time/time_slope;
		}
	}
	//std::cout << "Wave Amplitude: " << waveAmplitude << std::endl;

	return wave_force;
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
	double dt=0.01;
	double df = (*pFrequencies)(1)-(*pFrequencies)(0);
	double tmax = 1/df/2.0;
	IRFTime = arange(0, tmax, dt);
	
	// Allocate IRF matrix
	numPointsIRF = IRFTime.n_cols;
	pIRF = new arma::cube* [numBodies];
	pIRFPoints = new arma::mat* [numBodies];
	
	std::cout << "IRF Time: " << tmax << std::endl;
	std::cout << "IRF Time Points: " << numPointsIRF << std::endl;
	std::cout << "Frequency(0): " << (*pFrequencies)(0) << std::endl;
	std::cout << "Frequency(1): " << (*pFrequencies)(1) << std::endl;
	std::cout << "Frequency diff: " << df << std::endl;

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
				(*pIRF[ib])(0, i, j) = 2*trapz(dummy_mat, 2*M_PI*df)/M_PI;
				for (int k=1; k<IRFTime.n_cols; k++)
				{
					// Calculate new value of IRF
					dummy_mat = dampingFreq%cos(2*M_PI*(*pFrequencies)*IRFTime(0, k));
					(*pIRF[ib])(k, i, j) = 2*trapz(dummy_mat, 2*M_PI*df)/M_PI;
					
				}
			}
		}
	}

	(*pIRF[0]).save(arma::hdf5_name("IRF.h5", "irf"));

	std::cout << "Maximum retardation time: " << tmax << std::endl;
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time elapsed ComputeIRF: " << elapsed << std::endl;
}


int HydroDatabase::GetId(void)
{
	return id;
}


HydroDatabase::HydroDatabase(int incId, Body** incBody, Simulation* pIncSim): HydroForce()
{
	id = incId;
	pBodies = incBody;
	pSim = pIncSim;
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

	// Read position of the center of gravity
	std::stringstream cog_fn;
	cog_fn << "body_" << this->GetId() << "/cog";
	cog.load(arma::hdf5_name(filePath, cog_fn.str()));
	
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

	pTotalMass = new arma::mat(arma::size(*pStructuralMass), arma::fill::zeros);
	pTotalMass_inv = new arma::mat(arma::size(*pStructuralMass), arma::fill::zeros);
	*pTotalMass = *pStructuralMass;
	
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

		(*pTotalMass) = (*pTotalMass) + *pAddedMassHf[ii];
	}

	*pTotalMass_inv = arma::solve(*pTotalMass,eye(size(*pTotalMass)));
	
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
	
	// Generate starting pos time matrix
	pTimeStartPos = new arma::cube(numBodies, 6, 6, arma::fill::zeros);

	std::string pppPath = JoinPath(pSim->inputFolderPath, "periodNum.txt");
	char buffer_line [1000];
	FILE* p_file_pointer = fopen(pppPath.c_str(), "r");
	fscanf(p_file_pointer, "%lf %d %d %[^\n]\n", &waveAmplitude, &numPeriodExc, &numHeadingExc, buffer_line);
	fclose(p_file_pointer);
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