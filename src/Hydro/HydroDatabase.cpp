
#include <armadillo>
#include <iostream>
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
#include "../Waves/Wave.hpp"
#include "../Exceptions/Exception.hpp"


arma::mat HydroDatabase::CalculateHydrodynamicForces(double time)
{
	arma::mat F = arma::zeros(activeDofs,1);

	double rampa = std::min(1.0,time/10.0); // duración de la rampa harcodeado a 10s!!!!

	double yaw = pBodies[idBody]->pos(5,0);

	F = F + ComputeRadiationForces();

	if(pBodies[idBody]->firstOrderExcitationFlag==1)
	{
		F = F + ComputeFirstWaveExcForce(time)*rampa;
		//std::cout << "WARNING: Precomputed first order forces not implemented yet. \n" << std::endl;
	}
	if(pBodies[idBody]->firstOrderExcitationFlag==2)
	{
		F = F + ComputeFirstWaveExcForce(time)*rampa;
	}

	if(pBodies[idBody]->secondOrderExcitationFlag==1)
	{
		F = F + ComputeSecondWaveExcForce(time)*rampa;
		//std::cout << "WARNING: Precomputed second order forces not implemented yet. \n" << std::endl;
	}
	if(pBodies[idBody]->secondOrderExcitationFlag==2)
	{
		F = F + ComputeSecondWaveExcForce(time)*rampa;
	}
	if(pBodies[idBody]->secondOrderExcitationFlag==3)
	{
		F = F + F_meanDrift*rampa;
	}
	
	// Viscous drag forces
	arma::mat vv = pBodies[idBody]->vel;
	F = F - pBodies[idBody]->B_visc % vv
	      - pBodies[idBody]->B_visc2 % vv % arma::abs(vv);

	// Wind and current forces
	if (pMor->flag_wind)
	{
		F = F + pMor->ComputeWindForce(idBody, yaw, time);
	}
	if (pMor->flag_curr)
	{
		F = F + pMor->ComputeCurrForce(idBody, yaw, time);
	}

	return F;
}


arma::mat HydroDatabase::CalculateHydrostaticForces()
{	
	arma::mat hydrostatic_force = -(*pHydrostaticStiffness)*((*pBodies)[id].pos - (*pBodies)[id].pos_init);
	pBodies[idBody]->hydrostaticForces = hydrostatic_force;
	return hydrostatic_force;
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
	double df = (*pFrequencies)(1)-(*pFrequencies)(0);
	double tmax = 1/df/2.0;
	IRFTime = arange(0, 120.0, pSim->hydroTimeStep);
	
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

	//(*pIRF[0]).save(arma::hdf5_name("IRF.h5", "irf"));

	std::cout << "Maximum retardation time: " << tmax << std::endl;
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time elapsed ComputeIRF: " << elapsed << std::endl;
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
	double dt = pSim->hydroTimeStep;

	// Compute radiation forces for the 6DOFs
	for(int ib=0; ib<numBodies; ib++)
	{
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++)
			{
				if (pSim->timeBuffer(0, pSim->timeBufferCount) > 120.0)
				{	
					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);
					
					// Calculate begin and end indexes for matrix slicing
					idx_begin = pSim->timeBufferCount - (*pIRFPoints[ib])(i, j);
					idx_end = pSim->timeBufferCount;

					// Get velocity chunck from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, idx_begin, j, idx_end);

					// Calulate Duhamel integral
					vel_local_filter = arma::flipud(irf_local)%(vel_local.t());
					radiation_force(i) += trapz(vel_local_filter, dt);
					
				}
				else if (pSim->timeBufferCount > 0)
				{
					// Calculate begin and end indexes for matrix slicing
					idx_end = pSim->timeBufferCount;

					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, idx_end, i, j);

					// Get velocity chunk from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, 0, j, idx_end);

					// Calulate Duhamel integral
					vel_local_filter = arma::flipud(irf_local)%vel_local.t();
					radiation_force(i) += trapz(vel_local_filter, dt);

				}
			}
		}
	}

	radiation_force = -radiation_force;

	pBodies[idBody]->radiationForces = radiation_force;

	return radiation_force;
}


arma::mat HydroDatabase::GetCog()
{
	return this->cog;
}


int HydroDatabase::GetId(void)
{
	return id;
}


arma::mat HydroDatabase::GetInertiaMatrixInv(void)
{
	return *(this->pTotalMass_inv);
}


int HydroDatabase::GetNumBodies(void)
{
	return this->numBodies;
}


int HydroDatabase::GetNumPointsIrf(void)
{
	return this->numPointsIRF;
}


HydroDatabase::HydroDatabase(int incId, int incIdBody, Body** incBody, Simulation* pIncSim): HydroForce()
{
	id = incId;
	pBodies = incBody;
	pSim = pIncSim;
	idBody = incIdBody;
}


void HydroDatabase::LoadHydrodynamicData(std::string filePath)
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
	std::cout << "Reading frequencies...\n";
	std::stringstream frequencies_fn;
	pFrequencies = new arma::mat;
	frequencies_fn << "body_" << this->GetId() << "/frequencies";
	pFrequencies->load(arma::hdf5_name(filePath, frequencies_fn.str()));
	numFrequencies = pFrequencies->n_cols;
	
	// Read headings (radians)
	std::cout << "Reading headings...\n";
	std::stringstream headings_fn;
	pHeadings = new arma::mat;
	headings_fn << "body_" << this->GetId() << "/headings";
	pHeadings->load(arma::hdf5_name(filePath, headings_fn.str()));
	numHeadings = pHeadings->n_cols;
	
	// Read hydrostatic stiffness
	std::cout << "Reading hydrostatic matrix...\n";
	std::stringstream hydrostatic_stiffness_fn;
	pHydrostaticStiffness = new arma::mat;
	hydrostatic_stiffness_fn << "body_" << this->GetId() << "/hydstiffness";
	pHydrostaticStiffness->load(arma::hdf5_name(filePath, hydrostatic_stiffness_fn.str(), arma::hdf5_opts::trans));
	
	// Read structural mass properties
	std::cout << "Reading structural mass...\n";
	std::stringstream structural_mass_fn;
	pStructuralMass = new arma::mat;
	structural_mass_fn << "body_" << this->GetId() << "/mass";
	pStructuralMass->load(arma::hdf5_name(filePath, structural_mass_fn.str(), arma::hdf5_opts::trans));
	
	// Read Added Mass
	std::cout << "Reading Added Mass...\n";
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
	std::cout << "Reading high frequency...\n";
	std::stringstream added_mass_hf_fn;
	pAddedMassHf = new arma::mat* [numBodies];
	for (int ii=0; ii<numBodies; ii++)
	{
		pAddedMassHf[ii] = new arma::mat;
		added_mass_hf_fn.str("");
		added_mass_hf_fn << "body_" << this->GetId() << "/added_mass_hf/body_" << ii;
		pAddedMassHf[ii]->load(arma::hdf5_name(filePath, added_mass_hf_fn.str()));
		std::cout << "Applying matrix...\n";
		// (*pTotalMass)(arma::span(6*ii,6*(ii+1)-1), arma::span(6*ii,6*(ii+1)-1)) = (*pTotalMass)(arma::span(6*ii,6*(ii+1)-1), arma::span(6*ii,6*(ii+1)-1)) + *pAddedMassHf[ii];
	}
	
	// Read Low frequency asymptotic added mass
	std::cout << "Reading low frequency added mass...\n";
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
	std::cout << "Reading wave radiation damping...\n";
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
	std::cout << "Reading load frequency asymptotic wave radiation damping...\n";
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
	std::cout << "Reading mean drift coefficients...\n";
	std::stringstream mean_drift_fn;
	pMeanDrift = new arma::cube;
	mean_drift_fn << "body_" << this->GetId() << "/mean_drift";
	std::chrono::steady_clock::time_point begin_load = std::chrono::steady_clock::now();
	pMeanDrift->load(arma::hdf5_name(filePath, mean_drift_fn.str()));
	std::chrono::steady_clock::time_point end_load = std::chrono::steady_clock::now();
	
	// Read Wave exciting data
	std::cout << "Reading wave exciting data...\n";
	std::stringstream wave_exciting_mag_fn;
	wave_exciting_mag_fn << "body_" << this->GetId() << "/wave_exciting_mag";
	pWaveExcitingMag = new arma::cube;
	pWaveExcitingMag->load(arma::hdf5_name(filePath, wave_exciting_mag_fn.str()));
	
	std::stringstream wave_exciting_pha_fn;
	pWaveExcitingPha = new arma::cube;
	wave_exciting_pha_fn << "body_" << this->GetId() << "/wave_exciting_pha";
	pWaveExcitingPha->load(arma::hdf5_name(filePath, wave_exciting_pha_fn.str()));

	// Read QTF data
	std::cout << "Reading QTF data...\n";
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

	// Compute IRF function
	this->ComputeIRF();

	// Load Morison forces data
	pMor = new Morison(numBodies, pSim); pMor->ReadMorisonData();
}


arma::mat HydroDatabase::ComputeFirstWaveExcForce(double t)
{

	arma::mat Fe = arma::zeros(activeDofs,1);

	Wave* pWave = pSim->pWave;
	
	double x = pBodies[idBody]->pos(0,0);
	double y = pBodies[idBody]->pos(1,0); 
	double yaw = pBodies[idBody]->pos(5,0);

	if(pBodies[idBody]->secondOrderExcitationFlag==1)
	{
		x = 0; y = 0; yaw = 0;
	}

	//x = 0.0; y = 0.0; yaw = 0.0; // ELIMINO EN INSTANTO POSITION

	arma::cube H_Real = interp1((*pHeadings)+yaw,WE_Real_w,pWave->headings);

	//std::cout << "(*pHeadings)+yaw " << (*pHeadings)+yaw << std::endl;
	//std::cout << "pWave->headings " << pWave->headings << std::endl;
	//std::cout << "WE_Real_w " << WE_Real_w << std::endl;
	//std::cout << "H_Real " << H_Real << std::endl;

	arma::cube H_Imag = interp1((*pHeadings)+yaw,WE_Imag_w,pWave->headings);
	arma::cube H_Mag = arma::sqrt(arma::pow(H_Real,2)+arma::pow(H_Imag,2));
	arma::cube H_Pha = arma::atan2(H_Imag,H_Real);
	H_Mag = permute(H_Mag,213); H_Pha = permute(H_Pha,213);

	arma::mat H_Mag_loc, H_Pha_loc, Cm, Sm, Am, PHIm;

	//std::cout << "H_Real " << H_Real << std::endl;
	//std::cout << "H_Imag " << H_Imag << std::endl;
	//std::cout << "H_Mag " << H_Mag << std::endl;

	for(int ii=0; ii<activeDofs; ii++) // Por ahora generico para 6 dofs, esto habría que cambiarlo
	{
		H_Mag_loc = H_Mag(arma::span::all,arma::span::all,arma::span(ii));
		H_Pha_loc = H_Pha(arma::span::all,arma::span::all,arma::span(ii));

		Cm = arma::sum(H_Mag_loc % pWave->amplitudes % arma::cos(pWave->phases + H_Pha_loc + x*pWave->kx + y*pWave->ky),1);
		Sm = arma::sum(H_Mag_loc % pWave->amplitudes % arma::sin(pWave->phases + H_Pha_loc + x*pWave->kx + y*pWave->ky),1);

		Am = arma::sqrt(arma::pow(Cm,2)+arma::pow(Sm,2));
		PHIm = arma::atan2(Sm,Cm);

		Fe(ii,0) = arma::as_scalar(arma::sum(Am % arma::cos(t * pWave->ang_freqs - PHIm),0));
	}

	//std::cout << "t " << t << std::endl;
	//std::cout << "Fe " << Fe << std::endl;

	double Fx = arma::as_scalar(Fe(0,0)); double Fy = arma::as_scalar(Fe(1,0));
	Fe(0,0) = Fx*cos(yaw)-Fy*sin(yaw); Fe(1,0) = Fx*sin(yaw)+Fy*cos(yaw);

	pBodies[idBody]->excitationForces_1 = Fe;

	return Fe;
}


arma::mat HydroDatabase::ComputeSecondWaveExcForce(double t)
{

	arma::mat Fe = arma::zeros(activeDofs,1);

	Wave* pWave = pSim->pWave;
	
	double x = pBodies[idBody]->pos(0,0);
	double y = pBodies[idBody]->pos(1,0); 
	double yaw = pBodies[idBody]->pos(5,0);

	if(pBodies[idBody]->secondOrderExcitationFlag==1)
	{
		x = 0; y = 0; yaw = 0;
	}

	arma::cube temp1;
	arma::cube*** HDif = new arma::cube** [2];
	arma::cube*** HSum = new arma::cube** [2];
	for (int ii=0; ii<2; ii++)
	{
		HDif[ii] = new arma::cube* [activeDofs];
		HSum[ii] = new arma::cube* [activeDofs];
		for (int jj=0; jj<activeDofs; jj++)
		{
			HDif[ii][jj] = new arma::cube;
			HSum[ii][jj] = new arma::cube;

			temp1 = *QtfDiff_w[ii][jj]; 
			temp1 = interp1((*pHeadings)+yaw, permute(temp1,312), pWave->headings_1D);
			*HDif[ii][jj] = permute(temp1,231);

			temp1 = *QtfSum_w[ii][jj];
			temp1 = interp1((*pHeadings)+yaw, permute(temp1,312), pWave->headings_1D);
			*HSum[ii][jj] = permute(temp1,231);
		}
	}

	arma::mat temp2, temp3;
	for(int ii=0; ii<activeDofs; ii++)
	{
		temp2 = (*HDif[0][ii]).slice(0);
		temp3 = temp2%ampP*arma::cos(wD*t+phD+kxD*x+kyD*y);
		Fe(ii,0) = Fe(ii,0) + 0.5*arma::as_scalar(arma::sum(arma::sum(temp3,1),0));

		temp2 = (*HDif[1][ii]).slice(0);
		temp3 = temp2%ampP*arma::sin(wD*t+phD+kxD*x+kyD*y);
		Fe(ii,0) = Fe(ii,0) + 0.5*arma::as_scalar(arma::sum(arma::sum(temp3,1),0));

		temp2 = (*HSum[0][ii]).slice(0);
		temp3 = temp2%ampP*arma::cos(wS*t+phS+kxS*x+kyS*y);
		Fe(ii,0) = Fe(ii,0) + 0.5*arma::as_scalar(arma::sum(arma::sum(temp3,1),0));

		temp2 = (*HSum[1][ii]).slice(0);
		temp3 = temp2%ampP*arma::sin(wS*t+phS+kxS*x+kyS*y);
		Fe(ii,0) = Fe(ii,0) + 0.5*arma::as_scalar(arma::sum(arma::sum(temp3,1),0));
	}

	pBodies[idBody]->excitationForces_2 = Fe;
	double Fx = arma::as_scalar(Fe(0,0)); double Fy = arma::as_scalar(Fe(1,0));
	Fe(0,0) = Fx*cos(yaw)-Fy*sin(yaw); Fe(1,0) = Fx*sin(yaw)+Fy*cos(yaw);

	return Fe;
}


void HydroDatabase::SetUp(void)
{
	// Add the zero and infinite frequencies data to de matrices and cubes

	/*
		NOT PROPERLY IMPLEMENTED
	*/


	// Interpolate from the hidrodatabase frequencies to the wave frequencies
	Wave* pWave = pSim->pWave;
	arma::cube WE_Real = (*pWaveExcitingMag) % arma::cos((*pWaveExcitingPha));
	arma::cube WE_Imag = (*pWaveExcitingMag) % arma::sin((*pWaveExcitingPha));

	WE_Real_w = interp1(*pFrequencies,permute(WE_Real,231), pWave->freqs);
	WE_Imag_w = interp1(*pFrequencies,permute(WE_Imag,231), pWave->freqs);
	WE_Real_w = permute(WE_Real_w,213); WE_Imag_w = permute(WE_Imag_w,213);

	arma::cube temp;
	QtfDiff_w = new arma::cube** [2];
	QtfSum_w = new arma::cube** [2];
	for (int ii=0; ii<2; ii++)
	{
		QtfDiff_w[ii] = new arma::cube* [activeDofs];
		QtfSum_w[ii] = new arma::cube* [activeDofs];
		for (int jj=0; jj<activeDofs; jj++)
		{
			QtfDiff_w[ii][jj] = new arma::cube;
			QtfSum_w[ii][jj] = new arma::cube;

			temp = *pQtfDiff[ii][jj]; 
			*QtfDiff_w[ii][jj] = interp2(*pFrequencies, *pFrequencies, temp, pWave->freqs, pWave->freqs);

			temp = *pQtfSum[ii][jj]; 
			*QtfSum_w[ii][jj] = interp2(*pFrequencies, *pFrequencies, temp, pWave->freqs, pWave->freqs);
		}
	}

	ampP = arma::zeros(pWave->num_comps,pWave->num_comps);
	wS = arma::zeros(pWave->num_comps,pWave->num_comps);
	phS = arma::zeros(pWave->num_comps,pWave->num_comps);
	kxS = arma::zeros(pWave->num_comps,pWave->num_comps);
	kyS = arma::zeros(pWave->num_comps,pWave->num_comps);
	wD = arma::zeros(pWave->num_comps,pWave->num_comps);
	phD = arma::zeros(pWave->num_comps,pWave->num_comps);
	kxD = arma::zeros(pWave->num_comps,pWave->num_comps);
	kyD = arma::zeros(pWave->num_comps,pWave->num_comps);

	for (int ii=0; ii<pWave->num_comps; ii++)
	{
		for (int jj=0; jj<pWave->num_comps; jj++)
		{
			ampP(ii,jj) = pWave->amplitudes_1D(ii,0)*pWave->amplitudes_1D(ii,0);
			wS(ii,jj) = pWave->ang_freqs(ii,0)+pWave->ang_freqs(ii,0);
			phS(ii,jj) = pWave->phases_1D(ii,0)+pWave->phases_1D(ii,0);
			kxS(ii,jj) = pWave->kx_1D(ii,0)+pWave->kx_1D(ii,0);
			kyS(ii,jj) = pWave->ky_1D(ii,0)+pWave->ky_1D(ii,0);
			wD(ii,jj) = pWave->ang_freqs(ii,0)-pWave->ang_freqs(ii,0);
			phD(ii,jj) = pWave->phases_1D(ii,0)-pWave->phases_1D(ii,0);
			kxD(ii,jj) = pWave->kx_1D(ii,0)-pWave->kx_1D(ii,0);
			kyD(ii,jj) = pWave->ky_1D(ii,0)-pWave->ky_1D(ii,0);
		}
	}

	// Include viscous added mass
	for (int ii=0; ii<numBodies; ii++)
	{
		(*pTotalMass)(arma::span(6*ii,6*(ii+1)-1), arma::span(6*ii,6*(ii+1)-1)) = 
		(*pTotalMass)(arma::span(6*ii,6*(ii+1)-1), arma::span(6*ii,6*(ii+1)-1)) 
		+ (*pAddedMassHf[ii])%(arma::diagmat(pBodies[ii]->A_visc));
	}
	*pTotalMass_inv = arma::solve(*pTotalMass,eye(size(*pTotalMass)));

	// Compute mean drift force
	arma::cube temp_mD = interp2(*pFrequencies, *pHeadings, permute(*pMeanDrift,231), pWave->freqs, pWave->headings);
	for (int ii=0; ii<activeDofs; ii++)
	{
		temp_mD.slice(ii) = temp_mD.slice(ii)%pWave->amplitudes%pWave->amplitudes;
	}
	temp_mD = permute(temp_mD,312);
	F_meanDrift = arma::sum(arma::sum(temp_mD,1),2);
	if(pBodies[idBody]->secondOrderExcitationFlag==3)
	{
		pBodies[idBody]->excitationForces_2 = F_meanDrift;
	}	

	pBodies[idBody]->Fb = CalculateHydrodynamicForces(0.0);

	if(pBodies[idBody]->firstOrderExcitationFlag==1)
	{
		std::cout << "WARNING: Precomputed first order forces not implemented yet. \n" 
		             "         Using limited instanto position version instead. \n"<< std::endl;
	}
	if(pBodies[idBody]->secondOrderExcitationFlag==1)
	{
		std::cout << "WARNING: Precomputed second order forces not implemented yet. \n" 
		             "         Using limited instanto position version instead. \n"<< std::endl;
	}
}


void HydroDatabase::Print()
{
	std::cout << "Number of bodies associated: " << numBodies << std::endl;
	std::cout << "Number of frequencies: " << numFrequencies << std::endl;
	std::cout << "Number of headings: " << numHeadings << std::endl;
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