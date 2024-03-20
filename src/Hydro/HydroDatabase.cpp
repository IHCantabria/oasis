
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
	arma::mat F = arma::zeros(activeDofs, 1);

	double rampa = std::min(1.0, time / 50.0); // duración de la rampa harcodeado a 10s!!!!

	double yaw = pBodies[idBody]->pos(5, 0);

	F = F + ComputeRadiationForces();

	if (pBodies[idBody]->firstOrderExcitationFlag == 1)
	{
		F = F + ComputeFirstWaveExcForce(time) * rampa;
		// std::cout << "WARNING: Precomputed first order forces not implemented yet. \n" << std::endl;
	}
	if (pBodies[idBody]->firstOrderExcitationFlag == 2)
	{
		F = F + ComputeFirstWaveExcForce(time) * rampa;
	}

	if (pBodies[idBody]->secondOrderExcitationFlag == 1)
	{
		F = F + ComputeSecondWaveExcForce(time) * rampa;
		// std::cout << "WARNING: Precomputed second order forces not implemented yet. \n" << std::endl;
	}
	if (pBodies[idBody]->secondOrderExcitationFlag == 2)
	{
		F = F + ComputeSecondWaveExcForce(time) * rampa;
	}
	if (pBodies[idBody]->secondOrderExcitationFlag == 3 || pBodies[idBody]->secondOrderExcitationFlag == 4)
	{
		F = F + ComputeMeanDrift() * rampa;
	}

	// Viscous drag forces
	arma::mat vv = pBodies[idBody]->vel;
	F = F - pBodies[idBody]->B_visc % vv - pBodies[idBody]->B_visc2 % vv % arma::abs(vv);

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

arma::mat HydroDatabase::CalculateHydrostaticForces(double time)
{

	arma::mat hydrostatic_force;
	if (pBodies[idBody]->flag_hidrostatics == 0)
	{

		// Linear hydrostatic forces
		hydrostatic_force = -(*pHydrostaticStiffness) * (pBodies[idBody]->pos - pBodies[idBody]->pos_eq);
		if (pBodies[idBody]->filling_mass > 0)
		{
			arma::mat Fg = arma::zeros(3, 1);
			Fg(2, 0) = -pSim->gravity * pBodies[idBody]->filling_mass;
			hydrostatic_force.rows(0, 2) = hydrostatic_force.rows(0, 2) + Fg;
			hydrostatic_force.rows(3, 5) = hydrostatic_force.rows(3, 5) + arma::cross(pBodies[idBody]->pos_filling_cog, pBodies[idBody]->rotMat.t() * Fg);
		}
	}
	else if (pBodies[idBody]->flag_hidrostatics > 0)
	{

		// Non-linear hydrostatic forces
		pBodies[idBody]->pNLHSMesh->TransformMesh();
		pBodies[idBody]->pNLHSMesh->CutMesh(time);
		pBodies[idBody]->pNLHSMesh->IntegrateMesh();
		arma::mat pressure = CalculateHydrostaticPressure(time);

		arma::mat radius = pBodies[idBody]->pNLHSMesh->nodes -
						   arma::ones(pBodies[idBody]->pNLHSMesh->numNodes, 1) * pBodies[idBody]->pos.rows(0, 2).t();

		hydrostatic_force = arma::zeros(6, 1);
		arma::uvec i1_vec = {1, 2, 0}, i2_vec = {2, 0, 1};
		arma::mat weightsJacNormal = pBodies[idBody]->pNLHSMesh->weightsJacNormal;

		for (int i = 0; i < 3; i++)
		{
			hydrostatic_force(i, 0) = arma::dot(-pressure, weightsJacNormal.col(i));
			arma::uword i1 = i1_vec(i);
			arma::uword i2 = i2_vec(i);
			hydrostatic_force(i + 3, 0) = arma::dot(-pressure, (radius.col(i1) % weightsJacNormal.col(i2) - radius.col(i2) % weightsJacNormal.col(i1)));
		}

		// Add gravity force
		hydrostatic_force(2, 0) = hydrostatic_force(2, 0) - pSim->gravity * pBodies[idBody]->structuralMass;

		// Transform the moments into the local frame
		hydrostatic_force.rows(3, 5) = pBodies[idBody]->rotMat.t() * hydrostatic_force.rows(3, 5);
	}

	pBodies[idBody]->hydrostaticForces = hydrostatic_force;
	return hydrostatic_force;
}

void HydroDatabase::ComputeIRF(std::string HDBname)
{
	// Declare local variables
	int count_max;
	int count_zero;
	int max_consec;
	arma::mat dummy_mat;
	arma::mat max_position;
	arma::mat zero_cross;
	arma::mat dampingFreq;
	arma::mat frequencies_trapz;
	arma::mat dampingFreq_trapz;

	// Check input arguments
	if ((*pFrequencies)(1) < (*pFrequencies)(0))
	{
		perror("Frequencies does not increase monotonically...\n");
	}

	// Calculate maximum time allowed
	IRFTime = arange(0, IRFTotalTime, pSim->hydroTimeStep);

	// Allocate IRF matrix
	numPointsIRF = IRFTime.n_cols;
	pIRF = new arma::cube *[numBodies];
	pIRFPoints = new arma::mat *[numBodies];

	std::cout << "    IRF Time Points: " << numPointsIRF << std::endl;
	// std::cout << "    Frequency(0): " << (*pFrequencies)(0) << std::endl;
	// std::cout << "    Frequency(1): " << (*pFrequencies)(1) << std::endl;

	double df = 1.0 / (2.0 * IRFTotalTime);

	frequencies_trapz = arma::regspace(std::max((*pFrequencies).min(), df), df, std::min(1.0 / (2.0 * pSim->hydroTimeStep), (*pFrequencies).max())).t();

	// Loop to find the IRF value for each body influence and DOF
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	for (int ib = 0; ib < numBodies; ib++)
	{
		pIRFPoints[ib] = new arma::mat(6, 6, arma::fill::zeros);
		pIRF[ib] = new arma::cube(IRFTime.n_cols, 6, 6, arma::fill::zeros);
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				// Clear previous results
				count_max = 1;
				count_zero = -1;
				max_consec = 0;
				max_position = arma::zeros(1, IRFTime.n_cols);
				zero_cross = arma::zeros(1, IRFTime.n_cols);

				// Start new Dof data
				dampingFreq = (*pDampingRadiation[ib]).subcube(i, j, 0, i, j, numFrequencies - 1);
				arma::interp1(*pFrequencies, dampingFreq, frequencies_trapz, dampingFreq_trapz);
				dummy_mat = dampingFreq_trapz % cos(2 * arma::datum::pi * frequencies_trapz * IRFTime(0, 0));
				(*pIRF[ib])(0, i, j) = 2 * trapzi(2 * arma::datum::pi * frequencies_trapz, dummy_mat) / arma::datum::pi;
				(*pIRFPoints[ib])(i, j) = IRFTime.n_cols - 1;
				for (int k = 1; k < IRFTime.n_cols; k++)
				{
					// Calculate new value of IRF
					dummy_mat = dampingFreq_trapz % cos(2 * arma::datum::pi * frequencies_trapz * IRFTime(0, k));
					(*pIRF[ib])(k, i, j) = 2 * trapzi(2 * arma::datum::pi * frequencies_trapz, dummy_mat) / arma::datum::pi;
				}
			}
		}
	}

	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "    Time elapsed ComputeIRF: " << elapsed << std::endl;

	char buffer[50];
	for (int ib = 0; ib < pSim->numBodies; ib++)
	{
		// sprintf(buffer,"IRF_Body_%d_fromBody_%d.dat", pBodies[idBody]->GetId(), pBodies[ib]->GetId());
		std::string filename = JoinPath(pSim->outputFolderPath, HDBname);
		filename = filename + "_IRF_Body_" + std::to_string(pBodies[idBody]->GetId()) +
				   "_fromBody_" + std::to_string(pBodies[ib]->GetId()) + ".dat";
		pIRF[ib]->save(filename, arma::arma_ascii);
	}
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
	for (int ib = 0; ib < numBodies; ib++)
	{
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (pSim->timeBuffer(0, pSim->timeBufferCount) > IRFTotalTime)
				{
					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF - 1, i, j);

					// Calculate begin and end indexes for matrix slicing
					idx_begin = pSim->timeBufferCount - (*pIRFPoints[ib])(i, j);
					idx_end = pSim->timeBufferCount;

					// Get velocity chunck from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, idx_begin, j, idx_end);

					// Calulate Duhamel integral
					vel_local_filter = arma::flipud(irf_local) % (vel_local.t());
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
					vel_local_filter = arma::flipud(irf_local) % vel_local.t();
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

int HydroDatabase::GetNumBodies(void)
{
	return this->numBodies;
}

int HydroDatabase::GetNumPointsIrf(void)
{
	return this->numPointsIRF;
}

arma::mat HydroDatabase::GetTotalMass(void)
{
	return (*pTotalMass);
}

void HydroDatabase::UpdateStructuralMass(arma::mat newStructuralMass)
{
	*pStructuralMass = newStructuralMass;
}

void HydroDatabase::UpdateTotalMass(void)
{
	(*pTotalMass)(arma::span(0, 5), arma::span(6 * (pBodies[idBody]->hydroDatabaseIndex), 6 * (pBodies[idBody]->hydroDatabaseIndex + 1) - 1)) = (*pStructuralMass);
	for (int ii = 0; ii < numBodies; ii++)
	{
		(*pTotalMass)(arma::span(0, 5), arma::span(6 * ii, 6 * (ii + 1) - 1)) += (*pAddedMassHf[ii]);
	}
	(*pTotalMass)(arma::span(0, 5), arma::span(6 * (pBodies[idBody]->hydroDatabaseIndex), 6 * (pBodies[idBody]->hydroDatabaseIndex + 1) - 1)) +=
		(*pAddedMassHf[pBodies[idBody]->hydroDatabaseIndex]) % (arma::diagmat(pBodies[idBody]->A_visc));
}

HydroDatabase::HydroDatabase(int incId, int incIdBody, Body **incBody, Simulation *pIncSim) : HydroForce()
{
	id = incId;
	pBodies = incBody;
	pSim = pIncSim;
	idBody = incIdBody;
	IRFTotalTime = pSim->timeIRF;
}

void HydroDatabase::LoadHydrodynamicData(std::string filePath)
{
	// Read number of bodies
	std::cout << "  Reading number of bodies...\n";
	arma::mat num_bodies_mat;
	num_bodies_mat.load(arma::hdf5_name(filePath, "num_bodies"));
	numBodies = num_bodies_mat(0);

	// Read position of the center of gravity
	std::cout << "  Reading position of the center of gravity...\n";
	std::stringstream cog_fn;
	cog_fn << "body_" << this->GetId() << "/cog";
	cog.load(arma::hdf5_name(filePath, cog_fn.str()));

	// Read frequencies
	std::cout << "  Reading frequencies...\n";
	std::stringstream frequencies_fn;
	pFrequencies = new arma::vec;
	frequencies_fn << "body_" << this->GetId() << "/frequencies";
	pFrequencies->load(arma::hdf5_name(filePath, frequencies_fn.str()));
	numFrequencies = pFrequencies->n_rows;

	// Read headings (radians)
	std::cout << "  Reading headings...\n";
	std::stringstream headings_fn;
	pHeadings = new arma::vec;
	headings_fn << "body_" << this->GetId() << "/headings";
	pHeadings->load(arma::hdf5_name(filePath, headings_fn.str()));
	numHeadings = pHeadings->n_rows;

	// Read hydrostatic stiffness
	std::cout << "  Reading hydrostatic matrix...\n";
	std::stringstream hydrostatic_stiffness_fn;
	pHydrostaticStiffness = new arma::mat;
	hydrostatic_stiffness_fn << "body_" << this->GetId() << "/hydstiffness";
	pHydrostaticStiffness->load(arma::hdf5_name(filePath, hydrostatic_stiffness_fn.str(), arma::hdf5_opts::trans));

	// Read structural mass properties
	std::cout << "  Reading structural mass...\n";
	std::stringstream structural_mass_fn;
	pStructuralMass = new arma::mat;
	structural_mass_fn << "body_" << this->GetId() << "/mass";
	pStructuralMass->load(arma::hdf5_name(filePath, structural_mass_fn.str(), arma::hdf5_opts::trans));

	// Load the structural mass intr
	pBodies[idBody]->structuralMass = (*pStructuralMass)(0, 0);

	// Create total mass matrix and fill with structural data
	std::cout << "  Creating total mass matrix...\n";
	pTotalMass = new arma::mat(6, 6 * numBodies, arma::fill::zeros);
	(*pTotalMass)(arma::span(0, 5), arma::span(6 * (pBodies[idBody]->hydroDatabaseIndex), 6 * (pBodies[idBody]->hydroDatabaseIndex + 1) - 1)) = (*pStructuralMass);

	// Read Added Mass
	std::cout << "  Reading Added Mass...\n";
	std::stringstream added_mass_fn;
	pAddedMass = new arma::cube *[numBodies];
	for (int ii = 0; ii < numBodies; ii++)
	{
		pAddedMass[ii] = new arma::cube;
		added_mass_fn.str("");
		added_mass_fn << "body_" << this->GetId() << "/added_mass/body_" << ii;
		pAddedMass[ii]->load(arma::hdf5_name(filePath, added_mass_fn.str()));
	}

	// Read High frequency asymptotic added mass
	std::cout << "  Reading high frequency...\n";
	std::stringstream added_mass_hf_fn;
	pAddedMassHf = new arma::mat *[numBodies];
	for (int ii = 0; ii < numBodies; ii++)
	{
		pAddedMassHf[ii] = new arma::mat;
		added_mass_hf_fn.str("");
		added_mass_hf_fn << "body_" << this->GetId() << "/added_mass_hf/body_" << ii;
		pAddedMassHf[ii]->load(arma::hdf5_name(filePath, added_mass_hf_fn.str()));
		std::cout << "    Applying matrix...\n";
		std::cout << "    " << 6 * ii << " - " << 6 * (ii + 1) - 1 << "\n";
		(*pTotalMass)(arma::span(0, 5), arma::span(6 * ii, 6 * (ii + 1) - 1)) += (*pAddedMassHf[ii]);
	}
	(*pTotalMass)(arma::span(0, 5), arma::span(6 * (pBodies[idBody]->hydroDatabaseIndex), 6 * (pBodies[idBody]->hydroDatabaseIndex + 1) - 1)) +=
		(*pAddedMassHf[pBodies[idBody]->hydroDatabaseIndex]) % (arma::diagmat(pBodies[idBody]->A_visc));

	// Read Low frequency asymptotic added mass
	std::cout << "  Reading low frequency added mass...\n";
	std::stringstream added_mass_lf_fn;
	pAddedMassLf = new arma::mat *[numBodies];
	for (int ii = 0; ii < numBodies; ii++)
	{
		pAddedMassLf[ii] = new arma::mat;
		added_mass_lf_fn.str("");
		added_mass_lf_fn << "body_" << this->GetId() << "/added_mass_lf/body_" << ii;
		pAddedMassLf[ii]->load(arma::hdf5_name(filePath, added_mass_lf_fn.str()));
	}

	// Read wave radiation damping coeffients
	std::cout << "  Reading wave radiation damping...\n";
	std::stringstream damping_radiation_fn;
	pDampingRadiation = new arma::cube *[numBodies];
	for (int ii = 0; ii < numBodies; ii++)
	{
		pDampingRadiation[ii] = new arma::cube;
		damping_radiation_fn.str("");
		damping_radiation_fn << "body_" << this->GetId() << "/damping_radiation/body_" << ii;
		pDampingRadiation[ii]->load(arma::hdf5_name(filePath, damping_radiation_fn.str()));
	}

	// Read low frequency asymptotic wave radiation damping
	std::cout << "  Reading load frequency asymptotic wave radiation damping...\n";
	std::stringstream damping_radiation_lf_fn;
	pDampingRadiationLf = new arma::mat *[numBodies];
	for (int ii = 0; ii < numBodies; ii++)
	{
		pDampingRadiationLf[ii] = new arma::mat;
		damping_radiation_lf_fn.str("");
		damping_radiation_lf_fn << "body_" << this->GetId() << "/damping_radiation_lf/body_" << ii;
		pDampingRadiationLf[ii]->load(arma::hdf5_name(filePath, damping_radiation_lf_fn.str()));
	}

	// Read Wave exciting data
	std::cout << "  Reading wave exciting data...\n";
	std::stringstream wave_exciting_mag_fn;
	if (pBodies[idBody]->flag_hidrostatics == 2)
	{
		wave_exciting_mag_fn << "body_" << this->GetId() << "/wave_diffraction_mag";
	}
	else
	{
		wave_exciting_mag_fn << "body_" << this->GetId() << "/wave_exciting_mag";
	}
	pWaveExcitingMag = new arma::cube;
	pWaveExcitingMag->load(arma::hdf5_name(filePath, wave_exciting_mag_fn.str()));

	std::stringstream wave_exciting_pha_fn;
	pWaveExcitingPha = new arma::cube;
	if (pBodies[idBody]->flag_hidrostatics == 2)
	{
		wave_exciting_pha_fn << "body_" << this->GetId() << "/wave_diffraction_pha";
	}
	else
	{
		wave_exciting_pha_fn << "body_" << this->GetId() << "/wave_exciting_pha";
	}
	pWaveExcitingPha->load(arma::hdf5_name(filePath, wave_exciting_pha_fn.str()));

	if (pBodies[idBody]->secondOrderExcitationFlag > 0 && pBodies[idBody]->secondOrderExcitationFlag < 4)
	{
		// Read QTF data
		std::cout << "  Reading QTF data...\n";
		std::stringstream qtf_diff_fn;
		std::stringstream qtf_sum_fn;
		pQtfDiff = new arma::cube **[2];
		pQtfSum = new arma::cube **[2];
		for (int ii = 0; ii < 2; ii++)
		{
			pQtfDiff[ii] = new arma::cube *[activeDofs];
			pQtfSum[ii] = new arma::cube *[activeDofs];
			for (int jj = 0; jj < activeDofs; jj++)
			{
				pQtfDiff[ii][jj] = new arma::cube;
				pQtfSum[ii][jj] = new arma::cube;

				qtf_diff_fn.str("");
				qtf_diff_fn << "body_" << this->GetId() << "/qtf_diff"
							<< "/part_" << ii << "/dof_" << jj;
				pQtfDiff[ii][jj]->load(arma::hdf5_name(filePath, qtf_diff_fn.str()));

				qtf_sum_fn.str("");
				qtf_sum_fn << "body_" << this->GetId() << "/qtf_sum"
						   << "/part_" << ii << "/dof_" << jj;
				pQtfSum[ii][jj]->load(arma::hdf5_name(filePath, qtf_sum_fn.str()));
			}
		}

		pMeanDrift = new arma::cube(activeDofs, numFrequencies, numHeadings, arma::fill::zeros);
		int my_count = 0;
		double temp_value = 0.0;

		for (int ii = 0; ii < activeDofs; ii++)
		{
			for (int jj = 0; jj < numFrequencies; jj++)
			{
				for (int kk = 0; kk < numHeadings; kk++)
				{
					(*pMeanDrift)(ii, jj, kk) = (*pQtfDiff[0][ii])(jj, jj, kk);
				}
			}
		}
	}

	if (pBodies[idBody]->secondOrderExcitationFlag == 4)
	{
		// Read mean drift coefficients
		std::cout << "Reading mean drift coefficients...\n";
		std::stringstream mean_drift_fn;
		pMeanDrift = new arma::cube;
		mean_drift_fn << "body_" << this->GetId() << "/mean_drift";
		std::chrono::steady_clock::time_point begin_load = std::chrono::steady_clock::now();
		pMeanDrift->load(arma::hdf5_name(filePath, mean_drift_fn.str()));
		std::chrono::steady_clock::time_point end_load = std::chrono::steady_clock::now();
	}

	// Generate starting pos time matrix
	pTimeStartPos = new arma::cube(numBodies, 6, 6, arma::fill::zeros);

	// Compute IRF function
	std::cout << "  Computing IRF ...\n";
	std::string HDBname = filePath.substr(filePath.find_last_of("/") + 1);
	HDBname = HDBname.substr(0, HDBname.length() - 6);
	this->ComputeIRF(HDBname);

	// Load Morison forces data
	std::cout << "  Reading Morison forces data ...\n";
	pMor = new Morison(numBodies, pSim);
	pMor->ReadMorisonData();
}

arma::mat HydroDatabase::ComputeFirstWaveExcForce(double t)
{
	// Initialize excitation force vector
	arma::mat Fe = arma::zeros(activeDofs, 1);

	// Load wave pointer as local variable
	Wave *pWave = pSim->pWave;

	// Extract position and heading of the body
	double x = pBodies[idBody]->pos(0, 0);
	double y = pBodies[idBody]->pos(1, 0);
	double yaw = pBodies[idBody]->pos(5, 0);

	// TODO: Remove this when the precomputed first order excitation forces are implemented
	if (pBodies[idBody]->firstOrderExcitationFlag == 1)
	{
		x = 0;
		y = 0;
		yaw = 0;
	}

	// Interpolate transfer functions to the body heading
	arma::cube H_Real = interp1((*pHeadings) + yaw, WE_Real_w, pWave->headings_piece);
	arma::cube H_Imag = interp1((*pHeadings) + yaw, WE_Imag_w, pWave->headings_piece);

	// Convert interpolated transfer functions to magnitude and phase
	arma::cube H_Mag = arma::sqrt(arma::pow(H_Real, 2) + arma::pow(H_Imag, 2));
	arma::cube H_Pha = arma::atan2(H_Imag, H_Real);

	// Permute the transfer functions to the correct order (frequencies, headings, dofs)
	H_Mag = permute(H_Mag, 213);
	H_Pha = permute(H_Pha, 213);

	// Extract the index of the required wave piece
	arma::uvec ind_piece = arma::find((pWave->time_ini - pWave->time_gap) < t, 1, "last");
	// Check if current time is in a gap
	bool flag_gap = ((ind_piece(0) > 0) && (t < pWave->time_ini(ind_piece(0))));

	// Declare local variables
	arma::mat H_Mag_loc, H_Pha_loc, Cm, Sm, Am, PHIm;
	double F_piece, F_gap;

	// Loop over all active degrees of freedom
	// TODO: This only works for 6 DOFs, change it.
	for (int ii = 0; ii < activeDofs; ii++)
	{
		// Extract transfer functions for the current degree of freedom
		H_Mag_loc = H_Mag(arma::span::all, arma::span::all, arma::span(ii));
		H_Pha_loc = H_Pha(arma::span::all, arma::span::all, arma::span(ii));

		// Multiply the transfer functions with the wave amplitudes for the current wave piece
		Cm = arma::sum(H_Mag_loc % pWave->amplitudes_piece(ind_piece(0)) %
						   arma::cos(pWave->phases_piece(ind_piece(0)) +
									 H_Pha_loc +
									 x * pWave->kx_piece +
									 y * pWave->ky_piece),
					   1);
		Sm = arma::sum(H_Mag_loc % pWave->amplitudes_piece(ind_piece(0)) %
						   arma::sin(pWave->phases_piece(ind_piece(0)) +
									 H_Pha_loc +
									 x * pWave->kx_piece +
									 y * pWave->ky_piece),
					   1);

		// Compute amplitude and phase of the excitation force for the current wave piece
		Am = arma::sqrt(arma::pow(Cm, 2) + arma::pow(Sm, 2));
		PHIm = arma::atan2(Sm, Cm);

		// Compute excitation force for the current wave piece
		F_piece = arma::as_scalar(arma::sum(Am % arma::cos(t * pWave->ang_freqs_piece - PHIm), 0));

		if (flag_gap)
		{
			// Multiply the transfer functions with the wave amplitudes for the current wave piece
			Cm = arma::sum(H_Mag_loc % pWave->amplitudes_piece(ind_piece(0) - 1) %
							   arma::cos(pWave->phases_piece(ind_piece(0) - 1) +
										 H_Pha_loc +
										 x * pWave->kx_piece +
										 y * pWave->ky_piece),
						   1);
			Sm = arma::sum(H_Mag_loc % pWave->amplitudes_piece(ind_piece(0) - 1) %
							   arma::sin(pWave->phases_piece(ind_piece(0) - 1) +
										 H_Pha_loc +
										 x * pWave->kx_piece +
										 y * pWave->ky_piece),
						   1);

			// Compute amplitude and phase of the excitation force for the current wave piece
			Am = arma::sqrt(arma::pow(Cm, 2) + arma::pow(Sm, 2));
			PHIm = arma::atan2(Sm, Cm);

			// Compute excitation force for the current wave piece
			F_gap = arma::as_scalar(arma::sum(Am % arma::cos(t * pWave->ang_freqs_piece - PHIm), 0));

			// Store the excitation force mixing linearly the current and previous piece
			// TODO: Mix pieces with a qubic function instead of linear
			Fe(ii, 0) = F_piece * (t - pWave->time_end(ind_piece(0) - 1)) / pWave->time_gap +
						F_gap * (pWave->time_ini(ind_piece(0)) - t) / pWave->time_gap;
		}
		else
		{
			// Store the excitation force of the current piece
			Fe(ii, 0) = F_piece;
		}
	}

	// Rotate the excitation force to the fixed frame
	double Fx = arma::as_scalar(Fe(0, 0));
	double Fy = arma::as_scalar(Fe(1, 0));
	Fe(0, 0) = Fx * cos(yaw) - Fy * sin(yaw);
	Fe(1, 0) = Fx * sin(yaw) + Fy * cos(yaw);

	// Store the excitation force in the body variable for plotting
	pBodies[idBody]->excitationForces_1 = Fe;

	// Return the excitation force
	return Fe;
}

arma::mat HydroDatabase::ComputeSecondWaveExcForce(double t)
{
	// Initialize excitation force vector
	arma::mat Fe = arma::zeros(activeDofs, 1);

	// Load wave pointer as local variable
	Wave *pWave = pSim->pWave;

	// Extract position and heading of the body
	double x = pBodies[idBody]->pos(0, 0);
	double y = pBodies[idBody]->pos(1, 0);
	double yaw = pBodies[idBody]->pos(5, 0);

	// TODO: Remove this when the precomputed second order excitation forces are implemented
	if (pBodies[idBody]->secondOrderExcitationFlag == 1)
	{
		x = 0;
		y = 0;
		yaw = 0;
	}

	// Declare local variables
	arma::cube temp1;
	arma::mat temp2, temp3;
	double F_piece, F_gap;

	// Interpolate transfer functions to the body heading
	arma::cube ***HDif = new arma::cube **[2];
	arma::cube ***HSum = new arma::cube **[2];
	for (int ii = 0; ii < 2; ii++)
	{
		HDif[ii] = new arma::cube *[activeDofs];
		HSum[ii] = new arma::cube *[activeDofs];
		for (int jj = 0; jj < activeDofs; jj++)
		{
			HDif[ii][jj] = new arma::cube;
			HSum[ii][jj] = new arma::cube;

			temp1 = *QtfDiff_w[ii][jj];
			temp1 = interp1((*pHeadings) + yaw, permute(temp1, 312), pWave->headings_1D);
			*HDif[ii][jj] = permute(temp1, 231);

			temp1 = *QtfSum_w[ii][jj];
			temp1 = interp1((*pHeadings) + yaw, permute(temp1, 312), pWave->headings_1D);
			*HSum[ii][jj] = permute(temp1, 231);
		}
	}

	// Extract the index of the required wave piece
	arma::uvec ind_piece = arma::find((pWave->time_ini - pWave->time_gap) < t, 1, "last");
	// Check if current time is in a gap
	bool flag_gap = ((ind_piece(0) > 0) && (t < pWave->time_ini(ind_piece(0))));

	// Loop over all active degrees of freedom
	// TODO: This only works for 6 DOFs, change it.
	for (int ii = 0; ii < activeDofs; ii++)
	{
		// Initialize the excitation force for the current DOF
		F_piece = 0.0;

		// Accumulate the real part of the difference QTF for the current piece and DOF
		temp2 = (*HDif[0][ii]).slice(0);
		temp3 = temp2 % ampP(ind_piece(0)) * arma::cos(wD * t + phD(ind_piece(0)) + kxD * x + kyD * y);
		F_piece = F_piece + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

		// Accumulate the imaginary part of the difference QTF for the current piece and DOF
		temp2 = (*HDif[1][ii]).slice(0);
		temp3 = temp2 % ampP(ind_piece(0)) * arma::sin(wD * t + phD(ind_piece(0)) + kxD * x + kyD * y);
		F_piece = F_piece + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

		// Accumulate the real part of the sum QTF for the current piece and DOF
		temp2 = (*HSum[0][ii]).slice(0);
		temp3 = temp2 % ampP(ind_piece(0)) * arma::cos(wS * t + phS(ind_piece(0)) + kxS * x + kyS * y);
		F_piece = F_piece + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

		// Accumulate the imaginary part of the sum QTF for the current piece and DOF
		temp2 = (*HSum[1][ii]).slice(0);
		temp3 = temp2 % ampP(ind_piece(0)) * arma::sin(wS * t + phS(ind_piece(0)) + kxS * x + kyS * y);
		F_piece = F_piece + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

		if (flag_gap)
		{
			// Compute excitation force for the current wave piece
			F_gap = 0.0;

			// Accumulate the real part of the difference QTF for the current piece and DOF
			temp2 = (*HDif[0][ii]).slice(0);
			temp3 = temp2 % ampP(ind_piece(0) - 1) * arma::cos(wD * t + phD(ind_piece(0) - 1) + kxD * x + kyD * y);
			F_gap = F_gap + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

			// Accumulate the imaginary part of the difference QTF for the current piece and DOF
			temp2 = (*HDif[1][ii]).slice(0);
			temp3 = temp2 % ampP(ind_piece(0) - 1) * arma::sin(wD * t + phD(ind_piece(0) - 1) + kxD * x + kyD * y);
			F_gap = F_gap + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

			// Accumulate the real part of the sum QTF for the current piece and DOF
			temp2 = (*HSum[0][ii]).slice(0);
			temp3 = temp2 % ampP(ind_piece(0) - 1) * arma::cos(wS * t + phS(ind_piece(0) - 1) + kxS * x + kyS * y);
			F_gap = F_gap + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

			// Accumulate the imaginary part of the sum QTF for the current piece and DOF
			temp2 = (*HSum[1][ii]).slice(0);
			temp3 = temp2 % ampP(ind_piece(0) - 1) * arma::sin(wS * t + phS(ind_piece(0) - 1) + kxS * x + kyS * y);
			F_gap = F_gap + 0.5 * arma::as_scalar(arma::sum(arma::sum(temp3, 1), 0));

			// Store the excitation force mixing linearly the current and previous piece
			// TODO: Mix pieces with a qubic function instead of linear
			Fe(ii, 0) = F_piece * (t - pWave->time_end(ind_piece(0) - 1)) / pWave->time_gap +
						F_gap * (pWave->time_ini(ind_piece(0)) - t) / pWave->time_gap;
		}
		else
		{
			// Store the excitation force of the current piece
			Fe(ii, 0) = F_piece;
		}
	}

	// Rotate the excitation force to the fixed frame
	double Fx = arma::as_scalar(Fe(0, 0));
	double Fy = arma::as_scalar(Fe(1, 0));
	Fe(0, 0) = Fx * cos(yaw) - Fy * sin(yaw);
	Fe(1, 0) = Fx * sin(yaw) + Fy * cos(yaw);

	// Store the excitation force in the body variable for plotting
	pBodies[idBody]->excitationForces_2 = Fe;

	// Return the excitation force
	return Fe;
}

void HydroDatabase::SetUp(void)
{
	// Load the wave pointer as a local variable
	Wave *pWave = pSim->pWave;

	// Convert from real/imag to mag/pha the first order transfer functions
	arma::cube WE_Real = (*pWaveExcitingMag) % arma::cos((*pWaveExcitingPha));
	arma::cube WE_Imag = (*pWaveExcitingMag) % arma::sin((*pWaveExcitingPha));

	// Check that the transfer functions frequencies cover the wave frequencies, return a warning if not
	if (((*pFrequencies).min() > (pWave->freqs_piece).min()) || ((*pFrequencies).max() < (pWave->freqs_piece).max()))
	{
		std::cout << "WARNING: The frequencies provided in the hydrodinamic data base do not cover properly the wave!" << std::endl;
	}

	// Interpolate the first order transfer functions to the wave frequencies
	WE_Real_w = interp1(*pFrequencies, permute(WE_Real, 231), pWave->freqs_piece);
	WE_Imag_w = interp1(*pFrequencies, permute(WE_Imag, 231), pWave->freqs_piece);
	WE_Real_w = permute(WE_Real_w, 213);
	WE_Imag_w = permute(WE_Imag_w, 213);

	// Preprocess for QTFs
	if (pBodies[idBody]->secondOrderExcitationFlag > 0)
	{
		// Interpolate the second order transfer functions to the wave frequencies
		arma::cube temp;
		QtfDiff_w = new arma::cube **[2];
		QtfSum_w = new arma::cube **[2];
		for (int ii = 0; ii < 2; ii++)
		{
			QtfDiff_w[ii] = new arma::cube *[activeDofs];
			QtfSum_w[ii] = new arma::cube *[activeDofs];
			for (int jj = 0; jj < activeDofs; jj++)
			{
				QtfDiff_w[ii][jj] = new arma::cube;
				QtfSum_w[ii][jj] = new arma::cube;

				temp = *pQtfDiff[ii][jj];
				*QtfDiff_w[ii][jj] = interp2(*pFrequencies, *pFrequencies, temp, pWave->freqs_piece, pWave->freqs_piece);

				temp = *pQtfSum[ii][jj];
				*QtfSum_w[ii][jj] = interp2(*pFrequencies, *pFrequencies, temp, pWave->freqs_piece, pWave->freqs_piece);
			}
		}

		// Preprocess the matrices required for time domain QTF forces computation
		ampP = arma::field<arma::mat>(pWave->num_pieces);
		phS = arma::field<arma::mat>(pWave->num_pieces);
		phD = arma::field<arma::mat>(pWave->num_pieces);

		wS = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
		kxS = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
		kyS = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
		wD = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
		kxD = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
		kyD = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);

		for (int ii = 0; ii < pWave->num_comps_piece; ii++)
		{
			for (int jj = 0; jj < pWave->num_comps_piece; jj++)
			{
				wS(ii, jj) = pWave->ang_freqs_piece(ii) + pWave->ang_freqs_piece(jj);
				kxS(ii, jj) = pWave->kx_1D_piece(ii) + pWave->kx_1D_piece(jj);
				kyS(ii, jj) = pWave->ky_1D_piece(ii) + pWave->ky_1D_piece(jj);
				wD(ii, jj) = pWave->ang_freqs_piece(ii) - pWave->ang_freqs_piece(jj);
				kxD(ii, jj) = pWave->kx_1D_piece(ii) - pWave->kx_1D_piece(jj);
				kyD(ii, jj) = pWave->ky_1D_piece(ii) - pWave->ky_1D_piece(jj);
			}
		}

		for (int kk = 0; kk < pWave->num_pieces; kk++)
		{
			ampP(kk) = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
			phS(kk) = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
			phD(kk) = arma::zeros(pWave->num_comps_piece, pWave->num_comps_piece);
			for (int ii = 0; ii < pWave->num_comps_piece; ii++)
			{
				for (int jj = 0; jj < pWave->num_comps_piece; jj++)
				{
					ampP(kk)(ii, jj) = pWave->amplitudes_1D_piece(kk)(ii) * pWave->amplitudes_1D_piece(kk)(jj);
					phS(kk)(ii, jj) = pWave->phases_1D_piece(kk)(ii) + pWave->phases_1D_piece(kk)(jj);
					phD(kk)(ii, jj) = pWave->phases_1D_piece(kk)(ii) - pWave->phases_1D_piece(kk)(jj);
				}
			}
		}

		// Compute mean drift force
		// TODO: Review this!
		arma::cube temp_mD = interp2(*pFrequencies, *pHeadings, permute(*pMeanDrift, 231), pWave->freqs_piece, pWave->headings_piece);
		for (int ii = 0; ii < activeDofs; ii++)
		{
			temp_mD.slice(ii) = temp_mD.slice(ii) % pWave->amplitudes_piece(0) % pWave->amplitudes_piece(0);
		}
		temp_mD = permute(temp_mD, 312);
		F_meanDrift = arma::sum(arma::sum(temp_mD, 1), 2);
		if (pBodies[idBody]->secondOrderExcitationFlag == 3 || pBodies[idBody]->secondOrderExcitationFlag == 4)
		{
			pBodies[idBody]->excitationForces_2 = F_meanDrift;
			std::cout << "Computed mean drift: \n"
					  << F_meanDrift << "\n";
		}
	}

	// std::cout << "Calculate hyrodynamic forces at time 0.0...\n";
	pBodies[idBody]->Fb = CalculateHydrodynamicForces(0.0);

	if (pBodies[idBody]->firstOrderExcitationFlag == 1)
	{
		// TODO: Precompute first order forces
		if (pSim->simulationTime <= 0)
		{
			std::cout << "WARNING: Precomputed first order forces not implemented yet. \n"
						 "         Using limited instant position version instead. \n"
					  << std::endl;
		}
	}
	if (pBodies[idBody]->secondOrderExcitationFlag == 1)
	{
		// TODO: Precompute second order forces
		if (pSim->simulationTime <= 0)
		{
			std::cout << "WARNING: Precomputed second order forces not implemented yet. \n"
						 "         Using limited instant position version instead. \n"
					  << std::endl;
		}
	}
}

arma::mat HydroDatabase::ComputeMeanDrift(void)
{

	// // Compute mean drift force
	// Wave* pWave = pSim->pWave;
	// double yaw = pBodies[idBody]->pos(5,0);
	// // HARCODEO !!!!!!!!!!!! ------------------------------------------------------------ Implementar opcion de mean drift con y sin instant position
	// if(pBodies[idBody]->secondOrderExcitationFlag<5)
	// {
	// 	yaw = 0;
	// }
	// arma::cube temp_mD = interp2(*pFrequencies, *pHeadings+yaw, permute(*pMeanDrift,231), pWave->freqs, pWave->headings);
	// for (int ii=0; ii<activeDofs; ii++)
	// {
	// 	temp_mD.slice(ii) = temp_mD.slice(ii)%pWave->amplitudes%pWave->amplitudes;
	// }
	// temp_mD = permute(temp_mD,312);
	// arma::mat F = arma::sum(arma::sum(temp_mD,1),2);
	// double Fx = arma::as_scalar(F(0,0)); double Fy = arma::as_scalar(F(1,0));
	// F(0,0) = Fx*cos(yaw)-Fy*sin(yaw); F(1,0) = Fx*sin(yaw)+Fy*cos(yaw);
	// pBodies[idBody]->excitationForces_2 = F;

	arma::mat F = F_meanDrift;

	return F;
}

void HydroDatabase::Print()
{
	std::cout << "Number of bodies associated: " << numBodies << std::endl;
	std::cout << "Number of frequencies: " << numFrequencies << std::endl;
	std::cout << "Number of headings: " << numHeadings << std::endl;
}

void HydroDatabase::InterpolateHydro(HydroDatabase *pHydro1, HydroDatabase *pHydro2, double interpCoef)
{
	if (pHydro1->numBodies != pHydro1->numBodies)
	{
		std::stringstream ss;
		ss << "ERROR: The number of bodies is not the same in interpolated databases \n";
		throw ValueError(ss.str());
	}

	if (pHydro1->numFrequencies != pHydro1->numFrequencies)
	{
		std::stringstream ss;
		ss << "ERROR: The number of frequencies is not the same in interpolated databases \n";
		throw ValueError(ss.str());
	}

	if (pHydro1->numHeadings != pHydro1->numHeadings)
	{
		std::stringstream ss;
		ss << "ERROR: The number of headings is not the same in interpolated databases \n";
		throw ValueError(ss.str());
	}

	// std::cout << "            Interpolate pHydrostaticStiffness..." << std::endl;
	(*pHydrostaticStiffness) = *(pHydro1->pHydrostaticStiffness) * (1 - interpCoef) + *(pHydro2->pHydrostaticStiffness) * interpCoef;
	// std::cout << "            Interpolate pStructuralMass..." << std::endl;
	(*pStructuralMass) = *(pHydro1->pStructuralMass) * (1 - interpCoef) + *(pHydro2->pStructuralMass) * interpCoef;
	// std::cout << "            Interpolate pTotalMass..." << std::endl;
	(*pTotalMass) = *(pHydro1->pTotalMass) * (1 - interpCoef) + *(pHydro2->pTotalMass) * interpCoef;

	// std::cout << "            Interpolate pWaveExcitingMag..." << std::endl;
	arma::cube WaveExcitingReal_1 = *(pHydro1->pWaveExcitingMag) % arma::cos(*(pHydro1->pWaveExcitingPha));
	arma::cube WaveExcitingImag_1 = *(pHydro1->pWaveExcitingMag) % arma::sin(*(pHydro1->pWaveExcitingPha));
	arma::cube WaveExcitingReal_2 = *(pHydro2->pWaveExcitingMag) % arma::cos(*(pHydro2->pWaveExcitingPha));
	arma::cube WaveExcitingImag_2 = *(pHydro2->pWaveExcitingMag) % arma::sin(*(pHydro2->pWaveExcitingPha));
	arma::cx_cube WaveExcitingCx_1 = arma::cx_cube(WaveExcitingReal_1, WaveExcitingImag_1);
	arma::cx_cube WaveExcitingCx_2 = arma::cx_cube(WaveExcitingReal_2, WaveExcitingImag_2);
	arma::cx_cube WaveExcitingCx = WaveExcitingCx_1 * (1 - interpCoef) + WaveExcitingCx_2 * interpCoef;
	(*pWaveExcitingMag) = arma::abs(WaveExcitingCx);
	(*pWaveExcitingPha) = arma::arg(WaveExcitingCx);

	// std::cout << "            Interpolate added mass and damping..." << std::endl;
	for (int ii = 0; ii < numBodies; ii++)
	{
		*pAddedMass[ii] = *(pHydro1->pAddedMass[ii]) * (1 - interpCoef) + *(pHydro2->pAddedMass[ii]) * interpCoef;
		*pAddedMassHf[ii] = *(pHydro1->pAddedMassHf[ii]) * (1 - interpCoef) + *(pHydro2->pAddedMassHf[ii]) * interpCoef;
		*pAddedMassLf[ii] = *(pHydro1->pAddedMassLf[ii]) * (1 - interpCoef) + *(pHydro2->pAddedMassLf[ii]) * interpCoef;
		*pDampingRadiation[ii] = *(pHydro1->pDampingRadiation[ii]) * (1 - interpCoef) + *(pHydro2->pDampingRadiation[ii]) * interpCoef;
		*pDampingRadiationLf[ii] = *(pHydro1->pDampingRadiationLf[ii]) * (1 - interpCoef) + *(pHydro2->pDampingRadiationLf[ii]) * interpCoef;
		*pIRF[ii] = *(pHydro1->pIRF[ii]) * (1 - interpCoef) + *(pHydro2->pIRF[ii]) * interpCoef;
	}

	// std::cout << "            Interpolate pQtfDiff..." << std::endl;
	// Esto esta un poco feo, deberíamos tener un flag que nos diga si tenemos QTFs o no en la HDB y tirar de eso
	if (pBodies[idBody]->secondOrderExcitationFlag > 0 && pBodies[idBody]->secondOrderExcitationFlag < 4)
	{
		for (int ii = 0; ii < 2; ii++)
		{
			for (int jj = 0; jj < activeDofs; jj++)
			{
				*pQtfDiff[ii][jj] = *(pHydro1->pQtfDiff[ii][jj]) * (1 - interpCoef) + *(pHydro2->pQtfDiff[ii][jj]) * interpCoef;
				*pQtfSum[ii][jj] = *(pHydro1->pQtfSum[ii][jj]) * (1 - interpCoef) + *(pHydro2->pQtfSum[ii][jj]) * interpCoef;
			}
		}
	}

	// std::cout << "            Interpolate SetUp();" << std::endl;
	SetUp(); // A esta quizas habria que llamarla desde sinking e interpolar con las variables postprocesadas de setup
}

void HydroDatabase::UpdateHydroStiffness(arma::mat newHydrostaticStiffness)
{
	*pHydrostaticStiffness = newHydrostaticStiffness;
}

arma::mat HydroDatabase::CalculateHydrostaticPressure(double t)
{
	// std::cout << "--> Calculating Hydrostatic Pressure" << std::endl;

	arma::mat pressure;
	arma::mat z = pBodies[idBody]->pNLHSMesh->nodes.col(2);

	if (pBodies[idBody]->flag_hidrostatics == 1)
	{

		// Non-linear hydrostatic forces without wave
		pressure = -z * pSim->gravity * pSim->waterDensity;
	}
	else if (pBodies[idBody]->flag_hidrostatics == 2)
	{

		// Non-linear hydrostatic forces with wave
		arma::mat x = pBodies[idBody]->pNLHSMesh->nodes.col(0);
		arma::mat y = pBodies[idBody]->pNLHSMesh->nodes.col(1);
		pressure = pSim->pWave->GetPressure(t, x, y, z);
	}

	return pressure;
}
