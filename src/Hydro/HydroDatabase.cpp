
#include <armadillo>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <chrono>
#include "hdf5.h"
#include "H5Cpp.h"
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

	double rampa = std::min(1.0, time / 50.0); // TODO: the ramp time should not be hardcoded!
	double yaw = pBodies[idBody]->pos(5, 0);

	// std::cout << "--> Computing radiation forces..." << std::endl;
	if (pBodies[idBody]->radiationFlag == 1)
	{
		F = F + ComputeRadiationForces();
	}

	// std::cout << "--> Computing first order diffraction forces..." << std::endl;
	if (pBodies[idBody]->firstOrderExcitationFlag == 1)
	{
		F = F + ComputeFirstWaveExcForce(time) * rampa;
		// std::cout << "WARNING: Precomputed first order forces not implemented yet. \n" << std::endl;
	}
	if (pBodies[idBody]->firstOrderExcitationFlag == 2)
	{
		F = F + ComputeFirstWaveExcForce(time) * rampa;
	}

	// std::cout << "--> Computing second order diffraction forces..." << std::endl;
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

	// std::cout << "--> Computing viscous drag forces..." << std::endl;
	// TODO: Consider moving all fast computing components to the Hydrostatic forces class
	arma::mat vv = pBodies[idBody]->vel;
	F = F - pBodies[idBody]->B_visc % vv - pBodies[idBody]->B_visc2 % vv % arma::abs(vv);

	// std::cout << "--> Computing wind and current forces..." << std::endl;
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
	if (pBodies[idBody]->flag_hydrostatics == 0)
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
	else if (pBodies[idBody]->flag_hydrostatics > 0)
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

void HydroDatabase::ComputeAsymptoticAddedMass(std::string HDBname)
{
	// Get the index for the lowest frequency
	arma::uvec id_min_freq = arma::find((*pFrequencies) == (*pFrequencies).min());

	// Allocate the asymptotic added mass matrix
	pAddedMassLf = new arma::mat *[numBodies];
	pAddedMassHf = new arma::mat *[numBodies];
	for (int ib = 0; ib < numBodies; ib++)
	{
		pAddedMassLf[ib] = new arma::mat(6, 6, arma::fill::zeros);
		pAddedMassHf[ib] = new arma::mat(6, 6, arma::fill::zeros);
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				// Calculate the asymptotic low frequency added mass as an approximation
				// to the added mass at the lowest frequency
				(*pAddedMassLf[ib])(i, j) = (*pAddedMass[ib])(i, j, id_min_freq(0));
				// Calculate the asymptotic high frequency added mass as the integral of the IRF
				(*pAddedMassHf[ib])(i, j) += trapzi(IRFTime, (*pIRF[ib]).subcube(0, i, j, numPointsIRF - 1, i, j));
				// TODO: REVIEW THIS!!!!!!!!!!!!!
			}
		}
	}
}

void HydroDatabase::ComputeTotalMass(void)
{
	// Create total mass matrix and fill with structural data
	std::cout << "  Creating total mass matrix...\n";
	pTotalMass = new arma::mat(6, 6 * numBodies, arma::fill::zeros);
	(*pTotalMass)(arma::span(0, 5), arma::span(6 * (pBodies[idBody]->hydroDatabaseIndex), 6 * (pBodies[idBody]->hydroDatabaseIndex + 1) - 1)) = (*pStructuralMass);
	for (int ii = 0; ii < numBodies; ii++)
	{
		(*pTotalMass)(arma::span(0, 5), arma::span(6 * ii, 6 * (ii + 1) - 1)) += (*pAddedMassHf[ii]);
	}
	(*pTotalMass)(arma::span(0, 5), arma::span(6 * (pBodies[idBody]->hydroDatabaseIndex), 6 * (pBodies[idBody]->hydroDatabaseIndex + 1) - 1)) +=
		(*pAddedMassHf[pBodies[idBody]->hydroDatabaseIndex]) % (arma::diagmat(pBodies[idBody]->A_visc));
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
	// Check the end of the file name to determine the type of database
	if (filePath.find(".ehydb") != std::string::npos)
	{
		hydroDatabaseFlag = 0;
	}
	else if (filePath.find(".hydb.h5") != std::string::npos)
	{
		hydroDatabaseFlag = 1;
	}
	else
	{
		std::cout << "ERROR: File name does not have a valid extension. \n";
		std::cout << "       Valid extensions are: .hydb.h5 or .ehydb \n";
		throw std::exception();
	}

	if (hydroDatabaseFlag == 0)
	{
		std::cout << "  Loading hydrodynamic data from file: " << filePath << "...\n";
		LoadHydroDataEHYDB(filePath);
	}
	else if (hydroDatabaseFlag == 1)
	{
		std::cout << "  Loading hydrodynamic data from file: " << filePath << "...\n";
		LoadHydroDataH5(filePath);
	}

	// Load the structural mass into the body
	pBodies[idBody]->structuralMass = (*pStructuralMass)(0, 0);

	// Generate starting pos time matrix
	pTimeStartPos = new arma::cube(numBodies, 6, 6, arma::fill::zeros);

	// Assemble the HDF name in order to save the IRF data
	std::string HDBname = filePath.substr(filePath.find_last_of("/") + 1);
	if (hydroDatabaseFlag == 0)
	{
		HDBname = HDBname.substr(0, HDBname.length() - 6);
	}
	else if (hydroDatabaseFlag == 1)
	{
		HDBname = HDBname.substr(0, HDBname.length() - 8);
	}

	// Compute IRF function
	std::cout << "  Computing IRF ...\n";
	this->ComputeIRF(HDBname);
	std::cout << "  ... computing IRF done!\n";

	// Compute Asymptotic added mass if it was not loaded from the file
	if (hydroDatabaseFlag == 1)
	{
		std::cout << "  Computing asymptotic added mass...\n";
		this->ComputeAsymptoticAddedMass(HDBname);
		std::cout << "  ... computing asymptotic added mass done!\n";
	}

	// Compute total mass matrix
	std::cout << "  Computing total mass matrix...\n";
	this->ComputeTotalMass();
	std::cout << "  ... computing total mass matrix done!\n";

	// Load Morison forces data
	std::cout << "  Reading Morison forces data ...\n";
	pMor = new Morison(numBodies, pSim);
	pMor->ReadMorisonData();
}

void HydroDatabase::LoadHydroDataEHYDB(std::string filePath)
{
	// Read number of bodies
	std::cout << "Loading hydrodynamic data from file: " << filePath << "...\n";
	std::cout << "  Reading number of bodies...\n";
	arma::mat num_bodies_mat;
	num_bodies_mat.load(arma::hdf5_name(filePath, "num_bodies"));
	numBodies = num_bodies_mat(0);

	// Check if the number of bodies is correct
	if ((numBodies != pSim->numBodies) && (numBodies > 1))
	{
		std::stringstream ss;
		ss << " The number of bodies in datosBodies.dat does not match the number of bodies in the multibody hydrodatabase." << std::endl;
		throw IOError(ss.str());
	}

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
	pFrequencies->load(arma::hdf5_name(filePath, frequencies_fn.str(), arma::hdf5_opts::trans));
	numFrequencies = pFrequencies->n_rows;

	// Read headings (radians)
	std::cout << "  Reading headings...\n";
	std::stringstream headings_fn;
	pHeadings = new arma::vec;
	headings_fn << "body_" << this->GetId() << "/headings";
	pHeadings->load(arma::hdf5_name(filePath, headings_fn.str(), arma::hdf5_opts::trans));
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
	}

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
	if (pBodies[idBody]->flag_hydrostatics == 2)
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
	if (pBodies[idBody]->flag_hydrostatics == 2)
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

		std::cout << "  Computing mean drift coefficients...\n";
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
}

void HydroDatabase::LoadHydroDataH5(std::string filePath)
{
	// Declare local variables
	int num_values;

	// Open the HDF5 file
	H5::H5File file(filePath, H5F_ACC_RDONLY);

	// TODO: Load cog
	cog = arma::zeros(1, 3);

	// Get the number of frequencies from the "frequencies" dataset which is a 1D array
	H5::DataSet frequenciesDataset = file.openDataSet("/frequencies");
	H5::DataSpace frequenciesSpace = frequenciesDataset.getSpace();
	hsize_t numFrequencies_t;
	frequenciesSpace.getSimpleExtentDims(&numFrequencies_t, NULL);
	numFrequencies = numFrequencies_t;
	std::cout << "Number of frequencies: " << numFrequencies << std::endl;
	// Allocate memory for the frequencies buffer
	double *buffer_frequencies = new double[numFrequencies];
	// Read the frequencies data
	frequenciesDataset.read(buffer_frequencies, H5::PredType::NATIVE_DOUBLE);
	// Close the dataset
	frequenciesDataset.close();
	// Allocate memory for the frequencies array
	pFrequencies = new arma::vec(numFrequencies);
	// Copy the data from the buffer to the frequencies array
	for (int i = 0; i < numFrequencies; i++)
	{
		(*pFrequencies)(i) = buffer_frequencies[i];
	}
	// Free the buffer memory
	delete[] buffer_frequencies;

	// Get the number of headings from the "headings" dataset which is a 1D array
	H5::DataSet headingsDataset = file.openDataSet("/headings");
	H5::DataSpace headingsSpace = headingsDataset.getSpace();
	hsize_t numHeadings_t;
	headingsSpace.getSimpleExtentDims(&numHeadings_t, NULL);
	numHeadings = numHeadings_t;
	std::cout << "Number of headings: " << numHeadings << std::endl;
	// Allocate memory for the headings buffer
	double *buffer_headings = new double[numHeadings];
	// Read the headings data
	headingsDataset.read(buffer_headings, H5::PredType::NATIVE_DOUBLE);
	// Close the dataset
	headingsDataset.close();
	// Allocate memory for the headings array
	pHeadings = new arma::vec(numHeadings);
	// Copy the data from the buffer to the headings array
	for (int i = 0; i < numHeadings; i++)
	{
		(*pHeadings)(i) = buffer_headings[i];
	}
	// Free the buffer memory
	delete[] buffer_headings;

	// Get the number of bodies as the number of groups in the "mesh" group
	H5::Group meshGroup = file.openGroup("/mesh");
	numBodies = meshGroup.getNumObjs();
	if ((numBodies != pSim->numBodies) && (numBodies > 1))
	{
		std::stringstream ss;
		ss << " The number of bodies in datosBodies.dat does not match the number of bodies in the multibody hydrodatabase." << std::endl;
		throw IOError(ss.str());
	}
	std::cout << "Number of bodies: " << numBodies << std::endl;
	// Close the group
	meshGroup.close();

	// Read the "added_mass" and "damping_rad" datasets which are 5D arrays [numBodies, numBodies, numFrequencies, 6, 6]
	std::cout << "Reading added mass and damping radiation data...\n";
	H5::DataSet addedMassDataset = file.openDataSet("/added_mass");
	H5::DataSpace addedMassSpace = addedMassDataset.getSpace();
	hsize_t dims_am[5];
	addedMassSpace.getSimpleExtentDims(dims_am, NULL);
	// Open the "damping_rad" dataset
	H5::DataSet dampingRadDataset = file.openDataSet("/damping_rad");
	H5::DataSpace dampingRadSpace = dampingRadDataset.getSpace();
	hsize_t dims_dr[5];
	dampingRadSpace.getSimpleExtentDims(dims_dr, NULL);
	// Allocate memory for the added mass and damping radiation buffers [numBodies, numBodies, numFrequencies, 6, 6]
	num_values = numBodies * numBodies * numFrequencies * 6 * 6;
	double *buffer_added_mass = new double[num_values];
	double *buffer_damping_rad = new double[num_values];
	// Read the added mass and damping radiation data
	addedMassDataset.read(buffer_added_mass, H5::PredType::NATIVE_DOUBLE, addedMassSpace, addedMassSpace);
	dampingRadDataset.read(buffer_damping_rad, H5::PredType::NATIVE_DOUBLE, dampingRadSpace, dampingRadSpace);
	// Close the datasets
	addedMassDataset.close();
	dampingRadDataset.close();
	// Allocate memory for the added mass and damping radiation matrices
	pAddedMass = new arma::cube *[numBodies];
	pDampingRadiation = new arma::cube *[numBodies];
	for (int i = 0; i < numBodies; i++)
	{
		pAddedMass[i] = new arma::cube;
		pDampingRadiation[i] = new arma::cube;
		pAddedMass[i]->set_size(6, 6, numFrequencies);
		pDampingRadiation[i]->set_size(6, 6, numFrequencies);
		for (int j = 0; j < numFrequencies; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				for (int l = 0; l < 6; l++)
				{
					// Calculate the index for the 1D buffer [numBodies, numBodies, numFrequencies, 6, 6]
					int index = i * numBodies * numFrequencies * 6 * 6 + idBody * numFrequencies * 6 * 6 + j * 6 * 6 + k * 6 + l;
					// Copy data from buffer to arrays
					(*pAddedMass[i])(k, l, j) = buffer_added_mass[index];
					(*pDampingRadiation[i])(k, l, j) = buffer_damping_rad[index];
				}
			}
		}
	}
	// Free the buffer memory
	delete[] buffer_added_mass;
	delete[] buffer_damping_rad;

	// Read the "diffraction_force_mag", "diffraction_force_pha"
	// "froude_krylov_force_mag", "froude_krylov_force_pha"
	// "wave_exciting_mag" and "wave_exciting_pha"
	// datasets which are 4D arrays [numHeadings, numBodies, numFrequencies, 6]
	std::cout << "Reading diffraction force magnitude and phase data...\n";
	// Open "diffraction_force_mag"
	H5::DataSet diffractionForceMagDataset = file.openDataSet("/diffraction_force_mag");
	H5::DataSpace diffractionForceMagSpace = diffractionForceMagDataset.getSpace();
	hsize_t dims_dfm[4];
	diffractionForceMagSpace.getSimpleExtentDims(dims_dfm, NULL);
	// Open "diffraction_force_pha"
	H5::DataSet diffractionForcePhaDataset = file.openDataSet("/diffraction_force_pha");
	H5::DataSpace diffractionForcePhaSpace = diffractionForcePhaDataset.getSpace();
	hsize_t dims_dfp[4];
	diffractionForcePhaSpace.getSimpleExtentDims(dims_dfp, NULL);
	// Open "froude_krylov_force_mag"
	H5::DataSet froudeKrylovForceMagDataset = file.openDataSet("/froude_krylov_force_mag");
	H5::DataSpace froudeKrylovForceMagSpace = froudeKrylovForceMagDataset.getSpace();
	hsize_t dims_fkm[4];
	froudeKrylovForceMagSpace.getSimpleExtentDims(dims_fkm, NULL);
	// Open "froude_krylov_force_pha"
	H5::DataSet froudeKrylovForcePhaDataset = file.openDataSet("/froude_krylov_force_pha");
	H5::DataSpace froudeKrylovForcePhaSpace = froudeKrylovForcePhaDataset.getSpace();
	hsize_t dims_fkp[4];
	froudeKrylovForcePhaSpace.getSimpleExtentDims(dims_fkp, NULL);
	// Open "wave_exciting_mag"
	H5::DataSet waveExcitingMagDataset = file.openDataSet("/wave_exciting_mag");
	H5::DataSpace waveExcitingMagSpace = waveExcitingMagDataset.getSpace();
	hsize_t dims_wem[4];
	waveExcitingMagSpace.getSimpleExtentDims(dims_wem, NULL);
	// Open "wave_exciting_pha"
	H5::DataSet waveExcitingPhaDataset = file.openDataSet("/wave_exciting_pha");
	H5::DataSpace waveExcitingPhaSpace = waveExcitingPhaDataset.getSpace();
	hsize_t dims_wep[4];
	waveExcitingPhaSpace.getSimpleExtentDims(dims_wep, NULL);
	// Allocate memory for the diffraction force magnitude buffer [numHeadings, numBodies, numFrequencies, 6]
	num_values = numHeadings * numBodies * numFrequencies * 6;
	double *buffer_diffraction_force_mag = new double[num_values];
	double *buffer_diffraction_force_pha = new double[num_values];
	double *buffer_froude_krylov_force_mag = new double[num_values];
	double *buffer_froude_krylov_force_pha = new double[num_values];
	double *buffer_wave_exciting_mag = new double[num_values];
	double *buffer_wave_exciting_pha = new double[num_values];
	// Read the excitation force data
	diffractionForceMagDataset.read(buffer_diffraction_force_mag, H5::PredType::NATIVE_DOUBLE, diffractionForceMagSpace, diffractionForceMagSpace);
	diffractionForcePhaDataset.read(buffer_diffraction_force_pha, H5::PredType::NATIVE_DOUBLE, diffractionForcePhaSpace, diffractionForcePhaSpace);
	froudeKrylovForceMagDataset.read(buffer_froude_krylov_force_mag, H5::PredType::NATIVE_DOUBLE, froudeKrylovForceMagSpace, froudeKrylovForceMagSpace);
	froudeKrylovForcePhaDataset.read(buffer_froude_krylov_force_pha, H5::PredType::NATIVE_DOUBLE, froudeKrylovForcePhaSpace, froudeKrylovForcePhaSpace);
	waveExcitingMagDataset.read(buffer_wave_exciting_mag, H5::PredType::NATIVE_DOUBLE, waveExcitingMagSpace, waveExcitingMagSpace);
	waveExcitingPhaDataset.read(buffer_wave_exciting_pha, H5::PredType::NATIVE_DOUBLE, waveExcitingPhaSpace, waveExcitingPhaSpace);
	// Close the datasets
	diffractionForceMagDataset.close();
	diffractionForcePhaDataset.close();
	froudeKrylovForceMagDataset.close();
	froudeKrylovForcePhaDataset.close();
	waveExcitingMagDataset.close();
	waveExcitingPhaDataset.close();
	// Allocate memory for the diffraction force magnitude cube
	pWaveDiffMag = new arma::cube;
	pWaveDiffPha = new arma::cube;
	pWaveFKMag = new arma::cube;
	pWaveFKPha = new arma::cube;
	pWaveExcitingMag = new arma::cube;
	pWaveExcitingPha = new arma::cube;
	pWaveDiffMag->set_size(6, numFrequencies, numHeadings);
	pWaveDiffPha->set_size(6, numFrequencies, numHeadings);
	pWaveFKMag->set_size(6, numFrequencies, numHeadings);
	pWaveFKPha->set_size(6, numFrequencies, numHeadings);
	pWaveExcitingMag->set_size(6, numFrequencies, numHeadings);
	pWaveExcitingPha->set_size(6, numFrequencies, numHeadings);
	for (int i = 0; i < numHeadings; i++)
	{
		for (int j = 0; j < numFrequencies; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				// Calculate the index for the 1D buffer [numHeadings, numBodies, numFrequencies, 6]
				int index = i * numBodies * numFrequencies * 6 + idBody * numFrequencies * 6 + j * 6 + k;
				// Copy data from buffer to arrays
				(*pWaveDiffMag)(k, j, i) = buffer_diffraction_force_mag[index];
				(*pWaveDiffPha)(k, j, i) = buffer_diffraction_force_pha[index];
				(*pWaveFKMag)(k, j, i) = buffer_froude_krylov_force_mag[index];
				(*pWaveFKPha)(k, j, i) = buffer_froude_krylov_force_pha[index];
				(*pWaveExcitingMag)(k, j, i) = buffer_wave_exciting_mag[index];
				(*pWaveExcitingPha)(k, j, i) = buffer_wave_exciting_pha[index];
			}
		}
	}
	// Change pWaveExciting to pWaveDiff if FK forces are computed as nonlinear forces
	// TODO: This should be done in a more consistent way with the ehydb format
	if (pBodies[idBody]->flag_hydrostatics == 2)
	{
		pWaveExcitingMag = pWaveDiffMag;
		pWaveExcitingPha = pWaveDiffPha;
	}
	// Free the buffer memory
	delete[] buffer_diffraction_force_mag;
	delete[] buffer_diffraction_force_pha;
	delete[] buffer_froude_krylov_force_mag;
	delete[] buffer_froude_krylov_force_pha;
	delete[] buffer_wave_exciting_mag;
	delete[] buffer_wave_exciting_pha;

	// Read the "hydstiffness" and "mass" datasets which are 3D arrays [numBodies, 6, 6]
	std::cout << "Reading hydrostatic stiffness and mass data...\n";
	// Open "hydstiffness"
	H5::DataSet hydrostaticStiffnessDataset = file.openDataSet("/hydstiffness");
	H5::DataSpace hydrostaticStiffnessSpace = hydrostaticStiffnessDataset.getSpace();
	hsize_t dims_hs[3];
	hydrostaticStiffnessSpace.getSimpleExtentDims(dims_hs, NULL);
	// Open "mass"
	H5::DataSet massDataset = file.openDataSet("/mass");
	H5::DataSpace massSpace = massDataset.getSpace();
	hsize_t dims_m[3];
	massSpace.getSimpleExtentDims(dims_m, NULL);
	// Allocate memory for the hydrostatic stiffness buffer [numBodies, 6, 6]
	num_values = numBodies * 6 * 6;
	double *buffer_hydrostatic_stiffness = new double[num_values];
	double *buffer_mass = new double[num_values];
	// Read the hydrostatic stiffness and mass data
	hydrostaticStiffnessDataset.read(buffer_hydrostatic_stiffness, H5::PredType::NATIVE_DOUBLE, hydrostaticStiffnessSpace, hydrostaticStiffnessSpace);
	massDataset.read(buffer_mass, H5::PredType::NATIVE_DOUBLE, massSpace, massSpace);
	// Close the datasets
	hydrostaticStiffnessDataset.close();
	massDataset.close();
	// Allocate memory for the hydrostatic stiffness and mass matrix
	pHydrostaticStiffness = new arma::mat;
	pStructuralMass = new arma::mat;
	pHydrostaticStiffness->set_size(6, 6);
	pStructuralMass->set_size(6, 6);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			// Calculate the index for the 1D buffer [numBodies, 6, 6]
			int ind = idBody * 6 * 6 + i * 6 + j;
			// Copy data from buffer to arrays
			(*pHydrostaticStiffness)(i, j) = buffer_hydrostatic_stiffness[ind];
			(*pStructuralMass)(i, j) = buffer_mass[ind];
		}
	}
	// Free the buffer memory
	delete[] buffer_hydrostatic_stiffness;
	delete[] buffer_mass;

	// Read "qtf_diff_mag", "qtf_diff_pha", "qtf_sum_mag" and "qtf_sum_pha"
	// datasets which are 6D arrays [numBodies, numHeadings, numHeadings, numFrequencies, numFrequencies, 6]
	// Check for existence of QTF data
	if (file.nameExists("/qtf_diff_mag") == 0)
	{
		// TODO: Review 2nd order excitation flag usage
		if (pBodies[idBody]->secondOrderExcitationFlag > 0)
		{
			std::stringstream ss;
			ss << "QTF data not found in file." << std::endl;
			throw IOError(ss.str());
		}
		else
		{
			std::cout << "WARNING: QTF data not found in file. \n";
		}
	}
	else
	{
		if (pBodies[idBody]->secondOrderExcitationFlag > 0 && pBodies[idBody]->secondOrderExcitationFlag < 0)
		{
			std::cout << "Reading QTF data...\n";
			// Open "qtf_diff_mag"
			H5::DataSet qtfDiffMagDataset = file.openDataSet("/qtf_diff_mag");
			H5::DataSpace qtfDiffMagSpace = qtfDiffMagDataset.getSpace();
			hsize_t dims_qdfm[6];
			qtfDiffMagSpace.getSimpleExtentDims(dims_qdfm, NULL);
			// Open "qtf_diff_pha"
			H5::DataSet qtfDiffPhaDataset = file.openDataSet("/qtf_diff_pha");
			H5::DataSpace qtfDiffPhaSpace = qtfDiffPhaDataset.getSpace();
			hsize_t dims_qdfp[6];
			qtfDiffPhaSpace.getSimpleExtentDims(dims_qdfp, NULL);
			// Open "qtf_sum_mag"
			H5::DataSet qtfSumMagDataset = file.openDataSet("/qtf_sum_mag");
			H5::DataSpace qtfSumMagSpace = qtfSumMagDataset.getSpace();
			hsize_t dims_qsm[6];
			qtfSumMagSpace.getSimpleExtentDims(dims_qsm, NULL);
			// Open "qtf_sum_pha"
			H5::DataSet qtfSumPhaDataset = file.openDataSet("/qtf_sum_pha");
			H5::DataSpace qtfSumPhaSpace = qtfSumPhaDataset.getSpace();
			hsize_t dims_qsp[6];
			qtfSumPhaSpace.getSimpleExtentDims(dims_qsp, NULL);
			// Allocate memory for the QTF buffers [numBodies, numHeadings, numHeadings, numFrequencies, numFrequencies, 6]
			num_values = numBodies * numHeadings * numHeadings * numFrequencies * numFrequencies * 6;
			double *buffer_qtf_diff_mag = new double[num_values];
			double *buffer_qtf_diff_pha = new double[num_values];
			double *buffer_qtf_sum_mag = new double[num_values];
			double *buffer_qtf_sum_pha = new double[num_values];
			// Read the QTF data
			qtfDiffMagDataset.read(buffer_qtf_diff_mag, H5::PredType::NATIVE_DOUBLE, qtfDiffMagSpace, qtfDiffMagSpace);
			qtfDiffPhaDataset.read(buffer_qtf_diff_pha, H5::PredType::NATIVE_DOUBLE, qtfDiffPhaSpace, qtfDiffPhaSpace);
			qtfSumMagDataset.read(buffer_qtf_sum_mag, H5::PredType::NATIVE_DOUBLE, qtfSumMagSpace, qtfSumMagSpace);
			qtfSumPhaDataset.read(buffer_qtf_sum_pha, H5::PredType::NATIVE_DOUBLE, qtfSumPhaSpace, qtfSumPhaSpace);
			// Close the datasets
			qtfDiffMagDataset.close();
			qtfDiffPhaDataset.close();
			qtfSumMagDataset.close();
			qtfSumPhaDataset.close();
			// Allocate memory for the QTF cubes
			pQtfDiff = new arma::cube **[2]; // 2 for real and imaginary
			pQtfSum = new arma::cube **[2];	 // 2 for real and imaginary
			for (int ipart = 0; ipart < 2; ipart++)
			{
				pQtfDiff[ipart] = new arma::cube *[6];
				pQtfSum[ipart] = new arma::cube *[6];
				for (int idof = 0; idof < 6; idof++)
				{
					pQtfDiff[ipart][idof] = new arma::cube;
					pQtfSum[ipart][idof] = new arma::cube;
					pQtfDiff[ipart][idof]->set_size(numFrequencies, numFrequencies, numHeadings);
					pQtfSum[ipart][idof]->set_size(numFrequencies, numFrequencies, numHeadings);
					for (int if1 = 0; if1 < numFrequencies; if1++)
					{
						for (int if2 = 0; if2 < numFrequencies; if2++)
						{
							for (int ih = 0; ih < numHeadings; ih++)
							{
								// Calculate the index for the 1D buffer [numBodies, numHeadings, numHeadings, numFrequencies, numFrequencies, 6]
								int ind = idBody * numHeadings * numHeadings * numFrequencies * numFrequencies * 6 +
										  ih * numHeadings * numFrequencies * numFrequencies * 6 +
										  ih * numFrequencies * numFrequencies * 6 +
										  if1 * numFrequencies * 6 + if2 * 6 + idof;
								// Copy data from buffer to arrays
								// For the first part (real), use cos(phase)
								// For the second part (imaginary), use sin(phase)
								if (ipart == 0)
								{
									(*pQtfDiff[ipart][idof])(if1, if2, ih) = buffer_qtf_diff_mag[ind] * cos(buffer_qtf_diff_pha[ind]);
									(*pQtfSum[ipart][idof])(if1, if2, ih) = buffer_qtf_sum_mag[ind] * cos(buffer_qtf_sum_pha[ind]);
								}
								else
								{
									(*pQtfDiff[ipart][idof])(if1, if2, ih) = buffer_qtf_diff_mag[ind] * sin(buffer_qtf_diff_pha[ind]);
									(*pQtfSum[ipart][idof])(if1, if2, ih) = buffer_qtf_sum_mag[ind] * sin(buffer_qtf_sum_pha[ind]);
								}
							}
						}
					}
				}
			}
			// Free the buffer memory
			delete[] buffer_qtf_diff_mag;
			delete[] buffer_qtf_diff_pha;
			delete[] buffer_qtf_sum_mag;
			delete[] buffer_qtf_sum_pha;
		}
	}

	// Read "mean_drift_mag" as 4D array [numHeadings, numBodies, numFrequencies, 6]
	// Check for existence of mean drift data
	if (file.nameExists("/mean_drift_mag") == 0)
	{
		// TODO: Review 2nd order excitation flag usage
		if (pBodies[idBody]->secondOrderExcitationFlag == 4)
		{
			std::stringstream ss;
			ss << "Mean drift data not found in file." << std::endl;
			throw IOError(ss.str());
		}
		else
		{
			std::cout << "WARNING: Mean drift data not found in file. \n";
		}
	}
	else
	{
		std::cout << "Reading mean drift data...\n";
		H5::DataSet meanDriftMagDataset = file.openDataSet("/mean_drift_mag");
		H5::DataSpace meanDriftMagSpace = meanDriftMagDataset.getSpace();
		hsize_t dims_mdm[4];
		meanDriftMagSpace.getSimpleExtentDims(dims_mdm, NULL);
		// Allocate memory for the mean drift magnitude buffer [numHeadings, numBodies, numFrequencies, 6]
		num_values = numHeadings * numBodies * numFrequencies * 6;
		double *buffer_mean_drift_mag = new double[num_values];
		// Read the mean drift magnitude data
		meanDriftMagDataset.read(buffer_mean_drift_mag, H5::PredType::NATIVE_DOUBLE, meanDriftMagSpace, meanDriftMagSpace);
		// Close the dataset
		meanDriftMagDataset.close();
		// Allocate memory for the mean drift cube
		pMeanDrift = new arma::cube(activeDofs, numFrequencies, numHeadings, arma::fill::zeros);
		for (int ihd = 0; ihd < numHeadings; ihd++)
		{
			for (int ifr = 0; ifr < numFrequencies; ifr++)
			{
				for (int idof = 0; idof < activeDofs; idof++)
				{
					// Calculate the index for the 1D buffer [numHeadings, numBodies, numFrequencies, 6]
					int ind = ihd * numBodies * numFrequencies * 6 +
							  idBody * numFrequencies * 6 + ifr * 6 + idof;
					// Copy data from buffer to arrays
					(*pMeanDrift)(idof, ifr, ihd) = buffer_mean_drift_mag[ind];
				}
			}
		}
		// Free the buffer memory
		delete[] buffer_mean_drift_mag;
	}

	// Close the file
	file.close();
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
	// TODO: Implement different levels of instant position (none, xy, xy+yaw)
	// WE_Real_wt = interp1(wrapToPi((*pHeadings) + yaw), WE_Real_w, wrapToPi(pWave->headings_piece));
	// WE_Imag_wt = interp1(wrapToPi((*pHeadings) + yaw), WE_Imag_w, wrapToPi(pWave->headings_piece));

	// Convert interpolated transfer functions to magnitude and phase
	arma::cube H_Mag = arma::sqrt(arma::pow(WE_Real_wt, 2) + arma::pow(WE_Imag_wt, 2));
	arma::cube H_Pha = arma::atan2(WE_Imag_wt, WE_Real_wt);

	// Permute the transfer functions to the correct order (frequencies, headings, dofs)
	H_Mag = permute(H_Mag, 213);
	H_Pha = permute(H_Pha, 213);

	// Extract the index of the required wave piece
	arma::uword ind_piece;
	bool flag_gap;
	if (pWave->num_pieces > 1)
	{
		arma::uvec ind_piece_vec = arma::find((pWave->time_ini - pWave->time_gap) < t, 1, "last");
		ind_piece = ind_piece_vec(0);
		// Check if current time is in a gap
		flag_gap = ((ind_piece > 0) && (t < pWave->time_ini(ind_piece)));
	}
	else
	{
		ind_piece = 0;
		flag_gap = false;
	}

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
		Cm = arma::sum(H_Mag_loc % amplitudes_w(ind_piece(0)) %
						   arma::cos(phases_w(ind_piece(0)) +
									 H_Pha_loc +
									 x * kx_w +
									 y * ky_w),
					   1);
		Sm = arma::sum(H_Mag_loc % amplitudes_w(ind_piece(0)) %
						   arma::sin(phases_w(ind_piece(0)) +
									 H_Pha_loc +
									 x * kx_w +
									 y * ky_w),
					   1);

		// Compute amplitude and phase of the excitation force for the current wave piece
		Am = arma::sqrt(arma::pow(Cm, 2) + arma::pow(Sm, 2));
		PHIm = arma::atan2(Sm, Cm);

		// Compute excitation force for the current wave piece
		F_piece = arma::as_scalar(arma::sum(Am % arma::cos(t * ang_freqs_w - PHIm), 0));

		if (flag_gap)
		{
			// Multiply the transfer functions with the wave amplitudes for the current wave piece
			Cm = arma::sum(H_Mag_loc % amplitudes_w(ind_piece(0) - 1) %
							   arma::cos(phases_w(ind_piece(0) - 1) +
										 H_Pha_loc +
										 x * kx_w +
										 y * ky_w),
						   1);
			Sm = arma::sum(H_Mag_loc % amplitudes_w(ind_piece(0) - 1) %
							   arma::sin(phases_w(ind_piece(0) - 1) +
										 H_Pha_loc +
										 x * kx_w +
										 y * ky_w),
						   1);

			// Compute amplitude and phase of the excitation force for the current wave piece
			Am = arma::sqrt(arma::pow(Cm, 2) + arma::pow(Sm, 2));
			PHIm = arma::atan2(Sm, Cm);

			// Compute excitation force for the current wave piece
			F_gap = arma::as_scalar(arma::sum(Am % arma::cos(t * ang_freqs_w - PHIm), 0));

			// Store the excitation force mixing linearly the current and previous piece
			// TODO: Mix pieces with a qubic function instead of linear
			Fe(ii, 0) = F_piece * (t - pWave->time_end(ind_piece - 1)) / pWave->time_gap +
						F_gap * (pWave->time_ini(ind_piece) - t) / pWave->time_gap;
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
	arma::cube temp_cube;
	arma::mat temp_mat_1, temp_mat_2;
	double F_piece, F_gap;

	// Interpolate transfer functions to the body heading
	arma::mat ***HDif = new arma::mat **[2];
	// arma::mat ***HSum = new arma::mat **[2];

	// TODO: Remove this when the precomputed second order excitation forces are implemented
	if (pBodies[idBody]->secondOrderExcitationFlag == 1)
	{
		HDif = QtfDiff_wt;
		// HSum = QtfSum_wt;
	}
	else
	{
		HDif = QtfDiff_wt;
		// HSum = QtfSum_wt;

		// TODO: Implement different levels of instant position (none, xy, xy+yaw)
		// TODO: Review this interpolation
		// for (int ii = 0; ii < 2; ii++)
		// {
		// 	HDif[ii] = new arma::mat *[activeDofs];
		// 	HSum[ii] = new arma::mat *[activeDofs];
		// 	for (int jj = 0; jj < activeDofs; jj++)
		// 	{
		// 		HDif[ii][jj] = new arma::mat;
		// 		HSum[ii][jj] = new arma::mat;
		// 		temp_cube = *QtfDiff_w[ii][jj];
		// 		*HDif[ii][jj] = interp1((*pHeadings) + yaw, permute(temp_cube, 312), pWave->heading);
		// 		temp_cube = *QtfSum_w[ii][jj];
		// 		*HSum[ii][jj] = interp1((*pHeadings) + yaw, permute(temp_cube, 312), pWave->heading);
		// 	}
		// }
	}

	// Extract the index of the required wave piece
	arma::uword ind_piece;
	bool flag_gap;
	if (pWave->num_pieces > 1)
	{
		arma::uvec ind_piece_vec = arma::find((pWave->time_ini - pWave->time_gap) < t, 1, "last");
		ind_piece = ind_piece_vec(0);
		// Check if current time is in a gap
		flag_gap = ((ind_piece > 0) && (t < pWave->time_ini(ind_piece)));
	}
	else
	{
		ind_piece = 0;
		flag_gap = false;
	}

	// Loop over all active degrees of freedom
	// TODO: This only works for 6 DOFs, change it.
	for (int ii = 0; ii < activeDofs; ii++)
	{
		// Initialize the excitation force for the current DOF
		F_piece = 0.0;

		// Argument for dif term
		temp_mat_1 = wD * t - phD(ind_piece) - kxD * x - kyD * y;
		// Accumulate the real part of the difference QTF for the current piece and DOF
		temp_mat_2 = (*HDif[0][ii]) % ampP(ind_piece) % arma::cos(temp_mat_1);
		F_piece = F_piece + 0.5 * arma::accu(temp_mat_2);
		// Accumulate the imaginary part of the difference QTF for the current piece and DOF
		temp_mat_2 = (*HDif[1][ii]) % ampP(ind_piece) % arma::sin(temp_mat_1);
		F_piece = F_piece + 0.5 * arma::accu(temp_mat_2);

		// // Argument for sum term
		// temp_mat_1 = phS(ind_piece) + kxS * x + kyS * y - wS * t;
		// // Accumulate the real part of the sum QTF for the current piece and DOF
		// temp_mat_2 = (*HSum[0][ii]) % ampP(ind_piece) % arma::cos(temp_mat_1);
		// F_piece = F_piece + arma::accu(temp_mat_2);
		// // Accumulate the imaginary part of the sum QTF for the current piece and DOF
		// temp_mat_2 = (*HSum[1][ii]) % ampP(ind_piece) % arma::sin(temp_mat_1);
		// F_piece = F_piece - arma::accu(temp_mat_2);

		if (flag_gap)
		{
			// Compute excitation force for the current wave piece
			F_gap = 0.0;

			// Argument for dif term
			temp_mat_1 = wD * t - phD(ind_piece - 1) - kxD * x - kyD * y;
			// Accumulate the real part of the difference QTF for the current piece and DOF
			temp_mat_2 = (*HDif[0][ii]) % ampP(ind_piece - 1) % arma::cos(temp_mat_1);
			F_gap = F_gap + 0.5 * arma::accu(temp_mat_2);
			// Accumulate the imaginary part of the difference QTF for the current piece and DOF
			temp_mat_2 = (*HDif[1][ii]) % ampP(ind_piece - 1) % arma::sin(temp_mat_1);
			F_gap = F_gap + 0.5 * arma::accu(temp_mat_2);

			// // Argument for sum term
			// temp_mat_1 = phS(ind_piece - 1) + kxS * x + kyS * y - wS * t;
			// // Accumulate the real part of the sum QTF for the current piece and DOF
			// temp_mat_2 = (*HSum[0][ii]) % ampP(ind_piece - 1) % arma::cos(temp_mat_1);
			// F_gap = F_gap + arma::accu(temp_mat_2);
			// // Accumulate the imaginary part of the sum QTF for the current piece and DOF
			// temp_mat_2 = (*HSum[1][ii]) % ampP(ind_piece - 1) % arma::sin(temp_mat_1);
			// F_gap = F_gap - arma::accu(temp_mat_2);

			// Store the excitation force mixing linearly the current and previous piece
			// TODO: Mix pieces with a qubic function instead of linear
			Fe(ii, 0) = F_piece * (t - pWave->time_end(ind_piece - 1)) / pWave->time_gap +
						F_gap * (pWave->time_ini(ind_piece) - t) / pWave->time_gap;
		}
		else
		{
			// Store the excitation force of the current piece
			Fe(ii, 0) = F_piece;
		}
	}

	// Rotate the excitation force to the fixed frame
	// TODO: This should only be done if yaw instant position is used
	// double Fx = arma::as_scalar(Fe(0, 0));
	// double Fy = arma::as_scalar(Fe(1, 0));
	// Fe(0, 0) = Fx * cos(yaw) - Fy * sin(yaw);
	// Fe(1, 0) = Fx * sin(yaw) + Fy * cos(yaw);

	// Store the excitation force in the body variable for plotting
	pBodies[idBody]->excitationForces_2 = Fe;

	// Return the excitation force
	return Fe;
}

void HydroDatabase::SetUp(void)
{
	// TODO: Implement a logger with different levels of verbosity
	std::cout << "HydroDatabase::SetUp - At first" << std::endl;
	// Load the wave pointer as a local variable
	Wave *pWave = pSim->pWave;

	std::cout << "HydroDatabase::SetUp - Convert from real/imag to mag/pha the first order transfer functions" << std::endl;
	// Convert from real/imag to mag/pha the first order transfer functions
	arma::cube WE_Real = (*pWaveExcitingMag) % arma::cos((*pWaveExcitingPha));
	arma::cube WE_Imag = (*pWaveExcitingMag) % arma::sin((*pWaveExcitingPha));

	// Check that the transfer functions frequencies cover the wave frequencies, return a warning and crop frequencies outside if not
	arma::uvec ind_wave_freqs = arma::regspace<arma::uvec>(0, pWave->num_comps_piece - 1);
	if (((*pFrequencies).min() > (pWave->freqs_piece).min()) || ((*pFrequencies).max() < (pWave->freqs_piece).max()))
	{
		std::cout << "    --> WARNING: The frequencies provided in the hydrodinamic data base do not cover properly the wave!" << std::endl;
		ind_wave_freqs = arma::find((pWave->freqs_piece >= (*pFrequencies).min()) && (pWave->freqs_piece <= (*pFrequencies).max()));
	}
	numFrequencies_w = ind_wave_freqs.n_elem;
	freqs_w = pWave->freqs_piece.elem(ind_wave_freqs);
	ang_freqs_w = pWave->ang_freqs_piece.elem(ind_wave_freqs);
	kx_w = pWave->kx_piece.rows(ind_wave_freqs);
	ky_w = pWave->ky_piece.rows(ind_wave_freqs);
	amplitudes_w = arma::field<arma::mat>(pWave->num_pieces);
	phases_w = arma::field<arma::mat>(pWave->num_pieces);
	for (int kk = 0; kk < pWave->num_pieces; kk++)
	{
		amplitudes_w(kk) = pWave->amplitudes_piece(kk).rows(ind_wave_freqs);
		phases_w(kk) = pWave->phases_piece(kk).rows(ind_wave_freqs);
	}

	std::cout << "HydroDatabase::SetUp - Interpolate the first order transfer functions to the wave frequencies" << std::endl;
	// Interpolate the first order transfer functions to the wave frequencies
	WE_Real_w = interp1(*pFrequencies, permute(WE_Real, 231), freqs_w);
	WE_Imag_w = interp1(*pFrequencies, permute(WE_Imag, 231), freqs_w);
	WE_Real_w = permute(WE_Real_w, 213);
	WE_Imag_w = permute(WE_Imag_w, 213);

	// Interpolate the first order transfer functions to the wave headings
	WE_Real_wt = interp1(wrapToPi((*pHeadings)), WE_Real_w, wrapToPi(pWave->headings_piece));
	WE_Imag_wt = interp1(wrapToPi((*pHeadings)), WE_Imag_w, wrapToPi(pWave->headings_piece));

	// Preprocess for QTFs
	// TODO: QTFs should only be interpolated if they will be used
	if (pBodies[idBody]->secondOrderExcitationFlag > 0 && pBodies[idBody]->secondOrderExcitationFlag < 4)
	{
		std::cout << "HydroDatabase::SetUp - Interpolate the second order transfer functions to the wave frequencies" << std::endl;
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
				*QtfDiff_w[ii][jj] = interp2(*pFrequencies, *pFrequencies, temp, freqs_w, freqs_w);

				temp = *pQtfSum[ii][jj];
				*QtfSum_w[ii][jj] = interp2(*pFrequencies, *pFrequencies, temp, freqs_w, freqs_w);
			}
		}

		// Interpolate the second order transfer functions to the wave headings
		QtfDiff_wt = new arma::mat **[2];
		QtfSum_wt = new arma::mat **[2];
		for (int ii = 0; ii < 2; ii++)
		{
			QtfDiff_wt[ii] = new arma::mat *[activeDofs];
			QtfSum_wt[ii] = new arma::mat *[activeDofs];
			for (int jj = 0; jj < activeDofs; jj++)
			{
				QtfDiff_wt[ii][jj] = new arma::mat;
				QtfSum_wt[ii][jj] = new arma::mat;

				temp = *QtfDiff_w[ii][jj];
				*QtfDiff_wt[ii][jj] = interp1(*pHeadings, permute(temp, 312), pWave->heading);

				temp = *QtfSum_w[ii][jj];
				*QtfSum_wt[ii][jj] = interp1(*pHeadings, permute(temp, 312), pWave->heading);
			}
		}


		std::cout << "HydroDatabase::SetUp - Define the matrices required for time domain QTF forces computation" << std::endl;
		// Preprocess the matrices required for time domain QTF forces computation
		ampP = arma::field<arma::mat>(pWave->num_pieces);
		phS = arma::field<arma::mat>(pWave->num_pieces);
		phD = arma::field<arma::mat>(pWave->num_pieces);

		wS = arma::zeros(numFrequencies_w, numFrequencies_w);
		kxS = arma::zeros(numFrequencies_w, numFrequencies_w);
		kyS = arma::zeros(numFrequencies_w, numFrequencies_w);
		wD = arma::zeros(numFrequencies_w, numFrequencies_w);
		kxD = arma::zeros(numFrequencies_w, numFrequencies_w);
		kyD = arma::zeros(numFrequencies_w, numFrequencies_w);

		std::cout << "HydroDatabase::SetUp - Preprocess the matrices required for time domain QTF forces computation (w, k)" << std::endl;
		int i_wave, j_wave;
		for (int ii = 0; ii < numFrequencies_w; ii++)
		{
			i_wave = ind_wave_freqs(ii);
			for (int jj = 0; jj < numFrequencies_w; jj++)
			{
				j_wave = ind_wave_freqs(jj);
				wS(ii, jj) = ang_freqs_w(ii) + ang_freqs_w(jj);
				kxS(ii, jj) = pWave->kx_1D_piece(i_wave) + pWave->kx_1D_piece(j_wave);
				kyS(ii, jj) = pWave->ky_1D_piece(i_wave) + pWave->ky_1D_piece(j_wave);
				wD(ii, jj) = ang_freqs_w(ii) - ang_freqs_w(jj);
				kxD(ii, jj) = pWave->kx_1D_piece(i_wave) - pWave->kx_1D_piece(j_wave);
				kyD(ii, jj) = pWave->ky_1D_piece(i_wave) - pWave->ky_1D_piece(j_wave);
			}
		}

		std::cout << "HydroDatabase::SetUp - Preprocess the matrices required for time domain QTF forces computation (amp, pha)" << std::endl;
		for (int kk = 0; kk < pWave->num_pieces; kk++)
		{
			ampP(kk) = arma::zeros(numFrequencies_w, numFrequencies_w);
			phS(kk) = arma::zeros(numFrequencies_w, numFrequencies_w);
			phD(kk) = arma::zeros(numFrequencies_w, numFrequencies_w);
			for (int ii = 0; ii < numFrequencies_w; ii++)
			{
				i_wave = ind_wave_freqs(ii);
				for (int jj = 0; jj < numFrequencies_w; jj++)
				{
					j_wave = ind_wave_freqs(jj);
					ampP(kk)(ii, jj) = pWave->amplitudes_1D_piece(kk)(i_wave) * pWave->amplitudes_1D_piece(kk)(j_wave);
					phS(kk)(ii, jj) = pWave->phases_1D_piece(kk)(i_wave) + pWave->phases_1D_piece(kk)(j_wave);
					phD(kk)(ii, jj) = pWave->phases_1D_piece(kk)(i_wave) - pWave->phases_1D_piece(kk)(j_wave);
				}
			}
		}
	}

	if (pBodies[idBody]->secondOrderExcitationFlag > 2)
	{
		std::cout << "HydroDatabase::SetUp - Compute mean drift force" << std::endl;
		// Compute mean drift force
		// TODO: Review this!
		arma::cube temp_mD = interp2(*pFrequencies, *pHeadings, permute(*pMeanDrift, 231), freqs_w, pWave->headings_piece);
		for (int ii = 0; ii < activeDofs; ii++)
		{
			temp_mD.slice(ii) = temp_mD.slice(ii) % amplitudes_w(0) % amplitudes_w(0);
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

	std::cout << "HydroDatabase::SetUp - Set up the hydrodynamic forces at time zero" << std::endl;
	// Set up the hydrodynamic forces at time zero
	if (arma::accu(pWave->amplitudes) > 0)
	{

		pBodies[idBody]->Fb = CalculateHydrodynamicForces(0.0);
	}
	else
	{
		pBodies[idBody]->Fb = arma::zeros(activeDofs, 1);
	}

	std::cout << "HydroDatabase::SetUp - Precompute hydrodynamic forces time series" << std::endl;
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

	std::cout << "HydroDatabase::SetUp - At end" << std::endl;
}

arma::mat HydroDatabase::ComputeMeanDrift(void)
{

	// Compute mean drift force

	// TODO: Implement meandrift with instant position
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

	if (pBodies[idBody]->flag_hydrostatics == 1)
	{

		// Non-linear hydrostatic forces without wave
		pressure = -z * pSim->gravity * pSim->waterDensity;
	}
	else if (pBodies[idBody]->flag_hydrostatics == 2)
	{

		// Non-linear hydrostatic forces with wave
		arma::mat x = pBodies[idBody]->pNLHSMesh->nodes.col(0);
		arma::mat y = pBodies[idBody]->pNLHSMesh->nodes.col(1);
		pressure = pSim->pWave->GetPressure(t, x, y, z);
	}

	return pressure;
}
