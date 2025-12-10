
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <ctime>
#include "Sinking.hpp"
#include "../CommonTools.hpp"
#include "../Simulations/Simulation.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Bodies/Bodies.hpp"
#include "../Hydro/HydroForce.hpp"
#include "../MathTools.hpp"
#include "../Hydro/HydroDatabase.hpp"

Sinking::Sinking(int n, Simulation *pSim_inp)
{
	indBody = n;
	pSim = pSim_inp;
	rhoW = pSim->waterDensity;
	pSinkingBody = pSim->pBodies[indBody - 1];
}

void Sinking::ReadPropertiesASCII(FILE *pFile, std::string inputFolderPath)
{

	char buffer_line[1000];
	char cHydroDatabaseName[1000];
	double dtemp;
	int itemp;
	arma::mat mtemp;
	std::string file_path;

	// Read bodyReferenceMassMat
	for (int ii = 0; ii < 6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the reference mass matrix of sinking body: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		bodyReferenceMassMat(ii, ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	bodyReferenceMass = arma::as_scalar(bodyReferenceMassMat(0, 0));

	// Read number of HDBs
	if (fscanf(pFile, "%d %[^\n]\n", &numHDBs, buffer_line) != 2)
	{
		std::cout << "Hey 1" << std::endl;
		std::stringstream ss;
		ss << "An error occurred when trying to read the number of HDBs of sinking body: " << indBody << "\n";
		throw ValueError(ss.str());
	}

	if (numHDBs < 1)
	{
		std::stringstream ss;
		ss << "ERROR: The number of HDBs must be at least 1 " << indBody << "\n";
		throw ValueError(ss.str());
	}

	// Read InterpMasses
	InterpMasses = arma::zeros(numHDBs, 1);
	for (int ii = 0; ii < numHDBs; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the load status of sinking body: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		InterpMasses(ii, 0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	pHydro = new HydroDatabase *[numHDBs];

	// Read and load HDBs filenames
	for (int ii = 0; ii < numHDBs; ii++)
	{
		pHydro[ii] = new HydroDatabase(pSinkingBody->hydroDatabaseIndex, indBody - 1, pSim->pBodies, pSim);
		if (fscanf(pFile, "%s %[^\n]\n", cHydroDatabaseName, buffer_line) != 2)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the file name of HDBs of sinking body: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		file_path = JoinPath(inputFolderPath, cHydroDatabaseName);
		std::cout << "--> Loading Sinking HDB " << ii + 1 << "..." << std::endl;
		pHydro[ii]->LoadHydrodynamicData(file_path);
		std::cout << "--> Done loading Sinking HDB " << ii + 1 << std::endl;
	}

	// Read number of groups
	if (fscanf(pFile, "%d %[^\n]\n", &numGroups, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the number of groups times of sinking body: " << indBody << "\n";
		throw ValueError(ss.str());
	}

	if (numGroups < 1)
	{
		std::stringstream ss;
		ss << "ERROR: The number of groups must be at least 1 " << indBody << "\n";
		throw ValueError(ss.str());
	}

	groupsPoints.set_size(numGroups, 1);
	groupsCenters = arma::zeros(numGroups, 3);
	groupsAreas = arma::zeros(numGroups, 1);
	groupsIx = arma::zeros(numGroups, 1);
	groupsIy = arma::zeros(numGroups, 1);
	groupsFillingTimes.set_size(numGroups, 1);
	groupsFillingStates.set_size(numGroups, 1);

	for (int ii = 0; ii < numGroups; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile); // Ignore one line

		if (fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line) != 2)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the number of group polygon points: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		if (itemp < 1)
		{
			std::stringstream ss;
			ss << "ERROR: The number of group polygon points must be at least 1 " << indBody << "\n";
			throw ValueError(ss.str());
		}

		mtemp = arma::zeros(itemp + 1, 2);
		for (int jj = 0; jj < itemp; jj++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the group polygon X coordinates: " << indBody << "\n";
				throw ValueError(ss.str());
			}
			mtemp(jj, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		mtemp(itemp, 0) = mtemp(0, 0);
		for (int jj = 0; jj < itemp; jj++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the group polygon Y coordinates: " << indBody << "\n";
				throw ValueError(ss.str());
			}
			mtemp(jj, 1) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		mtemp(itemp, 1) = mtemp(0, 1);
		groupsPoints(ii, 0) = mtemp;

		for (int jj = 0; jj < itemp; jj++)
		{
			dtemp = (mtemp(jj, 0) * mtemp(jj + 1, 1) - mtemp(jj + 1, 0) * mtemp(jj, 1));
			groupsAreas(ii, 0) = groupsAreas(ii, 0) + dtemp;
			groupsIx(ii, 0) = groupsIx(ii, 0) + dtemp * (mtemp(jj, 1) * mtemp(jj, 1) + mtemp(jj, 1) * mtemp(jj + 1, 1) + mtemp(jj + 1, 1) * mtemp(jj + 1, 1));
			groupsIy(ii, 0) = groupsIy(ii, 0) + dtemp * (mtemp(jj, 0) * mtemp(jj, 0) + mtemp(jj, 0) * mtemp(jj + 1, 0) + mtemp(jj + 1, 0) * mtemp(jj + 1, 0));
			groupsCenters(ii, 0) = groupsCenters(ii, 0) + dtemp * (mtemp(jj, 0) + mtemp(jj + 1, 0));
			groupsCenters(ii, 1) = groupsCenters(ii, 1) + dtemp * (mtemp(jj, 1) + mtemp(jj + 1, 1));
		}
		groupsAreas(ii, 0) = groupsAreas(ii, 0) / 2.0;
		groupsIx(ii, 0) = groupsIx(ii, 0) / 12.0;
		groupsIy(ii, 0) = groupsIy(ii, 0) / 12.0;
		groupsCenters(ii, 0) = groupsCenters(ii, 0) / (6.0 * groupsAreas(ii, 0));
		groupsCenters(ii, 1) = groupsCenters(ii, 1) / (6.0 * groupsAreas(ii, 0));

		if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the group floor Z coordinate: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		groupsCenters(ii, 2) = dtemp;

		std::cout << "        For group " << ii << ", the area is: " << groupsAreas(ii, 0) << std::endl;
		std::cout << "                     the center is at: " << groupsCenters.row(ii);
		std::cout << "                     and the points read are: " << std::endl
				  << mtemp.t() << std::endl;

		if (fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line) != 2)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the number of filling times: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		if (itemp < 1)
		{
			std::stringstream ss;
			ss << "ERROR: The number of filling times must be at least 1 " << indBody << "\n";
			throw ValueError(ss.str());
		}

		mtemp = arma::zeros(itemp, 1);
		for (int jj = 0; jj < itemp; jj++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the group filling times: " << indBody << "\n";
				throw ValueError(ss.str());
			}
			mtemp(jj, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		groupsFillingTimes(ii, 0) = mtemp;

		mtemp = arma::zeros(itemp, 1);
		for (int jj = 0; jj < itemp; jj++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the group filling states: " << indBody << "\n";
				throw ValueError(ss.str());
			}
			mtemp(jj, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		groupsFillingStates(ii, 0) = mtemp;
	}
}

void Sinking::UpdateSinkingHydrodynamics(double t)
{

	// std::cout << "        UpdateGroupsFillingState(t);" << std::endl;
	UpdateGroupsFillingState(t);

	if (numHDBs > 1)
	{
		// std::cout << "        UpdateInterpHydro();" << std::endl;
		UpdateInterpHydro();
		// std::cout << "        pSinkingBody->pHydro->InterpolateHydro(pHydro[indHydro1], pHydro[indHydro2], hydroInterpCoef);" << std::endl;
		pSinkingBody->pHydro->InterpolateHydro(pHydro[indHydro1], pHydro[indHydro2], hydroInterpCoef);
	}

	// std::cout << "        UpdateBodyProperties();" << std::endl;
	UpdateBodyProperties();
}

void Sinking::UpdateSinkingHydrostatics(double t)
{

	UpdateGroupsFillingState(t);

	if (numHDBs > 1)
	{

		UpdateInterpHydro();

		arma::mat newHydrostaticStiffness = *(pHydro[indHydro1]->pHydrostaticStiffness) * (1 - hydroInterpCoef) +
											*(pHydro[indHydro2]->pHydrostaticStiffness) * hydroInterpCoef;
		pSinkingBody->pHydro->UpdateHydroStiffness(newHydrostaticStiffness);
	}

	pSinkingBody->filling_mass = totalFillingMass;
	pSinkingBody->pos_filling_cog = groupsCOG;
}

void Sinking::UpdateGroupsFillingState(double t)
{

	// std::cout << "Sinking::UpdateGroupsFillingState - Start" << std::endl;

	// std::cout << "Sinking::UpdateGroupsFillingState - Interpolate Filling State" << std::endl;
	arma::mat time = arma::ones(1, 1) * t;
	groupsMasses = arma::zeros(numGroups, 1);
	for (int ii = 0; ii < numGroups; ii++)
	{
		groupsMasses(ii, 0) = arma::as_scalar(interp1(groupsFillingTimes(ii, 0), groupsFillingStates(ii, 0), time));
	}
	totalFillingMass = arma::accu(groupsMasses);

	// std::cout << "Sinking::UpdateGroupsFillingState - Compute Heights and COGs" << std::endl;
	arma::mat masses_3 = arma::join_horiz(groupsMasses, groupsMasses, groupsMasses);
	arma::mat groupsHeights = (groupsMasses / rhoW) / groupsAreas;
	arma::mat centers = groupsCenters;
	centers.col(2) = centers.col(2) + groupsHeights * 0.5;
	groupsCOG = arma::trans(arma::sum(centers % masses_3) / totalFillingMass);

	// std::cout << "Sinking::UpdateGroupsFillingState - Compute Inertias" << std::endl;
	arma::mat temp = groupsMasses % (groupsHeights % groupsHeights / 3.0 + groupsHeights % centers.col(2) + centers.col(2) % centers.col(2));
	arma::mat Ix = rhoW * groupsHeights % groupsIx + temp;
	arma::mat Iy = rhoW * groupsHeights % groupsIy + temp;
	arma::mat Iz = rhoW * groupsHeights % (groupsIx + groupsIy);
	groupsInertia = arma::eye(6, 6) * totalFillingMass;
	groupsInertia(3, 3) = arma::accu(Ix);
	groupsInertia(4, 4) = arma::accu(Iy);
	groupsInertia(5, 5) = arma::accu(Iz);

	// std::cout << "Sinking::UpdateGroupsFillingState - End" << std::endl;
}

void Sinking::UpdateInterpHydro(void)
{

	arma::uvec ind = find(InterpMasses > totalFillingMass, 1);

	if (ind.is_empty())
	{

		indHydro2 = numHDBs - 1;
		indHydro1 = indHydro2 - 1;
		hydroInterpCoef = 1.0;
	}
	else
	{

		indHydro2 = arma::as_scalar(ind(0));
		indHydro1 = indHydro2 - 1;

		double mass1 = arma::as_scalar(InterpMasses(indHydro1, 0));
		double mass2 = arma::as_scalar(InterpMasses(indHydro2, 0));

		hydroInterpCoef = (totalFillingMass - mass1) / (mass2 - mass1);
	}
}

void Sinking::UpdateBodyProperties(void)
{

	arma::mat newStructuralMass = bodyReferenceMassMat + groupsInertia;

	pSinkingBody->filling_mass = totalFillingMass;
	pSinkingBody->pos_filling_cog = groupsCOG;
	pSinkingBody->inertia = newStructuralMass;
	pSinkingBody->pHydro->UpdateStructuralMass(newStructuralMass);
	pSinkingBody->pHydro->UpdateTotalMass();
}

void Sinking::OpenOutputFilesASCII(std::string path)
{
	char buffer1[50];

	int nn1 = sprintf(buffer1, "SinkingFillingCOG_Body_%d.txt", indBody - 1);

	std::string file_path1 = JoinPath(path, buffer1);

	pfile_FillingCOG = fopen(file_path1.c_str(), "w");
	if (pfile_FillingCOG == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn1 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}
}

void Sinking::CloseOutputFilesASCII(void)
{
	fclose(pfile_FillingCOG);
}

// Escibir datos a fichero
void Sinking::WriteOut(double t)
{
	fprintf(pfile_FillingCOG, "%f    %f    %f    %f    %f \n", t, groupsCOG(0, 0), groupsCOG(1, 0), groupsCOG(2, 0), totalFillingMass);
}