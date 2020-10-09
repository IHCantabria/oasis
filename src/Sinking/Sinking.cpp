
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


Sinking::Sinking(int n, Simulation* pSim_inp){
	indBody = n;
	pSim = pSim_inp;
	rhoW = pSim->waterDensity;
	pSinkingBody = pSim->pBodies[indBody-1];
}


void Sinking::ReadPropertiesASCII(FILE* pFile, std::string inputFolderPath){

	char buffer_line [1000];
	char cHydroDatabaseName [1000];
	double dtemp;
	std::string file_path;

	// Read bodyReferenceMassMat
	for (int ii=0; ii<6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1) 
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the reference mass matrix of sinking body: " << indBody << "\n";
			throw ValueError(ss.str());
		}
		bodyReferenceMassMat(ii,ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	bodyReferenceMass = arma::as_scalar(bodyReferenceMassMat(0,0));

	// Read number of HDBs
	if (fscanf(pFile, "%d %[^\n]\n", &numHDBs, buffer_line) != 2)
	{
		std::cout << "Hey 1" << std::endl;
		std::stringstream ss;
		ss << "An error ocurred when trying to read the number of HDBs of sinking body: " << indBody << "\n";
		throw ValueError(ss.str());
	}

	if (numHDBs<1){
		std::stringstream ss;
		ss << "ERROR: The number of HDBs must be at least 1 " << indBody <<"\n";
		throw ValueError(ss.str());
	}

	//Read InterpMasses
	InterpMasses = arma::zeros(numHDBs,1);
	for (int ii=0; ii<numHDBs; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1) 
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the load status of sinking body: " << indBody <<"\n";
			throw ValueError(ss.str());
		}
		InterpMasses(ii,0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	pHydro = new HydroDatabase*[numHDBs];

	// Read and load HDBs filenames
	for (int ii=0; ii<numHDBs; ii++)
	{
		pHydro[ii] = new HydroDatabase(ii, indBody-1 , pSim->pBodies, pSim);
		if (fscanf(pFile, "%s %[^\n]\n", cHydroDatabaseName, buffer_line) != 2)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the file name of HDBs of sinking body: " << indBody <<"\n";
			throw ValueError(ss.str());
		}
		file_path = JoinPath(inputFolderPath, cHydroDatabaseName);
		pHydro[ii]->LoadHydrodynamicData(file_path);
	}

	// Read number of filling times
	if (fscanf(pFile, "%d %[^\n]\n", &numTimes, buffer_line) != 2) 
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the number of filling times of sinking body: " << indBody <<"\n";
		throw ValueError(ss.str());
	}

	if (numTimes<1){
		std::stringstream ss;
		ss << "ERROR: The number of time points must be at least 1 " << indBody <<"\n";
		throw ValueError(ss.str());
	}

	// Read groupsFillingTimes
	groupsFillingTimes = arma::zeros(numTimes,1);
	for (int ii=0; ii<numTimes; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1) 
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the filling times of sinking body: " << indBody <<"\n";
			throw ValueError(ss.str());
		}
		groupsFillingTimes(ii,0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read number of groups
	if (fscanf(pFile, "%d %[^\n]\n", &numGroups, buffer_line) != 2) 
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the number of filling times of sinking body: " << indBody <<"\n";
		throw ValueError(ss.str());
	}

	if (numGroups<1){
		std::stringstream ss;
		ss << "ERROR: The number of groups must be at least 1 " << indBody <<"\n";
		throw ValueError(ss.str());
	}

	groupsSizes = arma::zeros(numGroups,2);
	groupsAreas = arma::zeros(numGroups,1);
	groupsCenters = arma::zeros(numGroups,3);
	groupsFillingStates = arma::zeros(numGroups,numTimes);
	for (int ii=0; ii<numGroups; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile); // Ignore one line

		if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2) 
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the group size of sinking body: " << indBody <<"\n";
			throw ValueError(ss.str());
		}
		groupsSizes(ii,0) = dtemp;

		if (fscanf(pFile, "%lf %[^\n]\n", &dtemp, buffer_line) != 2) 
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the group size of sinking body: " << indBody <<"\n";
			throw ValueError(ss.str());
		}
		groupsSizes(ii,1) = dtemp;

		for (int jj=0; jj<3; jj++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1) 
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the group centers of sinking body: " << indBody <<"\n";
				throw ValueError(ss.str());
			}
			groupsCenters(ii,jj) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);

		for (int jj=0; jj<numTimes; jj++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1) 
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the filling states sinking body: " << indBody <<"\n";
				throw ValueError(ss.str());
			}
			groupsFillingStates(ii,jj) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);

	}

	groupsAreas = prod(groupsSizes,1);
}


void Sinking::UpdateSinkingHydrodynamics(double t){

	UpdateGroupsFillingState(t);

	if (numHDBs>1) {
		UpdateInterpHydro();
		pSinkingBody->pHydro->InterpolateHydro(pHydro[indHydro1], pHydro[indHydro2], hydroInterpCoef);
	}

	UpdateBodyProperties();
}


void Sinking::UpdateSinkingHydrostatics(double t){

	UpdateGroupsFillingState(t);

	if (numHDBs>1) {

		UpdateInterpHydro();

		// No es necesario actualizar la posición de equilibrio, se toma la posicion de equilibrio inicial y se
		// aplica la fuerza de la masa del llenado.
		//pSinkingBody->pos_eq.rows(0, 2) = pHydro[indHydro1]->cog*(1-hydroInterpCoef) +
		//                                  pHydro[indHydro2]->cog*hydroInterpCoef;

		arma::mat newHydrostaticStiffness = *(pHydro[indHydro1]->pHydrostaticStiffness)*(1-hydroInterpCoef) +
		                                    *(pHydro[indHydro2]->pHydrostaticStiffness)*hydroInterpCoef;
		pSinkingBody->pHydro->UpdateHydroStiffness(newHydrostaticStiffness);
	}

	pSinkingBody->filling_mass = totalFillingMass;
	pSinkingBody->pos_filling_cog = groupsCOG;
}


void Sinking::UpdateGroupsFillingState(double t){

	arma::mat time = arma::ones(1,1)*t;
	arma::mat groupsMasses = interp1(groupsFillingTimes,groupsFillingStates,time); groupsMasses = groupsMasses.t();
	totalFillingMass = arma::accu(groupsMasses);

	arma::mat masses_3 = arma::join_horiz(groupsMasses,groupsMasses,groupsMasses);
	arma::mat groupsHeights = (groupsMasses/rhoW)/groupsAreas;
	arma::mat centers = groupsCenters;
	centers.col(2) = centers.col(2) + groupsHeights;
	groupsCOG = arma::trans(arma::sum(centers%masses_3)/totalFillingMass);

	arma::mat Ix = groupsMasses%(arma::pow(groupsSizes.col(1),2) + arma::pow(groupsHeights,2))/12 + 
	               groupsMasses%(arma::pow(centers.col(1),2)+arma::pow(centers.col(2),2));
	arma::mat Iy = groupsMasses%(arma::pow(groupsSizes.col(0),2) + arma::pow(groupsHeights,2))/12 + 
	               groupsMasses%(arma::pow(centers.col(0),2)+arma::pow(centers.col(2),2));
	arma::mat Iz = groupsMasses%(arma::pow(groupsSizes.col(0),2) + arma::pow(groupsSizes.col(1),2))/12 + 
	               groupsMasses%(arma::pow(centers.col(0),2)+arma::pow(centers.col(1),2));
	groupsInertia = arma::eye(6,6)*totalFillingMass;
	groupsInertia(3,3) = arma::accu(Ix); groupsInertia(4,4) = arma::accu(Iy); groupsInertia(5,5) = arma::accu(Iz);
}


void Sinking::UpdateInterpHydro(void){

	arma::uvec ind = find(InterpMasses >= totalFillingMass, 1);

	indHydro2 = arma::as_scalar(ind(0));
	indHydro1 = indHydro2 - 1;

	double mass1 = arma::as_scalar(InterpMasses(indHydro1,0));
	double mass2 = arma::as_scalar(InterpMasses(indHydro2,0));

	hydroInterpCoef = (totalFillingMass-mass1)/(mass2-mass1);
}


void Sinking::UpdateBodyProperties(void){

	arma::mat newStructuralMass = bodyReferenceMassMat + groupsInertia;

	pSinkingBody->filling_mass = totalFillingMass;
	pSinkingBody->pos_filling_cog = groupsCOG;
	pSinkingBody->inertia = newStructuralMass;
	pSinkingBody->pHydro->UpdateStructuralMass(newStructuralMass);
	pSinkingBody->pHydro->UpdateTotalMass();
}