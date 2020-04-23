/*
Libreria para las condiciones de contorno
*/

#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <armadillo>
#include "BCPs.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"


////////////////////////////////////////////////////////////////////////////
/////////////////////////// BCP CLASS DEFINITION ///////////////////////////
////////////////////////////////////////////////////////////////////////////
BCP::BCP(int incId)
{
	id = incId;
}


int BCP::GetId(void)
{
	return this->id;
}


int BCP::GetType(void)
{
	return this->typeBcp;
}


void BCP::Initialize(void)
{
	double a = 0.0;
}


void BCP::Print(void)
{
	printf("BCP: %d PROPERTIES:\n", this->GetId());
	printf("--> PosX: %f - PosY: %f - PosZ: %f\n", this->posG_BCP[0], this->posG_BCP[1], this->posG_BCP[2]);
	printf("--> Winch Id: %d\n\n", this->winchId);
}


void BCP::ReadPropertiesASCII(FILE* &pFilePointer)
{
	// Declare variables
	std::string header_check;
	char buffer_line [1000];
	char cActuatorFileName [1000];

	//Ignoro las tres primeras lineas, donde pone "New Body"

	for(int ii=0; ii<3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFilePointer);
		header_check = buffer_line;
		if (header_check.substr(0, 3).compare("///"))
		{
			std::stringstream ss;
			ss << "Error while parsing file: datosBCPs.dat - HINT BCP ID: " << this->GetId() << " - Please check that each type of BCP has its correct number of inputs.";
			throw IOError(ss.str());
		}
	}

	// Read BCP position
	fscanf(pFilePointer, "%lf %lf %lf %[^\n]\n", &pos(0, 0), &pos(1, 0), &pos(2, 0), buffer_line);

	// Read Wind ID
	fscanf(pFilePointer, "%d %[^\n]\n", &winchId, buffer_line);

	// Read actuator if any
	if (this->GetType() == 2)
	{
		fscanf(pFilePointer, "%s %[^\n]\n", &cActuatorFileName, buffer_line);
		actuatorFileName = cActuatorFileName;
		std::cout << actuatorFileName.c_str() << std::endl;
	}

	// Check Winch ID and joint coexistence
	if ((winchId !=0) && (this->GetType() == 3))
	{
		std::stringstream ss;
		ss << "Actuator and Winch boundary conditions defined at the same BCP --> BCP num: " << this->GetId()+1;
		throw ValueError(ss.str());
	}

	posWrtCdgLocal = pos;
	posG_BCP.rows(0,2) = pos;
	
}


void BCP::UpdateBoundary()
{
	double a = 0.0;
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// ANCHOR CLASS DEFINITION //////////////////////////
////////////////////////////////////////////////////////////////////////////
int AnchorBCP::GetType(void)
{
	return this->typeBcp;
}


void AnchorBCP::GetValues(double t)
{
	tBCP = t;
	vel = arma::zeros(3,1);
	acc = arma::zeros(3,1);
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// FAIRLEAD CLASS DEFINITION ////////////////////////
////////////////////////////////////////////////////////////////////////////
int FairleadBCP::GetType(void)
{
	return this->typeBcp;
}


void FairleadBCP::GetValues(double t)
{

	tBCP = t;
	ni = std::max(0,ni-10);
	do{
		ni = ni + 1;
	} while (tF(ni,0)<t);

	if(tF(ni,0)>t){
		ni = ni - 1;
	}

	dt = t - tF(ni,0);

	pos = (posF.row(ni) + dt * (posF.row(ni+1) - posF.row(ni)) / (tF(ni+1,0) - tF(ni,0))).t();
	vel = (velF.row(ni) + dt * (velF.row(ni+1) - velF.row(ni)) / (tF(ni+1,0) - tF(ni,0))).t();
	acc = (accF.row(ni) + dt * (accF.row(ni+1) - accF.row(ni)) / (tF(ni+1,0) - tF(ni,0))).t();

}


void FairleadBCP::Initialize(std::string folder_path)
{
	// Inicializo la variable donde guardar el numero de pasos temporales
	int nt;
	//Abro el fichero
 	std::string file_path = JoinPath(folder_path, actuatorFileName);
	std::ifstream datosPosF (file_path);

	// Leo el numero de pasos temporales a leer
	datosPosF >> nt;

	// Alocato la matriz que contiene la informacion
	tF = arma::zeros(nt,1);
	posF = arma::zeros(nt,3);
	velF = arma::zeros(nt,3);
	accF = arma::zeros(nt,3);
	// Leo toda la info
	for (int i=0; i<nt; i=i+1){
		datosPosF >> tF(i,0) >> posF(i,0) >> posF(i,1) >> posF(i,2) >> velF(i,0) >> velF(i,1) >> velF(i,2) >> accF(i,0) >> accF(i,1) >> accF(i,2);
	}
	//Cierro el fichero
	datosPosF.close();

	pos0 = pos;
	tBCP = -1.0;

	this->GetValues(0.0);
}


void FairleadBCP::ReadPropertiesASCII(FILE* &pFilePointer, std::string inputFilePath)
{	
	// Read properties from file
	BCP::ReadPropertiesASCII(pFilePointer);

	// Check if the actuator file exists
	std::string actuatorFilePath = JoinPath(inputFilePath, actuatorFileName);
	std::stringstream ss;
	ss << "ACTUATOR BCP NUMBER: " << this->GetId();
	CheckInputFile(actuatorFilePath, ss.str());
}


void FairleadBCP::Print(void)
{
	printf("BCP: %d PROPERTIES:\n");
	printf("--> PosX: %f - PosY: %f - PosZ: %f\n", this->posG_BCP[0], this->posG_BCP[1], this->posG_BCP[2]);
	printf("--> Winch Id: %d\n\n", this->winchId);
	printf("--> Actuator Filename: %s\n", actuatorFileName.c_str());

}


////////////////////////////////////////////////////////////////////////////
///////////////////////// JOINT CLASS DEFINITION ///////////////////////////
////////////////////////////////////////////////////////////////////////////
int JointBCP::GetType(void)
{
	return this->typeBcp;
}


void JointBCP::GetValues(double t)
{
	iLJ = 0;
	tBCP = t;
	pos = arma::mean(posLines).t();
	vel = arma::mean(velLines).t();
	acc = arma::mean(accLines).t();
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// BODYBCP CLASS DEFINITION /////////////////////////
////////////////////////////////////////////////////////////////////////////
int BodyBCP::GetType(void)
{
	return this->typeBcp;
}


void BodyBCP::GetValues(double t)
{
	tBCP = t;
	pos = posG_BCP.rows(0,2);
	vel = velG_BCP.rows(0,2);
	acc = accG_BCP.rows(0,2);
};

