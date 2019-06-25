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
#include "../os_tools.hpp"


////////////////////////////////////////////////////////////////////////////
/////////////////////////// BCP CLASS DEFINITION ///////////////////////////
////////////////////////////////////////////////////////////////////////////
BCP::BCP(int incId)
{
	id = incId;
}


int BCP::GetType(void)
{
	return this->typeBcp;
}


void BCP::Initialize()
{
	double a = 0.0;
}


void BCP::ReadPropertiesASCII(FILE* &pFilePointer)
{
	// Declare variables
	char buffer_line [1000];

	//Ignoro las tres primeras lineas, donde pone "New Body"
	for(int ii=0; ii<3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFilePointer);
	}

	//Read data
	//¡¡¡¡¡¡¡ IMPORTANT: FIRST THREEE LINES (CONTAINING DATA) OF EACH BCP ARE IGNORED. ITS FUNCTIONALITY IS DEPRECATED !!!!!
	fscanf(pFilePointer, "%[^\n]\n", buffer_line);
	fscanf(pFilePointer, "%[^\n]\n", buffer_line);
	fscanf(pFilePointer, "%[^\n]\n", buffer_line);
	//fscanf(pFilePointer, "%lf %[^\n]", &numLinesBcp, buffer_line);
	//pBcpLineIndex = new int[numLinesBcp];
	//pBcpLineNode = new int[numLinesBcp];
	fscanf(pFilePointer, "%lf %lf %lf %[^\n]\n", &pos(0, 0), &pos(1, 0), &pos(2, 0), buffer_line);
	if (this->GetType() == 2)
	{
		fscanf(pFilePointer, "%s %[^\n]\n", &actuatorFileName, buffer_line);
	}
	else
	{
		fscanf(pFilePointer, "%[\n]\n", buffer_line);
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
void AnchorBCP::GetValues(double t)
{
	vel = arma::zeros(3,1);
	acc = arma::zeros(3,1);
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// FAIRLEAD CLASS DEFINITION ////////////////////////
////////////////////////////////////////////////////////////////////////////
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

	this->GetValues(0.0);
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


////////////////////////////////////////////////////////////////////////////
///////////////////////// JOINT CLASS DEFINITION ///////////////////////////
////////////////////////////////////////////////////////////////////////////
void JointBCP::GetValues(double t)
{
	iLJ = 0;
	pos = arma::mean(posLines).t();
	vel = arma::mean(velLines).t();
	acc = arma::mean(accLines).t();
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// BODYBCP CLASS DEFINITION /////////////////////////
////////////////////////////////////////////////////////////////////////////
void BodyBCP::GetValues(double t)
{
	pos = posG_BCP.rows(0,2);
	vel = velG_BCP.rows(0,2);
	acc = accG_BCP.rows(0,2);
};

