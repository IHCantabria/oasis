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


////////////////////////////////////////////////////////////////////////////
/////////////////////////// BCP CLASS DEFINITION ///////////////////////////
////////////////////////////////////////////////////////////////////////////
void BCP::BCP(int incId)
{
	id = incId;
}


int BCP::GetValues(void)
{
	return this->typeBCP;
}

void BCP::ReadPropertiesASCII(FILE* pFilePointer)
{
	// Declare variables
	char buffer_line [1000];

	//Ignoro las tres primeras lineas, donde pone "New Body"
	for(int ii=0; ii<3; ii++)
	{
		fscanf(pFilePointer, "%[^\n]", buffer_line);
	}

	//Read data
	//¡¡¡¡¡¡¡ IMPORTANT: FIRST THREEE LINES (CONTAINING DATA) OF EACH BCP ARE IGNORED. ITS FUNCTIONALITY IS DEPRECATED !!!!!
	fscanf(pFilePointer, "%[^\n]", buffer_line);
	fscanf(pFilePointer, "%[^\n]", buffer_line);
	fscanf(pFilePointer, "%[^\n]", buffer_line);

	datosBCPs >> nLinesBCP; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	BCPLineIndex = new int[nLinesBCP];
	BCPLineNode = new int[nLinesBCP];
	fscanf(pFilePointer, "%lf %lf %lf %[^\n]", &pos(0, 0), &pos(1, 0), &pos(2, 0), buffer_line);

	if (this->typeBcp == 2)
	{
		fscanf(pFilePointer, "%s %[^\n]", &actuatorFileName, buffer_line);
	}
	else
	{
		fscanf(pFilePointer, "%[\n]", buffer_line);
	}

	posL = pos;
	posG_BCP.rows(0,2) = pos;
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// ANCHOR CLASS DEFINITION //////////////////////////
////////////////////////////////////////////////////////////////////////////
int AnchorBCP::GetValues(void)
{
	return this->typeBCP;
}


void AnchorBCP::GetValues(double t){
	if (t<1e-12){
		typeBCP = 1;
	}
	vel = arma::zeros(3,1);
	acc = arma::zeros(3,1);
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// FAIRLEAD CLASS DEFINITION ////////////////////////
////////////////////////////////////////////////////////////////////////////
int FairleadBCP::GetValues(void)
{
	return this->typeBCP;
}


void FairleadBCP::GetValues(double t){

	tBCP = t;

	if (t<1e-12){
		typeBCP = 2;
		// Inicializo la variable donde guardar el numero de pasos temporales
		int nt;
		//Abro el fichero
		std::ifstream datosPosF (fileName);
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
	}

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
int JointBCP::GetValues(void)
{
	return this->typeBCP;
}


void JointBCP::GetValues(double t){
	if (t<1e-12){
		typeBCP = 3;
		iLJ = 0;
		posLines = arma::zeros(nLinesBCP,3);
		velLines = arma::zeros(nLinesBCP,3);
		accLines = arma::zeros(nLinesBCP,3);
	} else{
		iLJ = 0;
		pos = arma::mean(posLines).t();
		vel = arma::mean(velLines).t();
		acc = arma::mean(accLines).t();
	}
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// BODYBCP CLASS DEFINITION /////////////////////////
////////////////////////////////////////////////////////////////////////////
int BodyBCP::GetValues(void)
{
	return this->typeBCP;
}


void BodyBCP::GetValues(double t)
{
	if (t<1e-12){
		typeBCP = 4;
	}

	pos = posG_BCP.rows(0,2);
	vel = velG_BCP.rows(0,2);
	acc = accG_BCP.rows(0,2);

};

