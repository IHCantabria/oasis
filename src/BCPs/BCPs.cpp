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
#include <armadillo>
#include "BCPs.hpp"

void BCP::leer_datosBCPs(void){
	int ii, jj; 
	std::string Dummy;
	const int nInored=8; // numero de lineas de texto que se leen para cada nuevo BCP

	//Abro el fichero
	std::ifstream datosBCPs ("input/datosBCPs.dat");

	//Ignoro las cuatro primeras lineas del fichero, que contiene el numero de BCPs a estudiar
	for(ii=1;ii<=4;ii=ii+1){
		datosBCPs >> Dummy; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea
	}

	//Ignoro las lineas que ya se han leido
	for(ii=1;ii<nBCP;ii=ii+1){
		for(jj=1;jj<=nInored;jj=jj+1){
			datosBCPs >> Dummy; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	//Ignoro las tres primeras lineas, donde pone "New BCP"
	for(ii=1;ii<=3;ii=ii+1){
		datosBCPs >> Dummy; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	}

	//Leo todo
	datosBCPs >> nLinesBCP; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	BCPLineIndex = new int[nLinesBCP];
	BCPLineNode = new int[nLinesBCP];
	for(ii=0;ii<nLinesBCP;ii=ii+1){
		datosBCPs >> BCPLineIndex[ii];
	}
	datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	for(ii=0;ii<nLinesBCP;ii=ii+1){
		datosBCPs >> BCPLineNode[ii];
	}
	datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
	datosBCPs >> pos(0,0); datosBCPs >> pos(1,0); datosBCPs >> pos(2,0);  datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');

	datosBCPs >> fileName; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');

	//Cierro el fichero
	datosBCPs.close();

	posL = pos;
}

void AnchorBCP::getValues(double t){
	if (t<1e-12){
		typeBCP = 1;
	}
	vel = arma::zeros(3,1);
	acc = arma::zeros(3,1);
}

void FairleadBCP::getValues(double t){

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

void JointBCP::getValues(double t){
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


void BodyBCP::getValues(double t){

	pos = posG_BCP.rows(0,2);
	vel = velG_BCP.rows(0,2);
	acc = accG_BCP.rows(0,2);

};

