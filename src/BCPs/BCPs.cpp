/*
Libreria para las condiciones de contorno
*/

//LIBRERIAS Y OTROS COMANDOS
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


/*
En el caso de BodyBCP, la funcion getValues puede hacer dos cosas distintas en función de como se llame.
	- Si se llama con t < 0.5, se actualiza la posición, velocidad y aceleración del BCP, y se toman
	  como nulas las fuerzas que actuan sobre este, para que los distintos elementos a los que esté
	  conectado el BCP añadan las fuerzas y momentos correspondientes.
	- Si se llama con t > 0.5, se calculan las fuerzas y momentos equivalentes en el centro de masas 
	  que provocan las fuerzas y momentos ejercidas sobre el BCP.
*/
void BodyBCP::getValues(double t){

	if (t < 0.5){
		posG =  RotMat*posL;
		posG_BCP.rows(0,2) = posG_body.rows(0,2) + posG;
		velG_BCP.rows(0,2) = velG_body.rows(0,2) + arma::cross(velG_body.rows(3,5),posG);
		accG_BCP.rows(0,2) = accG_body.rows(0,2) + arma::cross(accG_body.rows(3,5),posG);
		ForceBCP = 0.0;
	} else {
		double rx = posG(0,0);
		double ry = posG(1,0);
		double rz = posG(2,0);
		double Mx = ForceBCP(0,0);
		double My = ForceBCP(1,0);
		double Mz = ForceBCP(2,0);

		//[           Mz*rx, rx*ry, My*ry + Mz*rz]
		//[           Mz*ry,  ry^2,        -Mx*ry]
		//[ - Mx*rx - My*ry, ry*rz,        -Mx*rz] / (My*ry^2 + Mx*rx*ry + Mz*ry*rz)
		arma::mat matrix = arma::zeros(3,3);
		matrix(0,0) =           Mz*rx; matrix(0,1) = rx*ry; matrix(0,2) = My*ry + Mz*rz;
		matrix(1,0) =           Mz*ry; matrix(1,1) = ry*ry; matrix(1,2) =        -Mx*ry;
		matrix(2,0) = - Mx*rx - My*ry; matrix(2,1) = rx*rz; matrix(2,2) =        -Mx*rz;
		matrix = matrix/(My*ry*ry + Mx*rx*ry + Mz*ry*rz);

		arma::mat vector = ForceBCP.rows(3,5) - arma::dot(posG,ForceBCP.rows(3,5))/arma::dot(posG,posG);
		vector(1,0) = 0.0;

		ForceCDG.rows(0,2) = ForceBCP.rows(0,2) + matrix*vector;
		ForceCDG.rows(3,5) = arma::cross(posG,ForceBCP.rows(0,2)) + ForceBCP.rows(3,5);
	}

}