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

	//Ignoro las dos primeras lineas del fichero, que contiene el numero de BCPs a estudiar
	for(ii=1;ii<=2;ii=ii+1){
		datosBCPs >> Dummy; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea
	}

	//Ignoro las lineas que ya se han leido
	for(ii=1;ii<nBCP;ii=ii+1){
		for(jj=1;jj<=nInored;jj=jj+1){
			datosBCPs >> Dummy; datosBCPs.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	//Ignoro las tres primeras lineas, donde pone "New line"
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
	//std::cout << "ANCHOR" << std::endl;
	vel = arma::zeros(3,1);
	acc = arma::zeros(3,1);
}

void FairleadBCP::getValues(double t){

	//std::cout << "FAIRLEAD" << std::endl;
	if (t<1e-12){
		// Inicializo la variable donde guardar el numero de pasos temporales
		int nt;
		//Abro el fichero
		std::ifstream datosPosF (fileName);
		// Leo el numero de pasos temporales a leer
		datosPosF >> nt;
		// Alocato la matriz que contiene la informacion
		posF = arma::zeros(nt,4);
		// Leo toda la info
		for (int i=0; i<nt; i=i+1){
			datosPosF >> posF(i,0) >> posF(i,1) >> posF(i,2) >> posF(i,3);
		}
		//Cierro el fichero
		datosPosF.close();

		x_spl = spline(posF.col(0),posF.col(1),1);
		y_spl = spline(posF.col(0),posF.col(2),1);
		z_spl = spline(posF.col(0),posF.col(3),1);
	}

	arma::mat sol_x = x_spl.spl_eval(t);
	arma::mat sol_y = y_spl.spl_eval(t);
	arma::mat sol_z = z_spl.spl_eval(t);

	pos(0,0) = sol_x(0);
	pos(1,0) = sol_y(0);
	pos(2,0) = sol_z(0);

	vel(0,0) = sol_x(1);
	vel(1,0) = sol_y(1);
	vel(2,0) = sol_z(1);

	acc(0,0) = sol_x(2);
	acc(1,0) = sol_y(2);
	acc(2,0) = sol_z(2);


}