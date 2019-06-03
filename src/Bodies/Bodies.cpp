
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include "Bodies.hpp"

// Leer datos de los cuerpos
void Body::leer_datosBody(void){

	int ii, jj, kk;
	std::string Dummy;
	const int nInored=8; // numero de lineas que se leen para cada nueva linea

	//Abro el fichero
	std::ifstream datosBodies ("input/datosBodies.dat");

	//Ignoro la primera linea del fichero, que contiene el numero de lineas a estudiar
	datosBodies >> Dummy; datosBodies.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea

	//Ignoro las lineas que ya se han leido
	for(ii=1;ii<nBody;ii=ii+1){
		for(jj=1;jj<=nInored;jj=jj+1){
			datosBodies >> Dummy; datosBodies.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	//Ignoro las tres primeras lineas, donde pone "New line"
	for(ii=1;ii<=3;ii=ii+1){
		datosBodies >> Dummy; datosBodies.ignore(std::numeric_limits<int>::max(), '\n');
	}

	//Leo todo
	datosBodies >> nDOFs; datosBodies.ignore(std::numeric_limits<int>::max(), '\n');
	DOFs = new int[nDOFs];
	for(ii=0;ii<nDOFs;ii=ii+1){
		datosBodies >> DOFs[ii]; 
	}
	datosBodies.ignore(std::numeric_limits<int>::max(), '\n');

	datosBodies >> nBCPs; datosBodies.ignore(std::numeric_limits<int>::max(), '\n');
	index_BCPs = new int[nBCPs];
	BodyBCPs = new BCP*[nBCPs];
	for(ii=0;ii<nBCPs;ii=ii+1){
		datosBodies >> index_BCPs[ii]; 
	}
	datosBodies.ignore(std::numeric_limits<int>::max(), '\n');

	datosBodies >> pos(0,0); datosBodies >> pos(1,0); datosBodies >> pos(2,0); 
	datosBodies >> pos(3,0); datosBodies >> pos(4,0); datosBodies >> pos(5,0); 
	datosBodies.ignore(std::numeric_limits<int>::max(), '\n');

	//Cierro el fichero
	datosBodies.close();


	arma::cube temp_inertia;
	temp_inertia.load(arma::hdf5_name("input/flotante.h5","inertia"));

	inertia = temp_inertia.subcube(arma::span(nBody-1),arma::span::all,arma::span::all);

}

// Actualiza valores del BCP
void Body::updateBCPs(void){

	RotMat = arma::zeros(3,3); // Inicio la matriz de rotacion

	// Datos necesarios para la matriz de rotacion
	double cr = cos(pos(3,0)); double sr = sin(pos(3,0));
	double cp = cos(pos(4,0)); double sp = sin(pos(4,0));
	double cy = cos(pos(5,0)); double sy = sin(pos(5,0));

	// Matriz de rotacion calculada como Mz*My*Mx
	//[ cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)]
	//[ cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll)]
	//[         -sin(pitch),                               cos(pitch)*sin(roll),                               cos(pitch)*cos(roll)]
	RotMat(0,0) = cp*cy; RotMat(0,1) = cy*sp*sr-cr*sy; RotMat(0,2) = sr*sy+cr*cy*sp;
	RotMat(1,0) = cp*sy; RotMat(1,1) = cr*cy+sp*sr*sy; RotMat(1,2) = cr*sp*sy-cy*sr;
	RotMat(2,0) =   -sp; RotMat(2,1) =          cp*sr; RotMat(2,2) =          cp*cr;

	// Inversa de la matriz de rotacion
	invRotMat = arma::solve(RotMat,arma::eye(3,3));

	// Variable temporal
	arma::mat posG_temp;

	for(int ii=0;ii<nBCPs;ii=ii+1){	// Bucle sobre todos los BCPs
		posG_temp =  RotMat*(BodyBCPs[ii]->posL); // Brazo cdg-bcp en global
		BodyBCPs[ii]->posG = posG_temp; // Brazo cdg-bcp en global
		BodyBCPs[ii]->RotMat = RotMat; // Matriz de rotacion
		// Posicion del BCP en global, lo mismo para vel y acc.
		BodyBCPs[ii]->posG_BCP.rows(0,2) = pos.rows(0,2) + posG_temp; // Posicion del BCP en global
		BodyBCPs[ii]->velG_BCP.rows(0,2) = vel.rows(0,2) + arma::cross(vel.rows(3,5),posG_temp);
		BodyBCPs[ii]->accG_BCP.rows(0,2) = acc.rows(0,2) + arma::cross(acc.rows(3,5),posG_temp);
		BodyBCPs[ii]->posG_BCP.rows(3,5) = pos.rows(3,5);
		BodyBCPs[ii]->velG_BCP.rows(3,5) = vel.rows(3,5);
		BodyBCPs[ii]->accG_BCP.rows(3,5) = acc.rows(3,5);
		// Reseteo a cero la fuerza sobre el BCP
		BodyBCPs[ii]->ForceBCP = arma::zeros(6,1);
	}

	// Reseteo a cero la fuerza total de todos los BCPs
	BCPForces = arma::zeros(6,1);
}

// Calcula el efecto de las fuerzas sobre los BCPs sobre su CDG
void Body::computeBCPForces(void){

	// Variables locales necesarias usadas a continuacion
	arma::mat posG_temp;
	arma::mat ForceBCP_temp;
	arma::mat F_M;
	arma::mat M_F;
	double M_norm, det;
	double rx, ry, rz, Mx, My, Mz;

	for(int ii=0;ii<nBCPs;ii=ii+1){ // Bucle sobre los BCPs

		posG_temp = BodyBCPs[ii]->posG; // Guardo en variable temporal la posicion global del BCP
		ForceBCP_temp = BodyBCPs[ii]->ForceBCP; // Guardo en variable temporal las fuerzas y momentos sobre el BCP

		F_M = arma::zeros(3,1); // Inicio a cero al fuerza en el cdg causada por el momento sobre el BCP 
		M_norm = arma::norm(ForceBCP_temp.rows(3,5)); // Calculo la norma del momento sobre el BCP

		// Si hay momento en el BCP, calculo la fuerza que produce sobre le cdg
		if (M_norm > 1e-8) {

			// Datos necesarios usados a continuación
			rx = posG_temp(0,0); ry = posG_temp(1,0); rz = posG_temp(2,0);
			Mx = ForceBCP_temp(3,0); My = ForceBCP_temp(4,0); Mz = ForceBCP_temp(5,0);

			det = My*ry*ry + Mx*rx*ry + Mz*ry*rz;

			if (det > 1e-8) {
				// Inversa de la matriz del sistema para resover la ecuacion M = rxF, donde F es desconocido
				//[           Mz*rx, rx*ry, My*ry + Mz*rz]
				//[           Mz*ry,  ry^2,        -Mx*ry]
				//[ - Mx*rx - My*ry, ry*rz,        -Mx*rz] / (My*ry^2 + Mx*rx*ry + Mz*ry*rz)
				arma::mat matrix = arma::zeros(3,3);
				matrix(0,0) =           Mz*rx; matrix(0,1) = rx*ry; matrix(0,2) = My*ry + Mz*rz;
				matrix(1,0) =           Mz*ry; matrix(1,1) = ry*ry; matrix(1,2) =        -Mx*ry;
				matrix(2,0) = - Mx*rx - My*ry; matrix(2,1) = rx*rz; matrix(2,2) =        -Mx*rz;
				matrix = matrix/det;

				// Vector del sistema para resover la ecuacion M = rxF, donde F es desconocido (M - Mt)
				arma::mat vector = ForceBCP_temp.rows(3,5) - arma::dot(posG_temp,ForceBCP_temp.rows(3,5))/arma::dot(posG_temp,posG_temp);
				vector(1,0) = 0.0;

				// Obtengo la fuerza en el cdg causada por el momento en el bcp, ya en global, 
				// por que tanto el brazo posG, como la fuerza ForceBCP, estan en global.
				F_M = matrix*vector;
			}
		}

		M_F = arma::cross(posG_temp,ForceBCP_temp.rows(0,2)); // Momento sobre el cdg causado por la fuerza en el bcp, en global

		BCPForces.rows(0,2) = BCPForces.rows(0,2) + ForceBCP_temp.rows(0,2) + F_M; // Acumulo la fuerza total sobre el cdg en global.
		BCPForces.rows(3,5) = BCPForces.rows(3,5) + invRotMat * (ForceBCP_temp.rows(3,5) + M_F); // Acumulo el momento total sobre el cdg en local.
	}

}

// Escibir datos a fichero
void Body::write_out(double t){

	int ii, nn1, nn2, nn3, nn4, nn5, nn6;

	char buffer1[50], buffer2[50], buffer3[50], buffer4[50], buffer5[50], buffer6[50];

	if (t<1e-12){
		nn1=sprintf(buffer1,"output/DOF_1_Body_%d.txt", nBody);
		std::ofstream xpos(buffer1);
			xpos << t << "    " <<  this->pos(0,0) << "    "<<  this->vel(0,0) << "    "<<  this->acc(0,0) << std::endl;
		xpos.close();

		nn2=sprintf(buffer2,"output/DOF_2_Body_%d.txt", nBody);
		std::ofstream ypos(buffer2);
			ypos << t << "    " <<  this->pos(1,0) << "    "<<  this->vel(1,0) << "    "<<  this->acc(4,0) << std::endl;
		ypos.close();

		nn3=sprintf(buffer3,"output/DOF_3_Body_%d.txt", nBody);
		std::ofstream zpos(buffer3);
			zpos << t << "    " <<  this->pos(2,0) << "    "<<  this->vel(2,0) << "    "<<  this->acc(2,0) << std::endl;
		zpos.close();

		nn4=sprintf(buffer4,"output/DOF_4_Body_%d.txt", nBody);
		std::ofstream ropos(buffer4);
			ropos << t << "    " <<  this->pos(3,0) << "    "<<  this->vel(3,0) << "    "<<  this->acc(3,0) << std::endl;
		ropos.close();

		nn5=sprintf(buffer5,"output/DOF_5_Body_%d.txt", nBody);
		std::ofstream pipos(buffer5);
			pipos << t << "    " <<  this->pos(4,0) << "    "<<  this->vel(4,0) << "    "<<  this->acc(4,0) << std::endl;
		pipos.close();

		nn6=sprintf(buffer6,"output/DOF_6_Body_%d.txt", nBody);
		std::ofstream yapos(buffer6);
			yapos << t << "    " <<  this->pos(5,0) << "    "<<  this->vel(5,0) << "    "<<  this->acc(5,0) << std::endl;
		yapos.close();

	}else{
		nn1=sprintf(buffer1,"output/DOF_1_Body_%d.txt", nBody);
		std::ofstream xpos;
			xpos.open(buffer1, std::ios_base::app);
			xpos << t << "    " <<  this->pos(0,0) << "    "<<  this->vel(0,0) << "    "<<  this->acc(0,0) << std::endl;
		xpos.close();

		nn2=sprintf(buffer2,"output/DOF_2_Body_%d.txt", nBody);
		std::ofstream ypos;
			ypos.open(buffer2, std::ios_base::app);
			ypos << t << "    " <<  this->pos(1,0) << "    "<<  this->vel(1,0) << "    "<<  this->acc(1,0) << std::endl;
		ypos.close();

		nn3=sprintf(buffer3,"output/DOF_3_Body_%d.txt", nBody);
		std::ofstream zpos;
			zpos.open(buffer3, std::ios_base::app);
			zpos << t << "    " <<  this->pos(2,0) << "    "<<  this->vel(2,0) << "    "<<  this->acc(2,0) << std::endl;
		zpos.close();

		nn4=sprintf(buffer4,"output/DOF_4_Body_%d.txt", nBody);
		std::ofstream ropos;
			ropos.open(buffer4, std::ios_base::app);
			ropos << t << "    " <<  this->pos(3,0) << "    "<<  this->vel(3,0) << "    "<<  this->acc(3,0) << std::endl;
		ropos.close();

		nn5=sprintf(buffer5,"output/DOF_5_Body_%d.txt", nBody);
		std::ofstream pipos;
			pipos.open(buffer5, std::ios_base::app);;
			pipos << t << "    " <<  this->pos(4,0) << "    "<<  this->vel(4,0) << "    "<<  this->acc(4,0) << std::endl;
		pipos.close();

		nn6=sprintf(buffer6,"output/DOF_6_Body_%d.txt", nBody);
		std::ofstream yapos;
			yapos.open(buffer6, std::ios_base::app);
			yapos << t << "    " <<  this->pos(5,0) << "    "<<  this->vel(5,0) << "    "<<  this->acc(5,0) << std::endl;
		yapos.close();
	}

}
