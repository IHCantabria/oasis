#include <armadillo>
#include <string>
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
	temp_inertia.load(arma::hdf5_name("input/datosBodies.dat","inertia"));

	inertia = temp_inertia.subcube(arma::span(nBody-1),arma::span::all,arma::span::all);

}

// Obten la matriz de rotación y pasasela a los BCPs, junto con la posicion, velocidad y aceleración
void Body::getRotMat(void){

	RotMat = arma::zeros(3,3);

	double cr = cos(pos(3,0)); double sr = sin(pos(3,0));
	double cp = cos(pos(4,0)); double sp = sin(pos(4,0));
	double cy = cos(pos(5,0)); double sy = sin(pos(5,0));

	RotMat(0,0) = cp*cy; RotMat(0,1) = cy*sp*sr-cr*sy; RotMat(0,2) = sr*sy+cr*cy*sp;
	RotMat(1,0) = cp*sy; RotMat(1,1) = cr*cy+sp*sr*sy; RotMat(1,2) = cr*sp*sy-cy*sr;
	RotMat(2,0) =   -sp; RotMat(2,1) =          cp*sr; RotMat(2,2) =          cp*cr;

	//[ cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)]
	//[ cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll)]
	//[         -sin(pitch),                               cos(pitch)*sin(roll),                               cos(pitch)*cos(roll)]

	for(int ii=0;ii<nBCPs;ii=ii+1){	
		BodyBCPs[ii].posG_body = pos; 
		BodyBCPs[ii].velG_body = vel; 
		BodyBCPs[ii].accG_body = acc; 
		BodyBCPs[ii].RotMat = RotMat;
	}

}
