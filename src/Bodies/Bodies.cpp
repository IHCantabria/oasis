#include <armadillo>
#include <string>
#include "Bodies.hpp"

// Leer datos de los cuerpos
void Body::leer_datosBody(void){

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
