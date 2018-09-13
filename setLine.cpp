/*
Libreria para iniciar la linea, menos set_nLine que se define en la clase
*/

//LIBRERIAS Y OTROS COMANDOS
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include "CatLine.h"


extern double PI;
extern int nLines;
extern double g;
extern double rhoW;
extern double t_max;
extern double dt;


/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCION DE CatLine PARA LEER datosMoorings.dat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void CatLine::leer_datosMoorings () {

	int ii, jj; 
	std::string Dummy; 
	const int nInored=20; // numero de lineas que se llen para cada nueva linea

	//Abro el fichero
	std::ifstream datosMoorings ("datosMoorings.dat");

	//Ignoro la primera linea del fichero, que contiene el numero de lineas de mooring a estudiar
	datosMoorings >> Dummy; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n'); // El ignore sirve para ignorar el texto de la linea

	//Ignoro las lineas que ya se han leido
	for(ii=1;ii<nLine;ii++){
		for(jj=1;jj<=nInored;jj++){
			datosMoorings >> Dummy; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	//Ignoro las tres primeras lineas, donde pone "New line"
	for(ii=1;ii<=3;ii++){
		datosMoorings >> Dummy; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	//Leo todo
	datosMoorings >> nNodos; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> L;    datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> rho0; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> d;    datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> EA;   datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> beta; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> CB;   datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> Cmn;  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> Cdn;  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> Cdt;  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> GK;   datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> GC;   datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> Gmu;  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> Gvc;  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> Dz;   datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> posFair[0]; datosMoorings >> posFair[1]; datosMoorings >> posFair[2];  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> posAnch[0]; datosMoorings >> posAnch[1]; datosMoorings >> posAnch[2];  datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');

	//Cierro el fichero
	datosMoorings.close();

	A=PI*d*d*0.25;
	dL=L/nNodos;
	pos = new double[3*2*nNodos];
	vel = new double[3*2*nNodos];
	acc = new double[3*2*nNodos];
}

/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCION DE CatLine PARA INICIAR LA LINEA CON EL METODO QS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void CatLine::initLine (void) {

	int jj;

	xF=sqrt(pow((posFair[0]-posAnch[0]),2)+pow((posFair[1]-posAnch[1]),2));
	zF=(posFair[2]-posAnch[2]);

	qs_GetTen();

}