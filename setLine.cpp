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
#include <stdlib.h>
#include "CatLine.h"

extern double PI;
extern int nLines;
extern double g;
extern double rhoW;
extern double t;
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
	dL=L/(nNodos-1);

	pos = new double[3*nNodos];
	vel = new double[3*nNodos];
	acc = new double[3*nNodos];
	s = new double[nNodos];
	xc = new double[nNodos];
	zc = new double[nNodos];
	dxcds = new double[nNodos];
	dzcds = new double[nNodos];
	Te = new double[nNodos];


	for(ii=0;ii<nNodos;ii++) s[ii]=ii*dL;
}

/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCION DE CatLine PARA INICIAR LA LINEA CON EL METODO QS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void CatLine::initLine (void) {

	xF=sqrt(pow((posFair[0]-posAnch[0]),2)+pow((posFair[1]-posAnch[1]),2));
	zF=(posFair[2]-posAnch[2]);

	cosa=(posFair[0]-posAnch[0])/xF;
	sina=(posFair[1]-posAnch[1])/xF;

	this->qs_GetTen();
	this->qs_Solution();

	tenAnch[0] = cosa*HA;
	tenAnch[1] = sina*HA;
	tenAnch[2] = VA;

	tenFair[0] = cosa*HF;
	tenFair[1] = sina*HF;
	tenFair[2] = VF;

	for(int ii=0;ii<nNodos;ii++) {

		pos[3*ii] = posAnch[0]+cosa*xc[ii];
		pos[3*ii+1] = posAnch[1]+sina*xc[ii];
		pos[3*ii+2] = posAnch[2]+zc[ii];

		vel[3*ii] = 0.0;
		vel[3*ii+1] = 0.0;
		vel[3*ii+2] = 0.0;

	}

	acc = vel;

}


void CatLine::write_out (void) {

	int ii, nn1, nn2, nn3, nn4;

	char buffer1[50], buffer2[50], buffer3[50], buffer4[50];


	nn1=sprintf(buffer1,"NodePosX_%d.dat", nLine);
	std::ofstream xpos(buffer1);
		xpos << t << "    ";
		for(ii=0;ii<nNodos;ii++) xpos <<  posAnch[0]+cosa*xc[ii] << "    ";
		xpos << std::endl;
	xpos.close();

	nn2=sprintf(buffer2,"NodePosY_%d.dat", nLine);
	std::ofstream ypos(buffer2);
		ypos << t << "    ";
		for(ii=0;ii<nNodos;ii++) ypos <<  posAnch[1]+sina*xc[ii] << "    ";
		ypos << std::endl;
	ypos.close();

	nn3=sprintf(buffer3,"NodePosZ_%d.dat", nLine);
	std::ofstream zpos(buffer3);
		zpos << t << "    ";
		for(ii=0;ii<nNodos;ii++) zpos << posAnch[2]+zc[ii] << "    ";
		zpos << std::endl;
	zpos.close();

	nn4=sprintf(buffer4,"CatTen_%d.dat", nLine);
	std::ofstream ten(buffer4);
		ten << t << "    " << tenAnch[0] << "    " << tenAnch[1] << "    " << tenAnch[2] << "    "
		     << tenFair[0] << "    " << tenFair[1] << "    " << tenFair[2] << "    " << std::endl;
	ten.close();


}