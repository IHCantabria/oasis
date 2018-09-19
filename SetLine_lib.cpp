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
#include "classes.h"

extern double PI;
extern int nLines;
extern double g;
extern double rhoW;
extern double fondo;
extern double t;
extern double t_max;
extern double dt;


/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCION DE CatLine PARA LEER datosMoorings.dat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void MotherLine::leer_datosMoorings () {

	int ii, jj; 
	std::string Dummy; 
	const int nInored=20; // numero de lineas que se llen para cada nueva linea

	//Abro el fichero
	std::ifstream datosMoorings ("datosMoorings.dat");

	//Ignoro las tres primeras lineas del fichero, que contiene el numero de lineas a estudiar
	for(ii=1;ii<=3;ii++){
		datosMoorings >> Dummy; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea
	}

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

	if (posFair[2]<fondo || posAnch[2]<fondo) throw 0;

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


void MotherLine::print_out (void) {

		std::cout << "Para la linea " << this->nLine << " , se ha leido:" << std::endl << std::endl;
		std::cout << "nNodos   " << this->nNodos << std::endl;
		std::cout << "L        " << this->L << std::endl;
		std::cout << "rho0     " << this->rho0 << std::endl;
		std::cout << "d        " << this->d << std::endl;
		std::cout << "EA       " << this->EA << std::endl;
		std::cout << "beta     " << this->beta << std::endl;
		std::cout << "CB       " << this->CB << std::endl;
		std::cout << "Cmn      " << this->Cmn << std::endl;
		std::cout << "Cdn      " << this->Cdn << std::endl;
		std::cout << "Cdt      " << this->Cdt << std::endl;
		std::cout << "GK       " << this->GK << std::endl;
		std::cout << "GC       " << this->GC << std::endl;
		std::cout << "Gmu      " << this->Gmu << std::endl;
		std::cout << "Gvc      " << this->Gvc << std::endl;
		std::cout << "Dz       " << this->Dz << std::endl;
		std::cout << "posFair  " << this->posFair[0] << " " << this->posFair[1] << " " << this->posFair[2] << std::endl;
		std::cout << "posAnch  " << this->posAnch[0] << " " << this->posAnch[1] << " " << this->posAnch[2] << std::endl << std::endl;

	}

/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCION DE CatLine PARA INICIAR LA LINEA CON EL METODO QS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void TowingLine::initLine (void) {

	floor_flag = -1;

	xF=sqrt(pow((posFair[0]-posAnch[0]),2)+pow((posFair[1]-posAnch[1]),2));
	zF=(posFair[2]-posAnch[2]);

	cosa=(posFair[0]-posAnch[0])/xF;
	sina=(posFair[1]-posAnch[1])/xF;

	double Li= sqrt(pow(xF,2)+pow(zF,2));

	if ( Li>=L) {
		double cost=xF/Li;
		double sint=zF/Li;

		tenAnch[0] = 0;
		tenAnch[1] = 0;
		tenAnch[2] = 0;

		tenFair[0] = 0;
		tenFair[1] = 0;
		tenFair[2] = 0;

		for(int ii=0;ii<nNodos;ii++) {

			xc[ii]=cost*ii*Li/(nNodos-1);
			zc[ii]=sint*ii*Li/(nNodos-1);
		}

		for(int ii=0;ii<nNodos;ii++) {

			pos[3*ii] = posAnch[0]+cosa*xc[ii];
			pos[3*ii+1] = posAnch[1]+sina*xc[ii];
			pos[3*ii+2] = posAnch[2]+zc[ii];

		}


		std::cout << "WARNING: Towing line " << nLine << " tension is high. " << std::endl;

	}else{

		this->qs_GetTen();
		this->qs_Solution();

		if(isnan(xc[0])) throw 5;

		std::cout << std::endl << "xc    " << xc[0] << "  " << xc[1] << "  " << xc[2] << "  " << xc[3] << "  " << xc[4] << "  " << xc[5] << std::endl;
		std::cout << std::endl << "zc    " << zc[0] << "  " << zc[1] << "  " << zc[2] << "  " << zc[3] << "  " << zc[4] << "  " << zc[5] << std::endl << std::endl;


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

			if (pos[3*ii+2]<=fondo) throw 1;
		
		}
	}

	for(int ii=0;ii<nNodos;ii++) {
		vel[3*ii] = 0.0;
		vel[3*ii+1] = 0.0;
		vel[3*ii+2] = 0.0;
	}

	acc = vel;

}


void MooringLine::initLine (void) {

	int flag = 0;

	xF=sqrt(pow((posFair[0]-posAnch[0]),2)+pow((posFair[1]-posAnch[1]),2));
	zF=(posFair[2]-posAnch[2]);
	cosa=(posFair[0]-posAnch[0])/xF;
	sina=(posFair[1]-posAnch[1])/xF;

	double Li= sqrt(pow(xF,2)+pow(zF,2));

	if ( Li>=L) {
		double cost=xF/Li;
		double sint=zF/Li;

		tenAnch[0] = 0;
		tenAnch[1] = 0;
		tenAnch[2] = 0;

		tenFair[0] = 0;
		tenFair[1] = 0;
		tenFair[2] = 0;

		for(int ii=0;ii<nNodos;ii++) {

			xc[ii]=cost*ii*Li/(nNodos-1);
			zc[ii]=sint*ii*Li/(nNodos-1);
		}

		for(int ii=0;ii<nNodos;ii++) {

			pos[3*ii] = posAnch[0]+cosa*xc[ii];
			pos[3*ii+1] = posAnch[1]+sina*xc[ii];
			pos[3*ii+2] = posAnch[2]+zc[ii];

		}

		if ( Li==L && posAnch[2]==fondo && posFair[2]==fondo ) {
			std::cout << "WARNING: Mooring line " << nLine << " is laying on the floor " << std::endl;
		}else if (xF<1e-5){
			std::cout << "WARNING: Mooring line " << nLine << " is vertical. " << std::endl;
		}else{
			std::cout << "WARNING: Mooring line " << nLine << " tension is high. " << std::endl;
		}

	} else {

		if ( posAnch[2]==fondo && posFair[2]==fondo ) throw 2;
		if (xF<1e-5) throw 3;

		floor_flag = -1;

		this->qs_GetTen();
		this->qs_Solution();

		if(isnan(xc[0])) throw 5;

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
			if (pos[3*ii+2]<=fondo && ii>=1) flag = 1; 
		}

		if (flag == 1) {

			floor_flag = 1;

			this->qs_GetTen();
			this->qs_Solution();

			if(isnan(xc[0])) throw 5;


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
			}
		}
	}




	for(int ii=0;ii<nNodos;ii++) {
		vel[3*ii] = 0.0;
		vel[3*ii+1] = 0.0;
		vel[3*ii+2] = 0.0;
	}
	acc = vel;
}



void TensorLine::initLine (void) {

	xF=sqrt(pow((posFair[0]-posAnch[0]),2)+pow((posFair[1]-posAnch[1]),2));
	zF=(posFair[2]-posAnch[2]);

	double Li= sqrt(pow(xF,2)+pow(zF,2));

	if (Li<L) throw 4;

	cosa=(posFair[0]-posAnch[0])/xF;
	sina=(posFair[1]-posAnch[1])/xF;

	double cost=xF/Li;
	double sint=zF/Li;

	tenAnch[0] = 0;
	tenAnch[1] = 0;
	tenAnch[2] = 0;

	tenFair[0] = 0;
	tenFair[1] = 0;
	tenFair[2] = 0;

	for(int ii=0;ii<nNodos;ii++) {

		xc[ii]=cost*ii*Li/(nNodos-1);
		zc[ii]=sint*ii*Li/(nNodos-1);


		pos[3*ii] = posAnch[0]+cosa*xc[ii];
		pos[3*ii+1] = posAnch[1]+sina*xc[ii];
		pos[3*ii+2] = posAnch[2]+zc[ii];

		vel[3*ii] = 0.0;
		vel[3*ii+1] = 0.0;
		vel[3*ii+2] = 0.0;

	}

	acc = vel;

}


void MotherLine::write_out (void) {

	int ii, nn1, nn2, nn3, nn4;

	char buffer1[50], buffer2[50], buffer3[50], buffer4[50];


	nn1=sprintf(buffer1,"NodePosX_%d.txt", nLine);
	std::ofstream xpos(buffer1);
		xpos << t << "    ";
		for(ii=0;ii<nNodos;ii++) xpos <<  posAnch[0]+cosa*xc[ii] << "    ";
		xpos << std::endl;
	xpos.close();

	nn2=sprintf(buffer2,"NodePosY_%d.txt", nLine);
	std::ofstream ypos(buffer2);
		ypos << t << "    ";
		for(ii=0;ii<nNodos;ii++) ypos <<  posAnch[1]+sina*xc[ii] << "    ";
		ypos << std::endl;
	ypos.close();

	nn3=sprintf(buffer3,"NodePosZ_%d.txt", nLine);
	std::ofstream zpos(buffer3);
		zpos << t << "    ";
		for(ii=0;ii<nNodos;ii++) zpos << posAnch[2]+zc[ii] << "    ";
		zpos << std::endl;
	zpos.close();

	nn4=sprintf(buffer4,"CatTen_%d.txt", nLine);
	std::ofstream ten(buffer4);
		ten << t << "    " << tenAnch[0] << "    " << tenAnch[1] << "    " << tenAnch[2] << "    "
		     << tenFair[0] << "    " << tenFair[1] << "    " << tenFair[2] << "    " << std::endl;
	ten.close();


}