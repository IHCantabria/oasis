/*
MAIN DE LA IMPLEMENTACION DE NuevosFEM EN C++
*/

//LIBRERIAS Y OTROS COMANDOS
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include "CatLine.h"


double PI;
int nLines;
double g;
double rhoW;
double t;
double t_max;
double dt;


int main () {

	PI=acos(-1.0);
	t=0.0;

 	std::ifstream datosProblema ("datosProblema.dat");
	datosProblema >> g;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> rhoW;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> dt;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> t_max;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema.close();


	//LEO DE FICHERO Y PINTO EN PANTALLA CUANTAS LINEAS SE VAN A ESTUDIAR
	std::ifstream datosMoorings ("datosMoorings.dat");
	datosMoorings >> nLines;
	datosMoorings.close();
	std::cout << std::endl << "El numero de lineas a estudiar es: " << nLines << std::endl << std::endl;

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	CatLine * Line = new CatLine[nLines];

	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nLines; ii++){

		Line[ii].set_nLine(ii+1);
		Line[ii].leer_datosMoorings();

		std::cout << "Para la linea " << Line[ii].nLine << " , se ha leido:" << std::endl << std::endl;
		std::cout << "nNodos   " << Line[ii].nNodos << std::endl;
		std::cout << "L        " << Line[ii].L << std::endl;
		std::cout << "rho0     " << Line[ii].rho0 << std::endl;
		std::cout << "d        " << Line[ii].d << std::endl;
		std::cout << "EA       " << Line[ii].EA << std::endl;
		std::cout << "beta     " << Line[ii].beta << std::endl;
		std::cout << "CB       " << Line[ii].CB << std::endl;
		std::cout << "Cmn      " << Line[ii].Cmn << std::endl;
		std::cout << "Cdn      " << Line[ii].Cdn << std::endl;
		std::cout << "Cdt      " << Line[ii].Cdt << std::endl;
		std::cout << "GK       " << Line[ii].GK << std::endl;
		std::cout << "GC       " << Line[ii].GC << std::endl;
		std::cout << "Gmu      " << Line[ii].Gmu << std::endl;
		std::cout << "Gvc      " << Line[ii].Gvc << std::endl;
		std::cout << "Dz       " << Line[ii].Dz << std::endl;
		std::cout << "posFair  " << Line[ii].posFair[0] << " " << Line[ii].posFair[1] << " " << Line[ii].posFair[2] << std::endl;
		std::cout << "posAnch  " << Line[ii].posAnch[0] << " " << Line[ii].posAnch[1] << " " << Line[ii].posAnch[2] << std::endl << std::endl;

		Line[ii].initLine();
		Line[ii].write_out();
	}


	return 0;
}