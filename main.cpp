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
#include "classes.h"


double PI;
int nLines, nMoorLines, nTowLines, nTenLines;
double g;
double rhoW;
double fondo;
double t;
double t_max;
double dt;


int main () {

	PI=acos(-1.0);
	t=0.0;

 	std::ifstream datosProblema ("datosProblema.dat");
	datosProblema >> g;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> rhoW;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> fondo;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> dt;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> t_max;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema.close();


	//LEO DE FICHERO Y PINTO EN PANTALLA CUANTAS LINEAS SE VAN A ESTUDIAR
	std::ifstream datosMoorings ("datosMoorings.dat");
	datosMoorings >> nMoorLines; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> nTowLines; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings >> nTenLines; datosMoorings.ignore(std::numeric_limits<int>::max(), '\n');
	datosMoorings.close();

	nLines=nMoorLines+nTowLines+nTenLines;

	std::cout << std::endl << "El numero de lineas a estudiar es: " << nLines << std::endl << std::endl;

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	MooringLine * MoorLine = new MooringLine[nLines];

	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nMoorLines; ii++){

		MoorLine[ii].set_nLine(ii+1);
		MoorLine[ii].leer_datosMoorings();
		MoorLine[ii].print_out();
		MoorLine[ii].initLine();
		MoorLine[ii].write_out();
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	TowingLine * TowLine = new TowingLine[nLines];

	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nTowLines; ii++){

		TowLine[ii].set_nLine(ii+1+nMoorLines);
		TowLine[ii].leer_datosMoorings();
		TowLine[ii].print_out();
		TowLine[ii].initLine();
		TowLine[ii].write_out();
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	TensorLine * TenLine = new TensorLine[nLines];

	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nTenLines; ii++){

		TenLine[ii].set_nLine(ii+1+nMoorLines+nTowLines);
		TenLine[ii].leer_datosMoorings();
		TenLine[ii].print_out();
		TenLine[ii].initLine();
		TenLine[ii].write_out();
	}


	return 0;
}