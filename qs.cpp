/*
Libreria de método Quasi Static
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

void CatLine::qs_GetTen(void){

	double om = (rho0-rhoW*A);





	std::cout << "rhoW "  <<  rhoW << std::endl;
	std::cout << "nNodos " << this->nNodos << std::endl << std::endl;
}