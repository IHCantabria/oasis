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
#include <cmath>
#include "classes.h"
#include <armadillo>



double PI;
int nLines, nMoorLines, nTowLines, nTenLines;
int nNodosTotal, nSistema;
double g;
double rhoW;
double fondo;
double t;
double t_max;
double dt;


int main () {

	int flag_read_eq, flag_write_eq;

	PI=acos(-1.0);
	t=0.0;

	nNodosTotal=0;

 	std::ifstream datosProblema ("datosProblema.dat");
	datosProblema >> g;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> rhoW;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> fondo;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> dt;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> t_max;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_read_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
	datosProblema >> flag_write_eq;    datosProblema.ignore(std::numeric_limits<int>::max(), '\n');
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
	for(int ii=0; ii<nMoorLines; ii=ii+1){
		MoorLine[ii].set_nLine(ii+1);
		try
		{
		MoorLine[ii].leer_datosMoorings();
		MoorLine[ii].print_out();
		nNodosTotal=nNodosTotal+MoorLine[ii].N;
		if(flag_read_eq==0) MoorLine[ii].initLine();
		MoorLine[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << MoorLine[ii].nLine << " is under the floor level." << std::endl << std::endl;
			if (e==2) std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl << std::endl;
			if (e==3) std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " is not tense and vertical. It should be pretensed. " << std::endl << std::endl;
			if (e==6) std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " touches the seafloor althoug none of its ends are there. " << std::endl << std::endl;
			if (e==5) {std::cout<< "ERROR: Mooring line " << MoorLine[ii].nLine << " initial shape can't be computed with QS method. " << std::endl;
			std::cout<< "Try reaching the desired initial condition with the dynamic method from a different initial condition " << std::endl << std::endl;}
			return 0;
		}

		MoorLine[ii].SEM_computeA();	
		
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	TowingLine * TowLine = new TowingLine[nLines];

	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nTowLines; ii=ii+1){

		TowLine[ii].set_nLine(ii+1+nMoorLines);
		try
		{
		TowLine[ii].leer_datosMoorings();
		TowLine[ii].print_out();
		nNodosTotal=nNodosTotal+TowLine[ii].N;
		if(flag_read_eq==0) TowLine[ii].initLine();
		TowLine[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << TowLine[ii].nLine << " is under the floor level. " << std::endl << std::endl;
			if (e==1) std::cout<< "ERROR: Towing line " << TowLine[ii].nLine << " touches the seafloor. " << std::endl << std::endl;
			if (e==5) {std::cout<< "ERROR: Towing line " << MoorLine[ii].nLine << " initial shape can't be computed with QS method. " << std::endl;
			std::cout<< "Try reaching the desired initial condition with the dynamic method from a different initial condition " << std::endl << std::endl;}
			return 0;
		}

		TowLine[ii].SEM_computeA();		
	}

	//ALOCATO UN VECTOR DE OBJETOS, UNO PARA CADA LINEA
	TensorLine * TenLine = new TensorLine[nLines];

	//INICIO LAS LINEAS Y PINTO RESULTADOS POR PANTALLA
	for(int ii=0; ii<nTenLines; ii=ii+1){

		TenLine[ii].set_nLine(ii+1+nMoorLines+nTowLines);
		try
		{
		TenLine[ii].leer_datosMoorings();
		TenLine[ii].print_out();
		nNodosTotal=nNodosTotal+TenLine[ii].N;
		if(flag_read_eq==0) TenLine[ii].initLine();
		TenLine[ii].SEM_getBaseFunctions();
		}
		catch (int e) 
		{
			if (e==0) std::cout<< "ERROR: Line " << TenLine[ii].nLine << " is under the floor level. " << std::endl << std::endl;
			if (e==4) std::cout<< "ERROR: Tensor line " << TenLine[ii].nLine << " is not tense. " << std::endl << std::endl;
			return 0;
		}

		TenLine[ii].SEM_computeA();	
		
	}

	std::cout << "nNodosTotal   " << nNodosTotal << std::endl; 

	nSistema=3*2*nNodosTotal;


	arma::mat y = arma::zeros(nSistema,1);
	arma::mat yprime = arma::zeros(nSistema,1);
	int ini=0;

	if(flag_read_eq==0){
		for(int ii=0; ii<nMoorLines; ii=ii+1){
			for(int jj=0; jj<MoorLine[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = MoorLine[ii].pos.row(jj).t();
					ini=ini+3;
			}
		}
		for(int ii=0; ii<nTowLines; ii=ii+1){
			for(int jj=0; jj<TowLine[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = TowLine[ii].pos.row(jj).t();
					ini=ini+3;
			}		
		}
		for(int ii=0; ii<nTenLines; ii=ii+1){
			for(int jj=0; jj<TenLine[ii].N; jj=jj+1){
					y.rows(ini,ini+2) = TenLine[ii].pos.row(jj).t();
					ini=ini+3;
			}		
		}

	}else if (flag_read_eq==1){
	 	std::ifstream equi("Equilibrio.dat");
			for(int ii=0;ii<nSistema;ii=ii+1) {
				equi >> y(ii,0);    equi.ignore(std::numeric_limits<int>::max(), '\n');
			}
		equi.close();

		ini=0;
		for(int ii=0; ii<nMoorLines; ii=ii+1){
			for(int jj=0; jj<MoorLine[ii].N; jj=jj+1){
				MoorLine[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}
		}
		for(int ii=0; ii<nTowLines; ii=ii+1){
			for(int jj=0; jj<TowLine[ii].N; jj=jj+1){
				TowLine[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}		
		}
		for(int ii=0; ii<nTenLines; ii=ii+1){
			for(int jj=0; jj<TenLine[ii].N; jj=jj+1){
				TenLine[ii].pos.row(jj)=y.rows(ini,ini+2).t();
				ini=ini+3;
			}		
		}
	}


	if (flag_write_eq == 1) {
		std::ofstream equi("Equilibrio.dat");
			for(int ii=0;ii<nSistema;ii=ii+1) equi << y(ii,0) << std::endl;
		equi.close();
	}

	for(int ii=0; ii<nMoorLines; ii=ii+1) MoorLine[ii].write_out();
	for(int ii=0; ii<nTowLines; ii=ii+1) TowLine[ii].write_out();
	for(int ii=0; ii<nTenLines; ii=ii+1) TenLine[ii].write_out();

	return 0;
}