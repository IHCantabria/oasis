
#include <armadillo>
#include <string>
#include "Winchies.hpp"

// Leer inputs
void Winchie::leer_datosWinchies(void){

	int ii, jj, kk; 
	std::string Dummy;
	const int nInored=7; // numero de lineas que se leen para cada nuevo winchie

	//Abro el fichero
	std::ifstream datosWinchies ("input/datosWinchies.dat");

	//Ignoro la primera linea del fichero, que contiene el numero de winchies a estudiar
	datosWinchies >> Dummy; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea

	//Ignoro las lineas que ya se han leido
	for(ii=1;ii<nLine;ii=ii+1){
		for(jj=1;jj<=nInored;jj=jj+1){
			datosWinchies >> Dummy; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	//Ignoro las tres primeras lineas, donde pone "New winchie"
	for(ii=1;ii<=3;ii=ii+1){
		datosWinchies >> Dummy; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');
	}

	//Leo todo
	datosWinchies >> nLine; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');
	datosWinchies >> LineBCP; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');
	datosWinchies >> inertia; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');
	datosWinchies >> radius; datosWinchies.ignore(std::numeric_limits<int>::max(), '\n');

	//Cierro el fichero
	datosWinchies.close();

}

// Obten las aceleraciones del winchie y el cambio de longitud de las lineas
void Winchie::computeWinchie(void){

	double F;

	if (LineBCP == 1){
		F = arma::norm(LineW->ten_1);
	} else if (LineBCP == 2){
		F = arma::norm(LineW->ten_N);
	} else {
		std::cout << std::endl << "ERROR: Options are 1 or 2." << std::endl;
		throw std::exception();
	}

	alpha = (radius*F - tau)/inertia;

	LineW->dL = LineW->dL0 * (LineW->L + radius * theta) / LineW->L;

	//std::cout << "    theta =  " << theta << std::endl << std::endl;
	//std::cout << "    tau =  " << tau << std::endl << std::endl;

}

