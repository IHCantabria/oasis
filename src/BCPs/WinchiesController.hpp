
#ifndef WINCHIE_CONTROLER_FLAG
#define WINCHIE_CONTROLER_FLAG
#include <armadillo>
#include <string>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <ctime>
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "Winchies.hpp"

class Simulation;

class WinchieController {
public:

	int k = 0; // contador
	int nWinchies; // Number of winchies
	int indBody; // Index of the body that is controlled with the winches
	Winchie** Winchies; // Pointers a los winchies

	Simulation* pSim; // Pointer to simulation instance

	double Ac,Bc,Cc,Dc; // Controller state-space model parameters
	arma::mat Kc = arma::zeros(3,1); // Lead controller gains
	arma::mat Ki = arma::zeros(3,1); // Integral time gains
	double Af,Bf,Cf; // Filter state-space model parameters
	double Ar,Br,Cr,Dr; // Reference state-space model parameters
	arma::mat ur = arma::zeros(3,1); // Vector of reference positions x(m), y(m), heading(rad)
	arma::mat yb = arma::zeros(3,1); // Vector of actual positions of the body
	arma::mat error; // Mat with timeseries with the error

	arma::mat xc = arma::zeros(3,1); // Variables de estado del controlador
	arma::mat yc = arma::zeros(3,1); // Output del controlador
	arma::mat xf = arma::zeros(3,1); // Variables de estado del filtro
	arma::mat yf = arma::zeros(3,1); // Output del filtro
	arma::mat xr = arma::zeros(3,1); // Variables de estado de la referencia
	arma::mat yr = arma::zeros(3,1); // Output de la referencia

	double Kw; // Ganancia de los winchies

	double T_max, T_min;
	arma::uvec ind_x_pos, ind_x_neg, ind_y_pos, ind_y_neg, ind_g_pos, ind_g_neg;
	arma::mat coef_x_pos, coef_x_neg, coef_y_pos, coef_y_neg, coef_g_pos, coef_g_neg;

	WinchieController(void);
	WinchieController(int n, Winchie** Ws, Simulation* pIncSim); // Inicializar el objeto de la clase winchie controler

	void ReadPropertiesASCII(FILE* file_pointer); // Leer inputs

	void controlWinchies(void); // Apply control
	
};


#endif
