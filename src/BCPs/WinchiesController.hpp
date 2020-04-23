
#ifndef WINCHIE_CONTROLER_FLAG
#define WINCHIE_CONTROLER_FLAG
#include <armadillo>
#include <string>
#include "Winchies.hpp"

class WinchieController {
public:

	int nWinchies; // Identificador del winchie
	Winchie** Winchies; // Pointer a la linea en la que esta el winchie

	double max_tau;
	double min_tau;
	double t_ini;

	void set_WinchieController(int n, Winchie** Ws){nWinchies=n; Winchies = Ws;} // Inicializar el objeto de la clase winchie controler

	void leer_datosWinchieController(void); // Leer inputs

	void controlWinchies(void); // Apply control
	
};


#endif
