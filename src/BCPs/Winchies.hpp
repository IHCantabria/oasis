
#ifndef winchiedef_hpp__
#define winchiedef_hpp__
#include <armadillo>
#include <cstdio>
#include <string>
#include "../Lines/Lines.hpp"

class Winchie {
public:

	int nWinchie; // Identificador del winchie

	Line * LineW; // Pointer a la linea en la que esta el winchie
	int nLine; // Indice de la linea en la que está el winchie
	int LineBCP; // Flag que indica en que extremo de la linea

	double inertia; // momento de inercia del winchie
	double radius; // radio del winchie
	double drag;

	double theta = 0.0; // Angulo de rotacion dle winchie
	double omega = 0.0; // Velocidad de rotacion del winchie
	double alpha = 0.0; // Aceleracion de rotación del winchie

	double tau = 0.0; // Momento que aplica el motor sobre el winchie

	Winchie(int n){nWinchie = n;} // Inicializar el objeto de la clase winchie
	void ReadPropertiesASCII(FILE* pFile); // Leer inputs
	void computeWinchie(void); // Obten las aceleraciones del winchie y el cambio de longitud de las lineas
	
};


#endif //winchiedef_hpp__