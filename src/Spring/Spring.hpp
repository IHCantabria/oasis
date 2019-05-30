
#ifndef SPRING_FLAG
#define SPRING_FLAG
#include <armadillo>
#include <string>
#include "../BCPs/BCPs.hpp"

class Spring {
public:

	int nSpring; // Spring number

	int stressModelFlag; // Flag que determina que modelo de stress se da, [1: matriz lineal, 2: curva stress-strain]
	int dampingFlag; // Flag que determina si se considera dampng en los muelles, [0: no se considera, 1: se considera]

	int BCP_1; // Indices de los BCPs a los que esta conectado el muelle
	int BCP_2;
	BCP* SpringBCP[2]; // Puntos de contorno a los que está unido el muelle

	arma::field<arma::mat> SpringVectors; // Vectores unitarios que definen la orientación del muelle en local. orden normal (x), tangente 1 (y), tangente 2 (z)
	arma::cube SpringMatrix_K = arma::zeros(6,6,2); // Matriz de rigided del muelle en los dos extremos.
	arma::cube SpringMatrix_D = arma::zeros(6,6,2); // Matriz de damping del muelle en los dos extremos.
	// Para 6 grados de libertad y para los dos BCPs, una matriz 6x2 con el numero de datos, n, en las cuervas stress strain
	int n_StressStrain[6][2];
	// Para 6 grados de libertad y para los dos BCPs, una matriz nx1 con los estados de strain y otra nx6 con los stresses correspondientes
	arma::field<arma::mat> data_StressStrain;

	void set_nSpring(int n){nSpring=n;} // Inicializa un objeto de clase muelle dandole el indice
	void leer_datosSprings(void); // Lee inputs de los muelles
	void computeSpringForces(void); // Calcula las fuerzas que aplica el muelle en los BCPs y las guarda en variables de los BCPs
};



#endif