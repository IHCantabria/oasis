
#ifndef SPRING_FLAG
#define SPRING_FLAG
#include <armadillo>
#include <string>
#include <yaml-cpp/yaml.h>
#include "../BCPs/BCPs.hpp"

class Spring
{
public:
	int nSpring; // Spring number

	int stressModelFlag; // Flag que determina que modelo de stress se da, [1: matriz lineal, 2: curva stress-strain]
	int dampingFlag;	 // Flag que determina si se considera damping en los muelles, [0: no se considera, 1: se considera]
	int frictionFlag;	 // Flag que determina si se considera friccion en los muelles, [0: no se considera, 1: se considera]
	int frameFlag;		 // Flag que indica que frame se utiliza para el calculo de las deformaciones y las fuerzas, [0: frame medio, 1: BCP1, 2: BCP2]

	int BCP_1; // Indices de los BCPs a los que esta conectado el muelle
	int BCP_2;
	int BCP_1_type; // Flags que eindican el tipo de BCP
	int BCP_2_type;
	BCP *SpringBCP[2]; // Puntos de contorno a los que está unido el muelle

	arma::field<arma::mat> SpringVectors;		  // Vectores unitarios que definen la orientación del muelle en local para cada BCP. orden normal (x), tangente 1 (y), tangente 2 (z)
	arma::mat SpringMatrix_K = arma::zeros(6, 6); // Matriz de rigided del muelle
	double mu_d, mu_s, vt, Dt;					  // Coeficiente de fricción dinamico, estatico y velocidad de friccion maxima
	int flagStickSlip = 0;						  // [0: Slip, 1: Stick]
	arma::mat SpringStrains_Stick;
	arma::mat a_1 = arma::zeros(4, 1), a_2 = arma::zeros(4, 1); // Variables para almacenar los coeficientes polinomicos del coeficiente de friccion
	arma::mat SpringMatrix_M = arma::zeros(6, 6);				// Matriz de friccion del muelle
	arma::mat SpringMatrix_D = arma::zeros(6, 1);				// Coeficientes de damping del muelle
	// Para 6 grados de libertad, el numero de datos, n, en las cuervas stress strain
	int n_StressStrain[6];
	// Para 6 grados de libertad una matriz nx1 con los estados de strain y otra nx6 con los stresses correspondientes
	arma::field<arma::mat> data_StressStrain;

	Spring(int n) { nSpring = n; }					// Inicializa un objeto de clase muelle dandole el indice
	void ReadPropertiesASCII(std::string filePath); // Lee inputs de los muelles
	void ReadPropertiesYAML(YAML::Node node);
	void computeSpringForces(void);					// Calcula las fuerzas que aplica el muelle en los BCPs y las guarda en variables de los BCPs
};

#endif