#include <armadillo>
#include <string>

class Waves {
public:

	double Hs; // Altura significativa
	double Tp; // Periodo pico

	int nComp; // Numero de componentes
	arma::mat periods; // Periodos del espectro
	arma::mat frequencies; // Frecuencias del espectro
	arma::mat heights; // Alturas del espectro
	arma::mat phases; // Fases del espectro

	void leer_datosWaves(void); // Lee inputs de oleaje
	void computeWaveSpectrum(void); // Calcula el espectro
	
};