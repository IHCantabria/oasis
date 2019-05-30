
#ifndef HYDRO_FLAG
#define HYDRO_FLAG
#include <armadillo>
#include <string>
#include "../Bodies/Bodies.hpp"

class Hydro {
public:

	int nBodies; // Numero de cuerpos
	Body * Bodies; // Array de pointers a los cuerpos

	arma::mat hydro; // matriz de hidrostatica de los cuerpos
	arma::mat inertia; // matriz de inercia de los cuerpos
	arma::mat Ainf; // Matriz de masas añadidas asintotica
	arma::mat invM; // inversa de la suma de la matriz de inercias y la matriz de masas añadidas asintotica

	int nAngles; // Numero de angulos
	arma::mat angles; // Angulos

	int nFreqs; // Numero de frecuencias
	arma::mat periods; // Periodos
	arma::mat frequencies; // Frecuencias

	arma::field<arma::mat> A; // Matriz de masas añadidas

	arma::field<arma::mat> B; // Matric de damping

	arma::field<arma::mat> F_mod; // Modulo del vector de fuerzas
	arma::field<arma::mat> F_phase; // Fase del vector de fuerzas
	arma::field<arma::mat> F_real; // Parte real del vector de fuerzas
	arma::field<arma::mat> F_imag; // Parte imaginaria del vector de fuerzas

	arma::field<arma::mat> QTFd_mod; // Modulo de la QTF de diferencia de frecuencias
	arma::field<arma::mat> QTFd_phase; // Fase de la QTF de diferencia de frecuencias
	arma::field<arma::mat> QTFd_real; // Parte real de la QTF de diferencia de frecuencias
	arma::field<arma::mat> QTFd_imag; // Parte imaginaria de la QTF de diferencia de frecuencias

	arma::field<arma::mat> QTFs_mod; // Modulo de la QTF de suma de frecuencias
	arma::field<arma::mat> QTFs_phase; // Fase de la QTF de suma de frecuencias
	arma::field<arma::mat> QTFs_real; // Parte real de la QTF de suma de frecuencias
	arma::field<arma::mat> QTFs_imag; // Parte imaginaria de la QTF de suma de frecuencias

	arma::mat drag; // Coeficientes de drag a corrientes

	arma::mat kl; // Coeficientes de fricción lineal para compensar ausencia de efectos viscosos
	arma::mat knl; // Coeficientes de fricción no lineal para compensar ausencia de efectos viscosos

	double t_relax; // Tiempo de relajación
	double dt_IRF; // Paso de tiempo con el que se guardan los datos
	int nt_IRF; // Numero de datos de velocidad que se guardan para convolucionar con IRF
	arma::mat vel_hist; // Velocidades guardadas en todos los dofs en los ultimos t_relax segundos
	arma::cube IRF; // Impulse response function

	double t_Fe; // Tiempo de reconstrucción de fuerzas de excitación
	double dt_Fe; // Paso temporal con el que se reconstruyen las fuerzas de excitacion
	int nt_Fe; // Numero de datos de fuerza de excitacion
	arma::mat Fe; // Fuerzas de excitación para el oleaje del caso

	double Hs; // Altura significativa
	double Tp; // Periodo pico

	int nComp; // Numero de componentes
	arma::mat wave_periods; // Periodos del espectro
	arma::mat wave_frequencies; // Frecuencias del espectro
	arma::mat wave_heights; // Alturas del espectro
	arma::mat wave_phases; // Fases del espectro

	arma::mat HydroForces;

	// Inicializar el objeto de la clase hydro dandole el numero de cuerpos y el vector de pointers a los cuerpos
	void set_Hydro(int n, Body * Bs){nBodies=n; Bodies = Bs;}

	void leer_datosHydro(void); // Leer inputs

	void computeIRF(void); // Calcula la impulse response function

	void computeWaveSpectrum(void); // Calcula el espectro del oleaje
	void computeFe(void); // Calcula la serie temporal de fuerzas de excitación

	void computeHydroForces(double t); // Obten las fuerzas hidroestaticas e hidrodinamicas en el tiempo deseado
	
};


#endif