#include <armadillo>
#include <string>
#include "../Bodies/Bodies.hpp"
#include "../Waves/Waves.hpp"

class Hydro {
public:

	int nBodies; // Numero de cuerpos
	Body * Bodies; // Array de pointers a los cuerpos

	arma::mat hydroS; // matriz de hidrostatica de los cuerpos
	arma::mat inertia; // matriz de inercia de los cuerpos

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

	arma::mat kl; // Coeficiontes de fricción lineal
	arma::mat knl; // Coeficiontes de fricción no lineal

	arma::mat Ainf; // Matriz de masas añadidas asintotica

	double t_relax; // Tiempo de relajación
	double dt; // Paso de tiempo con el que se guardan los datos
	int nt; // Numero de datos de velocidad que se guardan para convolucionar con IRF
	arma::mat vel_hist // Velocidades guardadas en todos los dofs en los ultimos t_relax segundos
	arma::mat IRF; // Impulse response function

	arma::mat Fe; // Fuerzas de excitación para el oleaje del caso

	Waves * Wave; // Pointer al oleaje del caso

	void leer_datosHydro(void); // Leer inputs

	void computeAinf(void); // Calcula la matriz de masa añadida asintotica
	void computeIRF(void); // Calcula la impulse response function
	void computeFe(void); // Calcula la serie temporal de fuerzas de excitación

	arma::mat computeHydroForce(double t); // Obten las fuerzas hidroestaticas e hidrodinamicas en el tiempo deseado
	
};