
#ifndef body_hpp__
#define body_hpp__
#include <armadillo>
#include <string>
#include <cstdio>
#include "../BCPs/BCPs.hpp"
#include "../Hydro/HydroDatabase.hpp"

class Body {
private:
	int id; // Body index. It is an unique number assigned to each body in the simulation
	int velBufferSize=1e4; // Velocity Buffer size;

public:
	int velBufferCount=0; // Stores the positon of the last columun of the velocity buffer matrix filled.
	int numDofs; // Numero de grados de libertad que se consideran en el cuerpo
	int* pDofs; // Array de grados de libertad que se consideran en el cuerpo

	int numBcps; // Numero de BCPs en el cuerpo
	int* pIndexBcps;
	HydroDatabase* hydro; // Hydrodynamnics body associated to the body
	BCP** pBodyBcps; // Array de pointers a los puntos de condicion de contorno

	arma::mat pos = arma::zeros(6,1); // Posicion del cuerpo
	arma::mat vel = arma::zeros(6,1); // Velocidad del cuerpo
	arma::mat velBuffer = arma::zeros(6, velBufferSize); // Buffer in order to store the body velocities and calculate Duhamel's integral term
	arma::mat acc = arma::zeros(6,1); // Aceleracion del punto

	arma::mat inertia = arma::zeros(6,6); // matriz de inercia del cuerpo

	arma::mat rotMat = arma::zeros(3,3); // Matriz de rotación
	arma::mat invRotMat = arma::zeros(3,3); // Matriz de rotación

	arma::mat bcpForces = arma::zeros(6,1); // Fuerzas que los BCPs ejercen sobre el cuerpo, en la referencia del CDG

	Body(int n); // Inicializa un objeto de clase cuerpo dandole el indice
	void ComputeBcpForces(void); // Calcula el efecto de las fuerzas sobre los BCPs sobre su CDG
	int GetId(void); // Returns the body identification number
	void ReadPropertiesASCII(FILE* filePointer); // Leer datos de los cuerpos
	void StoreVelocities(arma::mat pos); // Store last step velocity into the velocity buffer matrix
	void UpdateBcps(void); // Actualiza valores del BCP
	void WriteOut(double t); // Escribir datos a fichero
};


#endif // body_hpp__