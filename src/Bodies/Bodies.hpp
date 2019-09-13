
#ifndef bodydef_hpp__
#define bodydef_hpp__

#include <armadillo>
#include <string>
#include <cstdio>
#include "../BCPs/BCPs.hpp"
#include "../Hydro/HydroDatabase.hpp"

class solver_data;
class Simulation;

class Body {
private:
	int id; // Body index. It is an unique number assigned to each body in the simulation

public:
	// General and management variables
	BCP** pBodyBcps; // Array de pointers a los puntos de condicion de contorno
	int* pDofs; // Array to store the number of degrees of freedom considered in the body
	int* pIndexBcps; // Array to store the index of the boundary condition points
	HydroDatabase* pHydro; // Hydrodynamnics body associated to the body
	int hydroDatabaseIndex; // Index of the body in the associated hydrodynamic database, if any
	std::string hydroDatabaseName; // Stores the hydrodynamic database name
	Simulation* pSim; // Pointer to simulation instance
	int takeCOGHydroDatabase; // Stores if read the initial COG position from the hydrodynamic database
	int numBcps=0; // Number of Boundary Condition Points considered
	int numDofs=0; // // Number of Degrees of freedom considered for the body
	
	// kinematic and Dynamic properties attributes
	arma::mat acc = arma::zeros(6,1); // Body acceleration w.r.t the global reference system
	arma::mat bcpForces = arma::zeros(6,1); // Fuerzas que los BCPs ejercen sobre el cuerpo, en la referencia del CDG
	arma::mat hydrostaticForces = arma::zeros(6, 1); // Storage for the hydrostatic forces
	arma::mat inertia = arma::zeros(6,6); // matriz de inercia del cuerpo
	arma::mat invRotMat = arma::zeros(3,3); // Matriz de rotación
	arma::mat pos = arma::zeros(6,1); // Body position w.r.t the global reference system
	arma::mat pos_init = arma::zeros(6,1); // Body initial position w.r.t global reference system. This is used as a reference for hydrostatic force calculation
	arma::mat radiationForces = arma::zeros(6, 1); // Storage for the wave radiation forces
	arma::mat rotMat = arma::zeros(3,3); // Matriz de rotación
	arma::mat vel = arma::zeros(6,1); // Body velocity w.r.t the global reference system
	arma::mat velBuffer = arma::zeros(6, velBufferSize); // Velocity Buffer (global coords) in order to store the body velocities and calculate Duhamel's integral term
	int velBufferSize=1e2; // Velocity Buffer size;
	int velBufferCount=0; // Stores the positon of the last columun of the velocity buffer matrix filled.

	FILE* pfile_DOF_1;
	FILE* pfile_DOF_2;
	FILE* pfile_DOF_3;
	FILE* pfile_DOF_4;
	FILE* pfile_DOF_5;
	FILE* pfile_DOF_6;
	FILE* pfile_HSF;
	FILE* pfile_WRF;

	// Methods definition
	Body(void){};
	Body(int n, Simulation* pSim); // Inicializa un objeto de clase cuerpo dandole el indice
	void ComputeBcpForces(void); // Calcula el efecto de las fuerzas sobre los BCPs sobre su CDG
	int GetId(void); // Returns the body identification number
	void LoadDependencies(void); // Load additional files and data necessary for the body type
	void LoadHydrodynamicDatabase(void); // Loads the hydrodynamics database associated, if any 
	void ReadPropertiesASCII(FILE* filePointer); // Leer datos de los cuerpos
	void OpenOutputFilesASCII(std::string path);
	void CloseOutputFilesASCII (void);
	void StoreVelocities(void); // Store last step velocity into the velocity buffer matrix
	void UpdateBcps(void); // Actualiza valores del BCP
	void UpdateHydrostaticForces(void); // Update the value of the wave radiation forces the current step
	void UpdateRadiationForces(void); // Update the value of the wave radiation forces the current step
	void WriteOut(double t); // Escribir datos a fichero
};

#endif // bodydef_hpp__