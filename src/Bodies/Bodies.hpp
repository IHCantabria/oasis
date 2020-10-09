
#ifndef bodydef_hpp__
#define bodydef_hpp__

#include <armadillo>
#include <string>
#include <cstdio>
#include "../BCPs/BCPs.hpp"
#include "../Hydro/HydroForce.hpp"

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
	HydroForce* pHydro; // Hydrodynamnics body associated to the body
	int hydroDatabaseIndex; // Index of the body in the associated hydrodynamic database, if any
	std::string hydroDatabaseName; // Stores the hydrodynamic database name
	int firstOrderExcitationFlag; // flag for first order excitation force
	int secondOrderExcitationFlag; // flag for second order excitation force
	arma::mat A_visc = arma::zeros(6,1); // viscous added mass coefficients
	arma::mat B_visc = arma::zeros(6,1); // viscous linear damping coefficients
	arma::mat B_visc2 = arma::zeros(6,1); // viscous cuadratic damping coefficients

	//Total hidrodynamic force
	arma::mat Fb = arma::zeros(6,1);


	Simulation* pSim; // Pointer to simulation instance
	int takeCOGHydroDatabase; // Stores if read the initial COG position from the hydrodynamic database
	int numBcps=0; // Number of Boundary Condition Points considered
	int numDofs=0; // // Number of Degrees of freedom considered for the body

	arma::span sysMatSpan1, sysMatSpan2;
	
	// kinematic and Dynamic properties attributes
	int velBufferSize=0; // Velocity Buffer size;
	int velBufferCount=0; // Stores the positon of the last columun of the velocity buffer matrix filled.
	double filling_mass = 0; // Body filling mass for sinking
	arma::mat acc = arma::zeros(6,1); // Body acceleration w.r.t the global reference system
	arma::mat bcpForces = arma::zeros(6,1); // Fuerzas que los BCPs ejercen sobre el cuerpo, en la referencia del CDG
	arma::mat hydrostaticForces = arma::zeros(6, 1); // Storage for the hydrostatic forces
	arma::mat inertia = arma::zeros(6,6); // matriz de inercia del cuerpo
	arma::mat invRotMat = arma::zeros(3,3); // Matriz de rotación
	arma::mat pos = arma::zeros(6,1); // Body's COG position w.r.t the global reference system
	arma::mat pos_eq = arma::zeros(6,1); // This is used as a equilibrium reference for hydrostatic force calculation
	arma::mat pos_filling_cog = arma::zeros(3,1); // Center of gravity position of the body filling for sinking
	arma::mat radiationForces = arma::zeros(6, 1); // Storage for the wave radiation forces
	arma::mat excitationForces_1 = arma::zeros(6, 1); // Storage for the wave excitation forces - 1st order
	arma::mat excitationForces_2 = arma::zeros(6, 1); // Storage for the wave excitation forces - 2nd order
	arma::mat rotMat = arma::zeros(3,3); // Matriz de rotación
	arma::mat vel = arma::zeros(6,1); // Body velocity w.r.t the global reference system
	arma::mat velBuffer; // Velocity Buffer (global coords) in order to store the body velocities and calculate Duhamel's integral term


	FILE* pfile_DOF_1;
	FILE* pfile_DOF_2;
	FILE* pfile_DOF_3;
	FILE* pfile_DOF_4;
	FILE* pfile_DOF_5;
	FILE* pfile_DOF_6;
	FILE* pfile_HSF;
	FILE* pfile_WRF;
	FILE* pfile_BCPF;
	FILE* pfile_WEF;

	// Methods definition
	Body(void){};
	Body(int n, Simulation* pSim); // Inicializa un objeto de clase cuerpo dandole el indice
	void ComputeBcpForces(void); // Calcula el efecto de las fuerzas sobre los BCPs sobre su CDG
	int GetId(void); // Returns the body identification number
	void LoadHydrodynamicDatabase(Body** hydroDatabaseBodies); // Loads the hydrodynamics database associated, if any 
	void ReadPropertiesASCII(FILE* filePointer); // Leer datos de los cuerpos
	void OpenOutputFilesASCII(std::string path);
	void CloseOutputFilesASCII (void);
	void StoreVelocities(bool restoreMatrix); // Store last step velocity into the velocity buffer matrix
	void UpdateBcps(void); // Actualiza valores del BCP
	void ResetBcps(void); // Resetea fuerzas BCPs
	// void UpdateHydrostaticForces(void); // Update the value of the wave radiation forces the current step
	// void UpdateRadiationForces(void); // Update the value of the wave radiation forces the current step
	void WriteOut(double t); // Escribir datos a fichero
};

#endif // bodydef_hpp__