#ifndef sinkingdef_hpp__
#define sinkingdef_hpp__

#include <armadillo>
#include <string>
#include <cstdio>
#include "../Hydro/HydroForce.hpp"

// Attribute class objects forward declaration
class Body;
class Simulation;

class Sinking {
public:

	double rhoW;

	Simulation* pSim; // Pointer to simulation instance
	Body* pSinkingBody; // Array to store the sinking body object

	int indBody; // index of the body in the general database
	int numGroups, numTimes;
	arma::mat groupsSizes; // [numGroups, 2]
	arma::mat groupsAreas; // [numGroups, 1]
	arma::mat groupsCenters; // [numGroups, 3]
	arma::mat groupsFillingTimes; // [numTimes, 1]
	arma::mat groupsFillingStates; // [numTimes, numGroups]

	int numHDBs;
	arma::mat InterpMasses; // Vector of filling masses for which the hydroforce objects are provided
	HydroDatabase** pHydro; // Vector of hydroforce objects for the different filling masses
	int indHydro1, indHydro2;
	double hydroInterpCoef;

	double totalFillingMass;
	double bodyReferenceMass;
	arma::mat bodyReferenceMassMat = arma::zeros(6,6);
	arma::mat groupsCOG = arma::zeros(3,1);
	arma::mat groupsInertia = arma::zeros(6,6);

	// Declare constructors
    Sinking(int n, Simulation* pSim_inp);

    // Methods definition
	void ReadPropertiesASCII(FILE* filePointer, std::string inputFolderPath); // Read sinking properties
	void UpdateGroupsFillingState(double t);
	void UpdateInterpHydro(void);
	void UpdateBodyProperties(void);
	void UpdateSinkingHydrodynamics(double t);
	void UpdateSinkingHydrostatics(double t);


};


#endif // sinkingdef_hpp__