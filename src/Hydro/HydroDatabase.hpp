
#ifndef hydrodatabasedef_hpp__
#define hydrodatabasedef_hpp__
#include <armadillo>
#include <string>
#include "HydroForce.hpp"

// Attribute class objects forward declaration
class Body;
class Simulation;
class solver_data;

class HydroDatabase: public HydroForce
{
private:
	int id;
	int activeDofs=6;
	
public:
	// General Management variables
	Body** pBodies; // Array de pointers a los cuerpos
	Simulation* pSim; // Pointer to the simulation instance. It gives fast access to the necessary simulation variables

	// Declare Auxiliar variables
	arma::cube* pTimeStartPos;
	double waveAmplitude=0.0;
	int numPeriodExc=0;
	int numHeadingExc=0;

	// Declare Hydrodynamic Storage Variables
	arma::mat cog; // Position of the center of gravity for the hydrodynamic Radiation-Diffraction problem
	int numPointsIRF; // Maximum number of points to describe the IRF function
	arma::cube** pAddedMass; // Matrix components: [body, dof, dof, freqs]
	arma::mat** pAddedMassHf; // Matrix components: [body, dof, dof]
	arma::mat** pAddedMassLf; // Matrix components: [body, dof, dof]
	arma::cube** pDampingRadiation; // Matrix components: [body, dof, dof, freqs];
	arma::mat** pDampingRadiationLf; // Matrix components: [body, dof, dof]
	arma::mat* fidd;
	arma::mat* pFrequencies; // Vector: [1, freqs]
	arma::mat* pHeadings; // Vector: [1, headings]
	arma::mat* pHydrostaticStiffness; // Matrix components: [dof, dof]
	arma::cube** pIRF; // Matrix components: [body, numTime, dof, dof]
	arma::mat** pIRFPoints; // Matrix components: [body, dof, dof]
	arma::mat IRFTime; // Matrix size containing the IRF time. Matrix dims: (1, numPointsIRF)
	arma::mat* pStructuralMass; // Matrix components: [dof, dof]
	arma::mat* pTotalMass; // Matrix components: [dof, dof]
	arma::mat* pTotalMass_inv; // Matrix components: [dof, dof]
	arma::cube* pMeanDrift; // Matrix components: [dofs, freqs, headings];
	int numBodies;
	int numFrequencies;
	int numHeadings;
	arma::cube*** pQtfDiff; // Matrix componets: [parts, dofs, freq1, freq2, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube*** pQtfSum; // Matrix componets: [parts, dofs, freq1, freq2, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube* pWaveExcitingMag; // Matrix components: [dofs, freqs, headings];
	arma::cube* pWaveExcitingPha; // Matrix components: [dofs, freqs, headings];

	// Class Constructors
	HydroDatabase(int incId, Body** incBodies, Simulation* pIncSim);

	// Class Methods
	arma::mat ComputeHydrostaticForces(void);
	void ComputeIRF(void); // Calcula la impulse response function
	void ReadHydroMechanicsHDF5(std::string filePath); // Leer inputs
	int GetId(void);
	void Print(void);
	
	double CalculateHydrodynamicForces(double time) {};
    void Load(void) {};	
};


#endif // hydrodatabasedef_hpp__