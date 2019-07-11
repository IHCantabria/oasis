
#ifndef hydrodatabasedef_hpp__
#define hydrodatabasedef_hpp__
#include <armadillo>
#include <string>

// Attribute class objects forward declaration
class Body;
class Simulation;
class solver_data;

class HydroDatabase
{
private:
	int id;
	int activeDofs=6;
	
public:
	int numPointsIRF; // Maximum number of points to describe the IRF function
	//int nBodies; // Numero de cuerpos
	Body** pBodies; // Array de pointers a los cuerpos
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
	arma::cube* pMeanDrift; // Matrix components: [dofs, freqs, headings];
	int numBodies;
	int numFrequencies;
	int numHeadings;
	arma::cube*** pQtfDiff; // Matrix componets: [parts, dofs, freq1, freq2, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube*** pQtfSum; // Matrix componets: [parts, dofs, freq1, freq2, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube* pWaveExcitingMag; // Matrix components: [dofs, freqs, headings];
	arma::cube* pWaveExcitingPha; // Matrix components: [dofs, freqs, headings];

	// Inicializar el objeto de la clase hydro dandole el numero de cuerpos y el vector de pointers a los cuerpos
	//void set_Hydro(int n, Body** Bs){nBodies=n; Bodies = Bs;}
	HydroDatabase(int incId, Body** incBodies);
	arma::mat ComputeHydrostaticForces(void);
	arma::mat ComputeRadiationForces(double t, solver_data SD);
	arma::mat ComputeFirstWaveExcForce(double t);
	void ComputeIRF(void); // Calcula la impulse response function
	void ReadHydroMechanicsHDF5(std::string filePath); // Leer inputs
	int GetId(void);
	void Print();
	
	//void computeWaveSpectrum(void); // Calcula el espectro del oleaje
	//void computeFe(void); // Calcula la serie temporal de fuerzas de excitación

	//void computeHydroForces(double t); // Obten las fuerzas hidroestaticas e hidrodinamicas en el tiempo deseado
	
};


#endif // hydrodatabasedef_hpp__