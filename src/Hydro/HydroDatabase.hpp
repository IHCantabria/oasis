
#ifndef hydrodatabasedef_hpp__
#define hydrodatabasedef_hpp__
#include <armadillo>
#include <string>
#include "HydroForce.hpp"
#include "Morison.hpp"

// Attribute class objects forward declaration
class Body;
class Simulation;
class solver_data;

class HydroDatabase : public HydroForce
{
private:
	int id;
	int activeDofs = 6;

public:
	// General Management variables
	int idBody;		  // Index of the current body of study
	Body **pBodies;	  // Array ob bodies pointers
	Simulation *pSim; // Pointer to the simulation instance. It gives fast access to the necessary simulation variables
	Morison *pMor;	  // Pointer to the Morison forces class

	// Declare Auxiliar variables
	arma::cube *pTimeStartPos;
	double waveAmplitude = 0.0;
	int numPeriodExc = 0;
	int numHeadingExc = 0;
	// Flag to indicate the hydro database type
	int hydroDatabaseFlag = 0; // [0: EHYDB, 1: H5]

	// Declare Hydrodynamic Storage Variables
	arma::mat cog;					 // Position of the center of gravity for the hydrodynamic Radiation-Diffraction problem
	int numPointsIRF;				 // Maximum number of points to describe the IRF function
	arma::cube **pAddedMass;		 // Matrix components: [body, dof, dof, freqs]
	arma::mat **pAddedMassHf;		 // Matrix components: [body, dof, dof]
	arma::mat **pAddedMassLf;		 // Matrix components: [body, dof, dof]
	arma::cube **pDampingRadiation;	 // Matrix components: [body, dof, dof, freqs];
	arma::mat **pDampingRadiationLf; // Matrix components: [body, dof, dof]
	arma::mat *fidd;
	arma::vec *pFrequencies;		  // Vector: [1, freqs]
	arma::vec *pHeadings;			  // Vector: [1, headings]
	arma::mat *pHydrostaticStiffness; // Matrix components: [dof, dof]
	arma::cube **pIRF;				  // Matrix components: [body, numTime, dof, dof]
	arma::mat **pIRFPoints;			  // Matrix components: [body, dof, dof]
	double IRFTotalTime;
	arma::mat IRFTime;			// Matrix size containing the IRF time. Matrix dims: (1, numPointsIRF)
	arma::mat *pStructuralMass; // Matrix components: [dof, dof]
	arma::mat *pTotalMass;		// Contains the StructuralMass + AddedMassHf: [dof, numDbBodies*dof]
	arma::cube *pMeanDrift;		// Matrix components: [dofs, freqs, headings];
	int numBodies;
	int numFrequencies;
	int numHeadings;
	arma::cube ***pQtfDiff;		  // Matrix components: [parts, dofs, freq1, freq2, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube ***pQtfSum;		  // Matrix components: [parts, dofs, freq1, freq2, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube *pWaveExcitingMag; // Matrix components: [dofs, freqs, headings];
	arma::cube *pWaveExcitingPha; // Matrix components: [dofs, freqs, headings];
	arma::cube *pWaveDiffMag;	  // Matrix components: [dofs, freqs, headings];
	arma::cube *pWaveDiffPha;	  // Matrix components: [dofs, freqs, headings];
	arma::cube *pWaveFKMag;		  // Matrix components: [dofs, freqs, headings];
	arma::cube *pWaveFKPha;		  // Matrix components: [dofs, freqs, headings];

	arma::cube ***QtfDiff_w;			 // Matrix components after SetUp: [parts, dofs, freq1_w, freq2_w, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube ***QtfSum_w;				 // Matrix components after SetUp: [parts, dofs, freq1_w, freq2_w, headings] | The first two are pointers to cube variables and last three are arma::cube
	arma::cube WE_Real_w;				 // Matrix components after SetUp: [headings, freqs_w, dofs];
	arma::cube WE_Imag_w;				 // Matrix components after SetUp: [headings, freqs_w, dofs];
	int numFrequencies_w;				 // Number of frequencies in the wave spectrum inside the frequency range of hydrodynamic database coefficients
	arma::vec freqs_w;					 // Wave frequencies vector: [1, numFrequencies_w] (Hz)
	arma::vec ang_freqs_w;				 // Wave angular frequencies vector: [1, numFrequencies_w] (rad/s)
	arma::mat kx_w;						 // Wave number in the x-axis in the frequency range of hydrodynamic database coefficients
	arma::mat ky_w;						 // Wave number in the y-axis in the frequency range of hydrodynamic database coefficients
	arma::field<arma::mat> amplitudes_w; // Wave amplitudes in the frequency range of hydrodynamic database coefficients
	arma::field<arma::mat> phases_w;	 // Wave phases in the frequency range of hydrodynamic database coefficients

	// Preprocess matrices for QTF computation
	arma::field<arma::mat> ampP; // Amplitude products matrices for all wave pieces
	arma::mat wD;				 // Angular frequencies differences matrices for all wave pieces
	arma::field<arma::mat> phD;	 // Phases differences matrices for all wave pieces
	arma::mat kxD;				 // X-axis wave number differences matrices for all wave pieces
	arma::mat kyD;				 // Y-axis wave number differences matrices for all wave pieces
	arma::mat wS;				 // Angular frequencies sums matrices for all wave pieces
	arma::field<arma::mat> phS;	 // Phases sums matrices for all wave pieces
	arma::mat kxS;				 // X-axis wave number sums matrices for all wave pieces
	arma::mat kyS;				 // Y-axis wave number sums matrices for all wave pieces

	// Mean drift force
	arma::mat F_meanDrift = arma::zeros(6, 1);

	// Variables for precomputed excitation forces
	arma::vec time_exc;				// Time vector for the precomputed excitation forces
	arma::mat force_excFirstOrder;	// Precomputed forces for the first order wave excitation
	arma::mat force_excSecondOrder; // Precomputed forces for the second order wave excitation

	// Class Constructors
	HydroDatabase(int incId, int incIdBody, Body **incBodies, Simulation *pIncSim);
	~HydroDatabase() {};

	// Class Methods
	arma::mat CalculateHydrostaticForces(double time);
	arma::mat ComputeRadiationForces(void);
	void ComputeIRF(std::string HDBname);				  // Compute the impulse response function
	void ComputeAsymptoticAddedMass(std::string HDBname); // Compute the impulse response function
	void ComputeTotalMass(void);						  // Compute the total mass matrix
	arma::mat GetCog(void);								  // Interface method, it returns center of gravity
	int GetNumBodies(void);								  // Returns the number of bodies in the database
	int GetNumPointsIrf(void);							  // Interface methods, it returns the number of points of the IRF
	arma::mat GetTotalMass(void);						  // Interface method, it returns the total mass matrix: StruturalMass+AddedMass+ViscousAddedMass
	void UpdateStructuralMass(arma::mat newStructuralMass);
	void UpdateTotalMass(void);
	arma::mat ComputeFirstWaveExcForce(double t); // Interpolate First Order Wave Excitation forces using first order polynomial
	arma::mat ComputeSecondWaveExcForce(double t);
	arma::mat ComputeMeanDrift(void);
	arma::mat CalculateHydrostaticPressure(double t);
	void Refresh(void); // This method refresh the state of the object properties
	int GetId(void);
	arma::mat GetInertiaMatrixInv(void);
	void Print(void);
	void SetUp(void);

	void LoadHydroDataEHYDB(std::string file_path); // Loads the corresponding hydrodynamic data in EHYDB format
	void LoadHydroDataH5(std::string file_path);	// Loads the corresponding hydrodynamic data in HDF5 format

	// Declare inherited virutal methods
	arma::mat CalculateHydrodynamicForces(double time);
	void LoadHydrodynamicData(std::string file_path); // Loads the corresponding hydrodynamic data
	void InterpolateHydro(HydroDatabase *pHydro1, HydroDatabase *pHydro2, double interpCoef);
	void UpdateHydroStiffness(arma::mat newHydrostaticStiffness);
};

#endif // hydrodatabasedef_hpp__