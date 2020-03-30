
#ifndef cumminsirregulardef_hpp__
#define cumminsirregulardef_hpp__

#include <armadillo>
#include "HydroDatabase.hpp"


class CumminsIrregular: public HydroDatabase
{
public:
    // Declare local variables
    double waveHeight = 0.0;
    double wavePeriod = 0.0;
    double waveHeading = 0.0;
     
    // Declare class constructors
    CumminsIrregular(int incId, Body** incBodies, Simulation* pIncSim, double incWaveHeight, double incWavePeriod, double incWaveHeading);

    // Declare class methods
    arma::mat ComputeFirstWaveExcForce(void);
    arma::mat ComputeHydrostaticForces(void);
    arma::mat CalculateHydrodynamicForces(double time);
	arma::mat ComputeRadiationForces(void);

};

#endif // cumminsirregulardef_hpp__