
#ifndef cumminsregulardef_hpp__
#define cumminsregulardef_hpp__

#include <armadillo>
#include "HydroDatabase.hpp"

class CumminsRegular: public HydroDatabase
{
public:
    // Declare class variables
    double waveHeight = 0.0;
    double wavePeriod = 0.0;
    double waveHeading = 0.0;
     
    // Declare class constructors
    CumminsRegular(int incId, int incIdBody, Body** incBodies, Simulation* pIncSim, double incWaveHeight, double incWavePeriod, double incWaveHeading);

    // Declare class methods
    arma::mat ComputeFirstWaveExcForce();
    inline arma::mat ComputeWaveExcForce();
    arma::mat CalculateHydrodynamicForces(double time);
    void InterpolateValues(void);

};

#endif // cumminsregulardef_hpp__