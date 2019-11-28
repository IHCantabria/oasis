
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
    CumminsRegular();

    // Declare class methods
    inline arma::mat ComputeWaveExcForce();
    double CalculateHydrodynamicForces(double time);
	arma::mat ComputeRadiationForces();
    void InterpolateValues(void);

};

#endif // cumminsregulardef_hpp__