
#ifndef cumminsirregulardef_hpp__
#define cumminsirregulardef_hpp__

#include <armadillo>
#include "HydroDatabase.hpp"


class CumminsIrregular: public HydroDatabase
{
public:
    // Declare local variables
     
    // Declare class constructors
    CumminsIrregular();

    // Declare class methods
    arma::mat ComputeFirstWaveExcForce(void);
    arma::mat ComputeHydrostaticForces(void);
    double CalculateHydrodynamicForces(double time);
	arma::mat ComputeRadiationForces(void);
    
};

#endif // cumminsirregulardef_hpp__