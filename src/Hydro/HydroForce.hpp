
#ifndef hydroforcehpp_def__
#define hydroforcehpp_def__

#include <armadillo>

class HydroForce
{
public:
    // Declare class interface methods
    virtual arma::mat CalculateHydrodynamicForces(double time) = 0;
    virtual void LoadHydrodynamicData(void) = 0;
}


#endif // hydroforcehpp_def__