
#ifndef hydroforcehpp_def__
#define hydroforcehpp_def__

#include <armadillo>
#include <string>

class HydroForce
{
public:
    // Declare class interface methods
    virtual arma::mat CalculateHydrodynamicForces(double time) = 0;
    virtual arma::mat GetCog(void) = 0;
    virtual int GetNumPointsIrf(void) = 0;
    virtual void LoadHydrodynamicData(std::string filePath) = 0;
};


#endif // hydroforcehpp_def__