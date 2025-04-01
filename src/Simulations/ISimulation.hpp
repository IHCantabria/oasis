
#ifndef isimulationdef_hpp__
#define isimulationdef_hpp__

#include <string>
#include <armadillo>

class ISimulation
{
public:
    // Declare general purpose class methods
    virtual arma::mat CalculateSystemDynamics(double time, arma::mat y) = 0;
};

#endif // isimulation_hpp__
