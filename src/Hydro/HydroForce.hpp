
#ifndef hydroforcedef_hpp__
#define hydroforcedef_hpp__

#include <armadillo>

class HydroForce
{
public:

    // Create constructors and destructor
    HydroForce(){};
    virtual ~HydroForce() = 0;

    // Create class methods
    virtual arma::mat CalculateHydrodynamicForces(double time) = 0;
    virtual arma::mat GetCog(void)=0;
    // virtual arma::mat GetInertiaMatrixInv(void)=0;
    virtual int GetNumBodies(void)=0;
	virtual int GetNumPointsIrf(void)=0;
    virtual void LoadHydrodynamicData(std::string file_path)=0;
};

#endif // hydroforcedef_hpp__