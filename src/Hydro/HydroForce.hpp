
#ifndef hydroforcedef_hpp__
#define hydroforcedef_hpp__

#include <armadillo>

class HydroDatabase;

class HydroForce
{
public:
    // Create class methods
    virtual arma::mat CalculateHydrodynamicForces(double time) = 0;
    virtual arma::mat CalculateHydrostaticForces(double time) = 0;
    virtual arma::mat GetCog(void) = 0;
    virtual int GetNumBodies(void) = 0;
    virtual int GetNumPointsIrf(void) = 0;
    virtual arma::mat GetTotalMass(void) = 0;
    virtual void UpdateStructuralMass(arma::mat newStructuralMass) = 0;
    virtual void UpdateTotalMass(void) = 0;
    virtual void LoadHydrodynamicData(std::string file_path) = 0;
    virtual void SetUp(void) = 0;
    virtual void InterpolateHydro(HydroDatabase *pHydro1, HydroDatabase *pHydro2, double interpCoef) = 0;
    virtual void UpdateHydroStiffness(arma::mat newHydrostaticStiffness) = 0;
};

#endif // hydroforcedef_hpp__