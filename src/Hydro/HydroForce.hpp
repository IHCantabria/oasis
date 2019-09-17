
#ifndef hydroforcedef_hpp__
#define hydroforcedef_hpp__

class HydroForce
{
public:

    // Create constructors
    HydroForce(){};

    // Create class methods
    virtual double CalculateHydrodynamicForces() = 0;
    virtual void Load() = 0;
};

#endif // hydroforcedef_hpp__