
#ifndef winddef_hpp__
#define winddef_hpp__
#include <armadillo>
#include <string>
#include <cstdio>

class Wind
{
public:
    // Declare class variables
    double pi = arma::datum::pi;
    double speed;
    double heading;

    // Declare class constructors
    Wind(double S, double D);

};


#endif // wavedef_hpp__