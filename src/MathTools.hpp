
#ifndef mathtoolsdef_hpp__
#define mathtoolsdef_hpp__

#include <armadillo>

#define M_PI 3.141592653589793 // PI number value (taken from MATLAB)

arma::mat arange(double a, double b, double step);
arma::mat interp1(arma::mat x, arma::mat y, arma::mat xi);
arma::mat linspace(double a, double b, int numPoints);
double trapz(arma::mat y, double h);
double trapzi(arma::mat t, arma::mat y);


#endif // mathtoolsdef_hpp__