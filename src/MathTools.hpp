
#ifndef mathtoolsdef_hpp__
#define mathtoolsdef_hpp__

#include <armadillo>

#define M_PI 3.141592653589793 // PI number value (taken from MATLAB)

arma::mat arange(double a, double b, double step);
arma::mat interp1(arma::mat x, arma::mat y, arma::mat xi);
arma::cube interp1(arma::mat x, arma::cube y, arma::mat xi);
arma::cube interp2(arma::mat x, arma::mat y, arma::cube z, arma::mat xi, arma::mat yi);
arma::mat linspace(double a, double b, int numPoints);
arma::mat mod(arma::mat a, double x);
double trapz(arma::mat y, double h);
double trapzi(arma::mat t, arma::mat y);
std::tuple<arma::mat,arma::mat> upcrossing(arma::mat t, arma::mat u);
arma::cube permute(arma::cube x, int ind);



#endif // mathtoolsdef_hpp__