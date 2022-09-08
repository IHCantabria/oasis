
#ifndef mathtoolsdef_hpp__
#define mathtoolsdef_hpp__

#include <armadillo>

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
arma::mat wrapTo180(arma::mat x);
arma::mat wrapToPi(arma::mat x);
arma::uvec comp_ind(int n,arma::uvec ind);
arma::umat comb_n_k(int n,int k);
double step(double x, double x0, double h0, double x1, double h1);

std::tuple<arma::mat,arma::uvec> unique_rows(arma::mat& x);
std::tuple<arma::uvec,arma::uvec> unique(arma::uvec& v);
arma::mat sort_rows(arma::mat x, int icol);
arma::umat indMat(arma::uvec ind, arma::umat x);


#endif // mathtoolsdef_hpp__