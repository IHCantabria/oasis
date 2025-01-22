
#ifndef mathtoolsdef_hpp__
#define mathtoolsdef_hpp__

#include <armadillo>

arma::mat arange(double a, double b, double step);
arma::mat interp1(arma::vec x, arma::mat y, arma::vec xi);
arma::cube interp1(arma::vec x, arma::cube y, arma::vec xi);
arma::cube interp2(arma::vec x, arma::vec y, arma::cube z, arma::vec xi, arma::vec yi);
arma::mat linspace(double a, double b, int numPoints);
arma::mat mod(arma::mat a, double x);
double trapz(arma::mat y, double h);
double trapzi(arma::mat t, arma::mat y);
std::tuple<arma::vec, arma::vec> upcrossing(arma::vec t, arma::vec u);
arma::cube permute(arma::cube x, int ind);
arma::mat wrapTo180(arma::mat x);
arma::mat wrapToPi(arma::mat x);
arma::uvec comp_ind(int n, arma::uvec ind);
arma::umat comb_n_k(int n, int k);
double step(double x, double x0, double h0, double x1, double h1);
arma::mat triangleChangeFrame(arma::mat V0, arma::mat V1, arma::mat V2);

std::tuple<arma::mat, arma::uvec> unique_rows(arma::mat &x);
arma::uvec unique(arma::uvec &v);
arma::mat sort_rows(arma::mat x, int icol);
arma::umat indMat(arma::uvec ind, arma::umat x);

#endif // mathtoolsdef_hpp__