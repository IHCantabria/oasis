
#include <armadillo>

class spline{
public:
	int order;
	arma::mat T, COEF;
	spline(void){order=1;};
	spline(arma::mat, arma::mat, int pp);
	arma::mat spl_eval(double);

};