
#include <armadillo>
#include <string>
#include "../Lines/Lines.hpp"

struct solver_data{
	int nLines, nSistema, nSistema2;
	Line * Lines;
};

class BDF{
private:
	arma::mat F, F0;
	arma::mat y0;
	arma::mat M;
	arma::mat dy, ss, yy;
public:
	double dt_max   = 1e-2;
	double dt_min   = 1e-7;
	double dt_ini   = 1e-4;
	double eta_min  = 1e-3;
	double eta_max  = 1e-2;
	double atol     = 1e-6;
	double rho      = 1e+1;
	double sigma    = 1e-2;
	int    nIterMax = 100;

	int nSistema;
	double t, t_prev, tmax, dt, eta;
	arma::mat y, y_prev, y_jac, J;
	solver_data SD;
	arma::mat (*fun) (double, arma::mat, solver_data);
	arma::mat (*jac) (double, arma::mat, solver_data);
	BDF(void){nSistema=0;};
	BDF(
		double t_u,
		double tmax_u,
		arma::mat y_u, 
		arma::mat (*fun_u) (double, arma::mat, solver_data), 
		arma::mat (*jac_u) (double, arma::mat, solver_data),
		solver_data SD_u
		){
			t = t_u;
			tmax = tmax_u;
			y = y_u;
			nSistema = y_u.n_rows;
			fun = fun_u;
			jac = jac_u;
			SD = SD_u;
			dt = dt_ini;
	}
	void step(void);
	void get_next_dt(void);
};