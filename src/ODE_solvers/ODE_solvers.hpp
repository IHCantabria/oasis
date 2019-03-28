
#include <armadillo>
#include <string>
#include "../Lines/Lines.hpp"

struct solver_data{
	int nSistema, nSistema2;
	int nLines;
	Line * Lines;
};

class BDF{
private:
	double h_0, h_1, h_2;
	double err, zz, sigma;
	arma::mat y_0, y_1, y_2;
	arma::mat yprime, dy;
	arma::mat I, J, LTE, EWT;
	arma::mat M, F, invM;
	int k, q;
	bool status;

public:
	double dt_max   = 1e-2;
	double dt_min   = 1e-9;
	double dt_ini   = 1e-9;
	double atol     = 1e-6;
	double rtol     = 1e-3;
	int    nIterMax = 10;

	int iJ = 0;

	int nSistema;
	double t, tmax, dt_out;
	arma::mat y;
	solver_data SD;
	arma::mat (*fun) (double, arma::mat, solver_data);

	BDF(void){nSistema=0;};

	BDF(
		double t_u,
		double tmax_u,
		double dt_out_u,
		arma::mat y_u, 
		arma::mat (*fun_u) (double, arma::mat, solver_data),
		solver_data SD_u
		){
			t = t_u;
			tmax = tmax_u;
			dt_out = dt_out_u;
			y = y_u;
			fun = fun_u;
			SD = SD_u;

			h_0 = dt_ini;
			h_1 = dt_ini;
			h_2 = dt_ini;
			dt_max = std::min(0.5*dt_out,dt_max);

			nSistema = y_u.n_rows;

			y_0 = y;
			y_1 = y;
			y_2 = y;
			yprime = arma::zeros(size(y));
			F = arma::zeros(size(y));
			LTE = arma::zeros(size(y));
			EWT = arma::zeros(size(y));

			I = arma::eye(nSistema,nSistema);
			J = arma::zeros(nSistema,nSistema);

			F(0) = 2*atol;
			do{
				jac(t + h_0, y);
				M = I - h_0 * J;
				F = y - y_0 - h_0 * yprime;
				status = arma::solve(dy,M,F,arma::solve_opts::fast);
				if (!status){
					dy = arma::solve(M,F);
				}
				y = y - dy;
			} while((arma::norm(F,"inf") > atol));
			y_0 = y;
			t = t + h_0;

			jac(t + h_0, y_0);
	}
	void step(void);
	void jac(double tt, arma::mat yy);
};