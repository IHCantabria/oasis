
#include <armadillo>

class spline{
public:
	int order;
	arma::mat T, COEF;
	spline(void){order=1;};
	spline(arma::mat, arma::mat, int pp);
	arma::mat spl_eval(double);
};

class solver_data_base{
public:
	int dummy;
};


class BDF{
public:
	int nSistema;
	double t, dt_max, dt_min, dt_ini, dt, atol, rtol;
	arma::mat y, yprime;
	solver_data_base* SD;
	arma::mat (*fun) (double, arma::mat, solver_data_base*);
	arma::mat (*jac) (double, arma::mat, solver_data_base*);
	BDF(void){nSistema=0;};
	BDF(
		double t_u, 
		arma::mat y_u, 
		arma::mat (*fun_u) (double, arma::mat, solver_data_base*), 
		arma::mat (*jac_u) (double, arma::mat, solver_data_base*),
		solver_data_base* SD_u,
		double dt_max_u = 1e0, 
		double dt_min_u = 1e-8, 
		double dt_ini_u = 1e-4, 
		double atol_u = 1e-6, 
		double rtol_u = 1e-12
	){
		t = t_u;
		y = y_u;
		yprime = y_u*0.0;
		nSistema = y_u.n_rows;
		fun = fun_u;
		jac = jac_u;
		SD = SD_u;
		dt_max = dt_max_u; 
		dt_min = dt_min_u; 
		dt_ini = dt_ini_u; 
		atol = atol_u; 
		rtol = rtol_u;
		dt = dt_ini;
	}
	void step(void){
		y = y + dt * fun(t,y,SD);
		t = t + dt;
	}
};
