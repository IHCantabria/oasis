#ifndef odedef_hpp__
#define odedef_hpp__
#include <armadillo>
#include <string>
#include <tuple>
#include "../Simulations/ISimulation.hpp"

class ODE_solver
// Base class for the implementation of ODE solvers.
{
private:
	bool status = true;

public:
	// Solver parameters
	double dt_max = 1e-2;
	double dt_min = 1e-9;
	double dt_ini = 1e-9;
	double atol = 1e-6;
	double rtol = 1e-3;
	int nIterMax = 10;
	// Number of jacobian evaluations
	int iJ = 0;
	// System size
	int nSystem;
	// Time, maximum time, output time step
	double t, tmax, dt_out;
	// State vector
	arma::mat y;
	// Simulation pointer
	ISimulation *pSim;
	// Evaluate the system dynamics function
	virtual arma::mat fun(double t, arma::mat y) = 0;
	// Initialize the solver
	virtual void init(void) = 0;
	// Perform a time step
	virtual void step(void) = 0;
	// Compute the jacobian matrix
	virtual void jac(double t, arma::mat y) = 0;
};

class BDF2 : public ODE_solver
// Class for the implementation of the adaptive BDF2 scheme.
// Based on “Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s”
// by Celaya et al. (2014)
{
private:
	// ***Declare private attributes***
	// Last three time steps
	double h_0, h_1, h_2;
	// Last three states
	arma::mat y_0, y_1, y_2;
	// State derivative vector and state change vector
	arma::mat yprime, dy;
	// Identity matrix, Jacobian matrix, Local Truncation Error, Error Weighted Tolerance
	arma::mat I, J, LTE, EWT;
	// Nonlinear system matrix, vector and matrix inverse
	arma::mat M, F, invM;
	// Nonlinear system iteration counter and jacobean matrix recycling counter
	int k, q;
	// Solver status
	bool status;

public:
	// ***Redeclare public attributes if needed***
	// // Solver parameters
	// double dt_max = 1e-2;
	// double dt_min = 1e-9;
	// double dt_ini = 1e-9;
	// double atol = 1e-6;
	// double rtol = 1e-3;
	// int nIterMax = 10;

	// ***Declare constructor***
	BDF2(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, ISimulation *pIncSim);

	// ***Declare methods***
	// Evaluate the system dynamics function
	arma::mat fun(double, arma::mat);
	// Initialize the solver
	void init(void);
	// Perform a time step
	void step(void);
	// Compute the jacobian matrix
	void jac(double tt, arma::mat yy);
	// Evaluate the nonlinear BDF2 scheme function
	arma::mat BDF2_fun(double t, arma::mat y);
};

class BDFN : public ODE_solver
// Class for the implementation of the adaptive BDF scheme of generic order N.
// Based on BDF2, generalized.
{
private:
	// ***Declare private attributes***
	// Current time step
	double dt;
	// Last N+1 time steps
	arma::vec tau_i;
	// Previous time step for LTE computation
	double tau_prev;
	// Last N+1 states
	arma::mat y_i;
	// Previous state for LTE computation
	arma::mat y_prev;
	// State derivative vector and state change vector
	arma::mat yprime, dy;
	// Identity matrix, Jacobian matrix, Local Truncation Error, Error Weighted Tolerance
	arma::mat I, Jfun, LTE, EWT;
	// Nonlinear system matrix, vector and matrix inverse
	arma::mat JF, F, invJF;
	// Nonlinear system iteration counter and jacobean matrix recycling counter
	int k, q;
	// Solver status
	bool status;
	// Order of the BDF scheme
	int N;
	// Time adaptivity flag
	bool adaptivity;
	// BDFN constant time step coefficients
	arma::vec y_coefs;
	double f_coef;

	// ***Declare private methods if needed***
	std::tuple<arma::vec, double> set_coefs(int N_u);

public:
	// ***Redeclare public attributes if needed***
	// // Solver parameters
	// double dt_max = 1e-2;
	// double dt_min = 1e-9;
	// double dt_ini = 1e-9;
	// double atol = 1e-6;
	// double rtol = 1e-3;
	// int nIterMax = 10;

	// ***Declare constructor***
	BDFN(int N_u, bool a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u, ISimulation *pIncSim);

	// ***Declare methods***
	// Evaluate the system dynamics function
	arma::mat fun(double, arma::mat);
	// Initialize the solver
	void init(void);
	// Perform a time step
	void step(void);
	// Compute the jacobian matrix
	void jac(double tt, arma::mat yy);
	// Evaluate the nonlinear BDF2 scheme function
	arma::mat BDFN_fun(double t, arma::mat y);
};

class ESDIRK : public ODE_solver
// Class for the implementation of the adaptive ESDIRK scheme.
// Based on "Fourth-order Runge–Kutta schemes for fluid mechanics applications. (Carpenter, 2005)"
{
private:
	// ***Declare private attributes***
	// Current time step
	double dt;
	// Identity matrix, Jacobian matrix
	arma::mat I, Jfun;
	// System vector at each stage
	arma::mat ys;
	// System vector derivative at each stage
	arma::mat ys_prime;
	// System vector derivative
	arma::mat y_prime;
	// Embedded solution vector
	arma::mat yhat;
	// Time at each stage
	arma::vec ts;
	// Butcher coefficients for the ESDIRK scheme
	int s = 6;
	arma::vec beta = {0.0, 1.0 / 2.0, 83.0 / 250.0, 31.0 / 50.0, 17.0 / 20.0, 1.0};
	arma::mat a = {
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1.0 / 4.0, 1.0 / 4.0, 0.0, 0.0, 0.0, 0.0},
		{8611.0 / 62500.0, -1743.0 / 31250.0, 1.0 / 4.0, 0.0, 0.0, 0.0},
		{5012029.0 / 34652500.0, -654441.0 / 2911500.0, 174375.0 / 388108.0, 1.0 / 4.0, 0.0, 0.0},
		{15267082809.0 / 155376265600.0, -71443401.0 / 120774400.0, 730878875.0 / 902184768.0, 2285395.0 / 8070912.0, 1.0 / 4.0, 0.0},
		{82889.0 / 524892.0, 0.0, 15625.0 / 83664.0, 69875.0 / 102672.0, -2260.0 / 8211.0, 1.0 / 4.0}};
	arma::vec b = {4586570599.0 / 29645900160.0, 0.0, 178811875.0 / 945068544.0, 814220225.0 / 1159782912.0, -3700637.0 / 11593932.0, 61727.0 / 225920.0};
	// Solver status
	bool status;
	// Time steps counter
	int nSteps = 0;
	// Error ratios
	double error_ratio = 1.0;
	double error_ratio_old;
	// Previous time step
	double dt_old;
	// Time adaptivity flag
	bool adaptivity;
	// Local Truncation Error
	double LTE = 1.0;
	double LTE_old = 1.0;

public:
	// ***Redeclare public attributes if needed***
	// // Solver parameters
	// double dt_max = 1e-2;
	// double dt_min = 1e-9;
	// double dt_ini = 1e-9;
	// double atol = 1e-6;
	// double rtol = 1e-3;
	// int nIterMax = 10;

	// ***Declare constructor***
	ESDIRK(bool a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u, ISimulation *pIncSim);

	// ***Declare methods***
	// Evaluate the system dynamics function
	arma::mat fun(double, arma::mat);
	// Initialize the solver
	void init(void);
	// Perform a time step
	void step(void);
	// Compute the jacobian matrix
	void jac(double tt, arma::mat yy);
	// Evaluate the nonlinear ESDIRK scheme function
	arma::mat ESDIRK_fun(double t_i, arma::mat y_i, int ii);
};

#endif