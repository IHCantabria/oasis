#ifndef odedef_hpp__
#define odedef_hpp__
#include <armadillo>
#include <string>
#include "../Simulations/Simulation.hpp"

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
	Simulation *pSim;
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
	// // Number of jacobian evaluations
	// int iJ = 0;
	// // System size
	// int nSystem;
	// // Time, maximum time, output time step
	// double t, tmax, dt_out;
	// // State vector
	// arma::mat y;
	// // Simulation pointer
	// Simulation *pSim;

	// ***Declare constructor***
	BDF2(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation *pIncSim);

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

public:
	// ***Redeclare public attributes if needed***
	// // Solver parameters
	// double dt_max = 1e-2;
	// double dt_min = 1e-9;
	// double dt_ini = 1e-9;
	// double atol = 1e-6;
	// double rtol = 1e-3;
	// int nIterMax = 10;
	// // Number of jacobian evaluations
	// int iJ = 0;
	// // System size
	// int nSystem;
	// // Time, maximum time, output time step
	// double t, tmax, dt_out;
	// // State vector
	// arma::mat y;
	// // Simulation pointer
	// Simulation *pSim;

	// ***Declare constructor***
	BDFN(int N_u, double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation *pIncSim);

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

#endif