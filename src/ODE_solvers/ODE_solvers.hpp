#ifndef odedef_hpp__
#define odedef_hpp__
#include <armadillo>
#include <string>
#include "../Simulations/Simulation.hpp"

class BDF2
// Class for the implementation of the adaptative BDF2 scheme.
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
	// Nonlinear system iteration counter and jacobean matrix recyling counter
	int k, q;
	// Solver status
	bool status;

public:
	// ***Declare public attributes***
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
	int nSistema;
	// Time, maximum time, output time step
	double t, tmax, dt_out;
	// State vector
	arma::mat y;
	// Simulation pointer
	Simulation *pSim;

	// ***Declare constructor***
	BDF2(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation *pIncSim);

	// ***Declare methods***
	// Evaluate the system dynamics function
	arma::mat fun(double, arma::mat);
	// Initialize the solver
	void Initialize(void);
	// Perform a time step
	void step(void);
	// Compute the jacobian matrix
	void jac(double tt, arma::mat yy);
	// Evaluate the nonlinear BDF2 scheme function
	arma::mat BDF2_fun(double t, arma::mat y);
};

#endif