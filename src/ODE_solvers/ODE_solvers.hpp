// SPDX-License-Identifier: GPL-3.0-or-later
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
    // Number of times using the Newton method
    double nNewton = 0;
    // Number of newton iterations
    double nNewtonIter = 0;
    // Average Newton iterations per try
    double nNewtonIterAvg = 0.0;
    // Number of times convergence failed
    int nConvergenceFailed = 0;
    // Time, maximum time, output time step
    double t, tmax, dt_out;
    // State vector
    arma::mat y;
    // Simulation pointer
    ISimulation* pSim;

    // File pointer for the debugging output
    FILE* pfile;
    // Debug flag
    bool debug_flag;

    // Evaluate the system dynamics function
    virtual arma::mat fun(double t, arma::mat y) = 0;
    // Initialize the solver
    virtual void init(void) = 0;
    // Finalize the solver
    virtual void finalize(void) = 0;
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
    BDF2(double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u, ISimulation* pIncSim);

    // ***Declare methods***
    // Evaluate the system dynamics function
    arma::mat fun(double, arma::mat);
    // Initialize the solver
    void init(void);
    // Finalize the solver
    void finalize(void);
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
    BDFN(int N_u, int a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u,
         ISimulation* pIncSim);

    // ***Declare methods***
    // Evaluate the system dynamics function
    arma::mat fun(double, arma::mat);
    // Initialize the solver
    void init(void);
    // Finalize the solver
    void finalize(void);
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
    int s;
    arma::vec beta;
    arma::mat a;
    arma::vec b;
    // Solver status

    bool status;
    // Time steps counter
    int nSteps = 0;
    int nStepsMax;
    // Adaptivity type
    int adaptivity_type;            // 1: Basic, 2: Noventa 2018, 3: Ranocha 2024
    int adaptivity_smooth_flag;     // Multiplier smoothing [1: Min-Max, 2: atan]
    double adaptivity_smooth_param; // Smoothing parameter
    int LTE_norm_type;              // 1: L1, 2: L2, 3: Linf, 4: RMS
    double rho_min;                 // Minimum multiplier to accept a time step
    // Error ratios
    double error_ratio = 1.0;
    double error_ratio_old = 1.0;
    double error_ratio_old2;
    // Previous time steps
    double dt_old;
    double dt_old2;
    // Flag for truncated dt for writing output
    bool dt_truncated_flag = false;
    // Time step before truncation
    double dt_truncated;
    // LTE before truncation
    double LTE_truncated;
    // Time adaptivity flag
    bool adaptivity;
    // Error Weighted Tolerance
    double EWT;
    // Local Truncation Error
    double LTE;
    // Time step multiplier
    double rho;

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
    ESDIRK(bool a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, int nStepsMax_u, arma::mat y_u,
           ISimulation* pIncSim);

    // ***Declare methods***
    // Evaluate the system dynamics function
    arma::mat fun(double, arma::mat);
    // Initialize the solver
    void init(void);
    // Finalize the solver
    void finalize(void);
    // Perform a time step
    void step(void);
    // Compute the jacobian matrix
    void jac(double tt, arma::mat yy);
    // Evaluate the nonlinear ESDIRK scheme function
    arma::mat ESDIRK_fun(double t_i, arma::mat y_i, int ii);
    // Computes and sets the next time step size
    void set_dt(void);
    // Compute LTE, EWT, error ratio and time step multiplier
    void compute_dt(void);
};

#endif