#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include "ODE_solvers.hpp"
#include "../Simulations/Simulation.hpp"

// Methods for the BDF2 class
BDF2::BDF2(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation *pIncSim)
{
	// ***BDF2 constructor***
	// Set the simulation pointer
	pSim = pIncSim;
	// Set the initial time, maximum time, output time step and state vector
	t = t_u;
	tmax = tmax_u;
	dt_out = dt_out_u;
	y = y_u;
	// Initiallize the variables with the last three time steps
	h_0 = dt_ini;
	h_1 = dt_ini;
	h_2 = dt_ini;
	// Compute the maximum time step
	dt_max = std::min(dt_out, dt_max);
	// Compute the system size
	nSystem = y_u.n_rows;
	// Initialize the variables for the last three states
	y_0 = y;
	y_1 = y;
	y_2 = y;
	// Initialize the required vectors and matrices
	yprime = arma::zeros(size(y));
	F = arma::zeros(size(y));
	LTE = arma::zeros(size(y));
	EWT = arma::zeros(size(y));
	I = arma::eye(nSystem, nSystem);
	J = arma::zeros(nSystem, nSystem);
}

arma::mat BDF2::fun(double tt, arma::mat yy)
{
	// ***Interfaze method to the simulation method responsible for the state derivarive evaluation***
	return pSim->CalculateSystemDynamics(tt, yy);
}

void BDF2::jac(double tt, arma::mat yy)
{
	// ***Compute the jacobian matrix of the function that returns the state derivative***
	yprime = fun(tt, yy);
	for (int ii = 0; ii < nSystem; ii = ii + 1)
	{
		J.col(ii) = 1e12 * (fun(tt, yy + 1e-12 * I.col(ii)) - yprime);
	}
}

void BDF2::init()
{
	// ***Initialize the BDF2 solver with one iteration of a BDF1 scheme***
	// Initiallize the function to minimize to a value larger than the tolerance, and the iteration counter
	F(0) = 2 * atol;
	int nIter = 0;
	// Nonlinear system iterative solver, checking for convergence
	do
	{
		// Compute the jacobian matrix (this call also updates yprime)
		jac(t + h_0, y);
		// Increment the jacobian evaluation counter
		iJ = iJ + 1;
		// Compute the nonlinear system matrix
		M = I - h_0 * J;
		// Compute the nonlinear system vector
		F = y - y_0 - h_0 * yprime;
		// Solve for the state change vector
		status = arma::solve(dy, M, F, arma::solve_opts::fast);
		// If the solver fails, solve the system using the slower method
		if (!status)
		{
			dy = arma::solve(M, F);
		}
		// Update the state vector
		y = y - dy;
		// Update the number of iterations
		nIter = nIter + 1;
	} while ((arma::norm(F, "inf") > atol) & (nIter < nIterMax));
	// If the number of iterations is larger than the maximum, return an error
	if (nIter >= nIterMax)
	{
		std::cout << "ERROR: Convergence Failed!" << std::endl;
		throw std::exception();
	}
	// Otherwise, asign the state computed with the BDF1 scheme to the last state and update the solver time
	y_0 = y;
	t = t + h_0;
	// Compute the jacobian matrix at the new state, so it can be reused in the next iterations
	jac(t, y);
}

arma::mat BDF2::BDF2_fun(double t, arma::mat y)
{
	//***Evaluate the nonlinear BDF2 scheme function***
	F = (1.0 + h_0 / (h_1 + h_0)) * y - ((h_1 + h_0) / h_1) * y_0 + ((h_0 * h_0 / h_1) / (h_1 + h_0)) * y_1 - h_0 * fun(t + h_0, y);
	return F;
}

void BDF2::step(void)
{
	// ***Perform a time step with the BDF2 scheme***
	int NN;				  // Declare local maximum number of iterations
	q = 0;				  // Initiallize jacobian recycling counter
	bool flag_nan = true; // Declare NaN flag to true so the solver can try again with a smaller step size
// Loop to try again with a smaller step size if NaN is detected
LOOP:
	// Compute the maximum number of iterations in terms of the number of times the jacobean matrix is recycled
	if (q < 2)
	{
		NN = 100;
	}
	else
	{
		NN = nIterMax;
	}
	// Initiallize the iteration counter
	k = 0;
	// Estimate the newstate with an explicit Euler step and compute the function to minimize
	y = y_0 + h_0 * (y_0 - y_1) / h_1;
	F = BDF2_fun(t + h_0, y);
	// Loop to solve the nonlinear system
	do
	{
		// If the jacobean matrix was recycled twice, compute it again
		if (q >= 2)
		{
			jac(t + h_0, y);
			iJ = iJ + 1;
			q = q + 1;
		}
		// Compute the nonlinear system matrix
		M = (1.0 + h_0 / (h_1 + h_0)) * I - h_0 * J;
		// Solve for the state change vector
		status = arma::solve(dy, M, -F, arma::solve_opts::fast);
		// If the solver fails, solve the system using the slower method
		if (!status)
		{
			dy = arma::solve(M, -F);
		}
		// Compute the line-search step
		// TODO: implement ARMIJO line search
		double rho = 1.0;
		// Update the state vector
		y = y + rho * dy;
		// Compute the function to minimize
		F = BDF2_fun(t + h_0, y);
		// Update the iteration counter
		k = k + 1;
	} while (((arma::norm(dy) > atol + rtol * arma::norm(y)) | (arma::norm(F) > atol)) & (k < NN));
	// Check if the maximum number of iterations was reached
	if (k >= NN)
	{
		if (q < 2)
		{
			// If the jacobean matrix was not updated before, reduce the time step and try again
			// TODO: review this (maybe h_0=dt_min?)
			h_0 = std::max(pow(10.0, -2 * q) * h_0, dt_min);
			jac(t + h_0, y_0 + h_0 * (y_0 - y_1) / h_1);
			iJ = iJ + 1;
			q = q + 1;
			std::cout << "WARNING: In BDF2, Maximum number of iterations reached! Trying again with smaller step size... " << std::endl;
			goto LOOP;
		}
		else
		{
			// Otherwise, return an error
			std::cout << "ERROR: Convergence Failed!" << std::endl;
			throw std::exception();
		}
	}

	// Compute the error weighted tolerance
	EWT = atol * arma::ones(size(y)) + rtol * arma::abs(y);
	// Check for NaN values in the state vector
	if (EWT.has_nan())
	{
		// If NaN is detected for the first time, try again with a smaller step size, otherwise return an error
		if (flag_nan)
		{
			// update the NaN flag
			flag_nan = false;
			// reduce the time step
			h_0 = dt_min;
			// recomput the jacobean matrix
			jac(t + h_0, y_0 + h_0 * (y_0 - y_1) / h_1);
			iJ = iJ + 1;
			std::cout << "	WARNING: In BDF2, NaN detected! Trying again with smaller step size... " << std::endl;
			goto LOOP;
		}
		else
		{
			std::cout << std::endl
					  << "ERROR: NaN detected after trying again!" << std::endl;
			throw std::exception();
		}
	}
	// Compute the local truncation error comparing with a BDF1 scheme
	arma::mat M1 = I - h_0 * J;
	arma::mat F1 = y - y_0 - h_0 * fun(t + h_0, y);
	status = arma::solve(LTE, M1, F1, arma::solve_opts::fast);
	if (!status)
	{
		LTE = arma::solve(M1, F1);
	}
	// Compute the adaptive time step multiplier
	double sigma = pow(0.5 * arma::norm(EWT) / arma::norm(LTE), 0.25);
	// If the error is too large, reduce the time step and try again
	if ((sigma < 0.9) & (h_0 > dt_min))
	{
		h_0 = h_0 * sigma;
		std::cout << "	WARNING: In BDF2, Error too large! Trying again with smaller step size... " << std::endl;
		goto LOOP;
	}
	// Update the time, the last three states and time steps, and the next time step size
	t = t + h_0;
	// Update the last three states
	y_2 = y_1;
	y_1 = y_0;
	y_0 = y;
	// Update the last three time steps
	h_2 = h_1;
	h_1 = h_0;
	// Compute the adaptive time step
	h_0 = sigma * h_0;
	h_0 = std::max(h_0, dt_min);
	h_0 = std::min(h_0, dt_max);
	h_0 = std::min(h_0, dt_out - std::fmod(t, dt_out) + dt_min);
	h_0 = std::min(h_0, tmax - t + h_0);
}

// Methods for the BDFN class
BDFN::BDFN(int N_u, double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation *pIncSim)
{
	// ***BDFN constructor***
	// Check for the order of the BDF scheme to be larger than 2
	if (N_u < 2)
	{
		std::cout << "ERROR: BDFN order must be larger than 1!" << std::endl;
		throw std::exception();
	}
	// Set the order of the BDF scheme
	N = N_u;
	// Set the simulation pointer
	pSim = pIncSim;
	// Set the initial time, maximum time, output time step and state vector
	t = t_u;
	tmax = tmax_u;
	dt_out = dt_out_u;
	y = y_u;
	// Initiallize the time memory vector with the last N times
	tau_i = arma::zeros(N);
	tau_prev = t;
	// Compute the maximum time step
	dt_max = std::min(dt_out, dt_max);
	// Compute the system size
	nSystem = y_u.n_rows;
	// Initialize the states memory vector for the last N states
	y_i = arma::zeros(nSystem, N);
	y_prev = y;
	// Initialize the required vectors and matrices
	yprime = arma::zeros(size(y));
	F = arma::zeros(size(y));
	LTE = arma::zeros(size(y));
	EWT = arma::zeros(size(y));
	I = arma::eye(nSystem, nSystem);
	Jfun = arma::zeros(nSystem, nSystem);
	JF = arma::zeros(nSystem, nSystem);
}

arma::mat BDFN::fun(double tt, arma::mat yy)
{
	// ***Interfaze method to the simulation method responsible for the state derivarive evaluation***
	return pSim->CalculateSystemDynamics(tt, yy);
}

void BDFN::jac(double tt, arma::mat yy)
{
	// ***Compute the jacobian matrix of the function that returns the state derivative***
	yprime = fun(tt, yy);
	for (int ii = 0; ii < nSystem; ii = ii + 1)
	{
		Jfun.col(ii) = 1e12 * (fun(tt, yy + 1e-12 * I.col(ii)) - yprime);
	}
}

void BDFN::init()
{
	// ***Initialize the BDFN solver with N iteration of a BDF1 scheme***
	// Set the initial time step size
	dt = dt_ini;
	// Perform N-1 iterations of a BDF1 scheme
	for (int ii = 0; ii < N; ii = ii + 1)
	{
		// Compute yprime at the current state
		yprime = fun(t, y);
		// Set the previous state
		arma::mat y_old = y;
		// Estimate the new state with an explicit Euler step and compute the function to minimize
		y = y_old + dt * yprime;
		F = y - y_old - dt * yprime;
		// Set the iteration counter to zero
		int nIter = 0;
		// Nonlinear system iterative solver, checking for convergence
		do
		{
			// Compute the jacobian matrix (this call also updates yprime)
			jac(t + dt, y);
			iJ = iJ + 1;
			// Compute the nonlinear system matrix
			JF = I - dt * Jfun;
			// Compute the nonlinear system vector
			F = y - y_old - dt * yprime;
			// Solve for the state change vector
			status = arma::solve(dy, JF, -F, arma::solve_opts::fast);
			// If the solver fails, solve the system using the slower method
			if (!status)
			{
				dy = arma::solve(JF, -F);
			}
			// Update the state vector
			y = y + dy;
			// Compute the function to minimize
			F = y - y_old - dt * fun(t + dt, y);
			// Update the number of iterations
			nIter = nIter + 1;
		} while ((arma::norm(F, "inf") > atol) & (nIter < nIterMax));
		// If the number of iterations is larger than the maximum, return an error
		if (nIter >= nIterMax)
		{
			std::cout << "ERROR: Convergence failed in BDFN during init!" << std::endl;
			throw std::exception();
		}
		// Otherwise, update the time and save the state computed with the BDF1 scheme to the corresponding state vector
		t = t + dt;
		y_i.col(ii) = y;
		tau_i(ii) = t;
	}
	// Compute the jacobian matrix at the next state estimate, so it can be reused in the next iterations
	jac(t + dt, y_i.col(N - 1) + dt * fun(t, y_i.col(N - 1)));
}

arma::mat BDFN::BDFN_fun(double t, arma::mat y)
{
	//***Evaluate the nonlinear BDFN scheme function***
	double tmp_coef;
	arma::mat FF = fun(t, y);
	tmp_coef = 0.0;
	for (int jj = 0; jj < N; jj = jj + 1)
	{
		tmp_coef = tmp_coef + 1.0 / (t - tau_i(jj));
	}
	FF = FF - y * tmp_coef;
	for (int ii = 0; ii < N; ii = ii + 1)
	{
		tmp_coef = 1.0 / (tau_i(ii) - t);
		for (int jj = 0; jj < N; jj = jj + 1)
		{
			if (jj != ii)
			{
				tmp_coef = tmp_coef * (t - tau_i(jj)) / (tau_i(ii) - tau_i(jj));
			}
		}
		FF = FF - y_i.col(ii) * tmp_coef;
	}
	return FF;
}

void BDFN::step(void)
{
	// ***Perform a time step with the BDF2 scheme***
	int NN;				  // Declare local maximum number of iterations
	double tmp_coef;	  // Declare local variable for the temporal coefficient
	q = 0;				  // Initialize jacobian recycling counter
	bool flag_nan = true; // Declare NaN flag to true so the solver can try again with a smaller step size
// Loop to try again with a smaller step size if NaN is detected
LOOP:
	// Compute the maximum number of iterations in terms of the number of times the jacobean matrix is recycled
	if (q < 2)
	{
		NN = 100;
	}
	else
	{
		NN = nIterMax;
	}
	// Initialize the iteration counter
	k = 0;
	// Estimate the newstate with an explicit Euler step and compute the function to minimize
	yprime = fun(t, y_i.col(N - 1)); // TODO: this call can be avoided
	y = y_i.col(N - 1) + dt * yprime;
	F = BDFN_fun(t + dt, y);
	dy = y;
	// Loop to solve the nonlinear system
	while ((arma::norm(F, "inf") > atol) & (k <= NN))
	{
		// If the jacobean matrix was recycled twice, compute it again
		if (q >= 2)
		{
			jac(t + dt, y);
			iJ = iJ + 1;
		}
		// Compute the nonlinear system matrix
		tmp_coef = 0.0;
		for (int jj = 0; jj < N; jj = jj + 1)
		{
			tmp_coef = tmp_coef + 1.0 / ((t + dt) - tau_i(jj));
		}
		JF = Jfun - I * tmp_coef;
		// Solve for the state change vector
		status = arma::solve(dy, JF, -F, arma::solve_opts::fast);
		// If the fast solver fails, solve the system using the slower method
		if (!status)
		{
			dy = arma::solve(JF, -F);
		}
		// Compute the line-search step
		double alpha = 1e-4;
		double lambda = 0.9;
		double rho_min = 1e-2;
		double rho = 1.0;
		arma::mat F_new = BDFN_fun(t + dt, y + rho * dy);
		while ((arma::norm(F_new) > (1.0 - alpha) * arma::norm(F)) & (rho > rho_min))
		{
			rho = rho * lambda;
			F_new = BDFN_fun(t + dt, y + rho * dy);
		}
		// Update the state vector
		y = y + rho * dy;
		// Recompute the function to minimize
		F = F_new;
		// Update the iteration counter
		k = k + 1;
	}
	// Check if the maximum number of iterations was reached
	if (k >= NN)
	{
		if (q < 2)
		{
			// If the jacobean matrix was not updated before, reduce the time step and try again
			if (dt > 100 * dt_min)
			{
				dt = std::max(dt / 100, dt_min);
				q = q + 1;
				std::cout << "WARNING: In BDFN, Maximum number of iterations reached! Trying again with smaller step size... " << std::endl;
				goto LOOP;
			}
			else
			{
				q = 2;
				dt = dt_min;
				std::cout << "WARNING: In BDFN, Maximum number of iterations reached! " << std::endl;
				std::cout << "         Trying again with minimum step size and Jacobean matrix recomputation... " << std::endl;
				goto LOOP;
			}
		}
		else
		{
			// Otherwise, return an error
			std::cout << "ERROR: Convergence Failed!" << std::endl;
			throw std::exception();
		}
	}

	// Compute the error weighted tolerance
	EWT = atol * arma::ones(size(y)) + rtol * arma::abs(y);
	// Check for NaN values in the state vector
	if (EWT.has_nan())
	{
		// If NaN is detected for the first time, try again with a smaller step size, otherwise return an error
		if (flag_nan)
		{
			// update the NaN flag
			flag_nan = false;
			// reduce the time step
			dt = dt_min;
			// recompute the jacobean matrix
			y = y_i.col(N - 1) + dt * (y_i.col(N - 1) - y_i.col(N - 2)) / (tau_i(N - 1) - tau_i(N - 2));
			jac(t + dt, y);
			iJ = iJ + 1;
			std::cout << "	WARNING: In BDF2, NaN detected! Trying again with smaller step size... " << std::endl;
			goto LOOP;
		}
		else
		{
			std::cout << std::endl
					  << "ERROR: NaN detected after trying again!" << std::endl;
			throw std::exception();
		}
	}
	// Compute the local truncation error comparing with a BDF(N-1) scheme
	LTE = arma::zeros(size(y));
	arma::mat y_LTE = arma::join_horiz(arma::join_horiz(y_prev, y_i), y);
	arma::vec tau_LTE = arma::join_vert(arma::join_vert(tau_prev * arma::ones(1), tau_i), (t + dt) * arma::ones(1));
	for (int ii = 0; ii < N + 2; ii = ii + 1)
	{
		tmp_coef = 1.0;
		for (int jj = 0; jj < N + 2; jj = jj + 1)
		{
			if (jj != ii)
			{
				tmp_coef = tmp_coef / (tau_LTE(ii) - tau_LTE(jj));
			}
		}
		LTE = LTE + y_LTE.col(ii) * tmp_coef;
	}
	LTE = LTE * pow(dt, N + 1);
	// Compute the adaptive time step multiplier
	double sigma = pow(0.5 * arma::norm(EWT) / arma::norm(LTE), 0.25);
	// If the error is too large, reduce the time step and try again
	if ((sigma < 0.9) & (dt > dt_min))
	{
		dt = std::max(dt * sigma, dt_min);
		std::cout << "	WARNING: In BDFN, Error too large! Trying again with smaller step size... " << std::endl;
		goto LOOP;
	}
	// Set an upper bound for sigma
	sigma = std::min(sigma, 2.0);
	// Update the time, the last three states and time steps, and the next time step size
	t = t + dt;
	// Update the last states and time steps in memory
	y_prev = y_i.col(0);
	tau_prev = tau_i(0);
	for (int ii = 0; ii < N - 1; ii = ii + 1)
	{
		y_i.col(ii) = y_i.col(ii + 1);
		tau_i(ii) = tau_i(ii + 1);
	}
	y_i.col(N - 1) = y;
	tau_i(N - 1) = t;
	// Compute the adaptive time step
	dt = sigma * dt;
	dt = std::max(dt, dt_min);
	dt = std::min(dt, dt_max);
	dt = std::min(dt, dt_out - std::fmod(t, dt_out) + dt_min);
	dt = std::min(dt, tmax - t + dt_min);
}
