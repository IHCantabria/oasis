#include <iostream>
#include <tuple>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include "ODE_solvers.hpp"
#include "../Simulations/ISimulation.hpp"

// Methods for the BDF2 class
BDF2::BDF2(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, ISimulation *pIncSim)
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
	// Initialize the iteration counter
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
BDFN::BDFN(int N_u, bool a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u, ISimulation *pIncSim)
{
	// ***BDFN constructor***
	// Check for the order of the BDF scheme to be larger than 2
	if (N_u < 1)
	{
		std::cout << "ERROR: BDFN order must be larger or equal than 1!" << std::endl;
		throw std::exception();
	}
	// Check for the order of the BDF scheme to be smaller than 7
	if (N_u > 6)
	{
		std::cout << "ERROR: BDFN order must be smaller or equal than 6!" << std::endl;
		throw std::exception();
	}
	// Check for the output time step to be larger than the maximum time step
	if (dt_out_u < dt_max_u)
	{
		std::cout << "ERROR: Output time step must be larger or equal than the maximum time step!" << std::endl;
		throw std::exception();
	}
	// If adaptivity is not used, check if the output time step is a multiple of the maximum time step
	if (!a_u & (std::fmod(dt_out_u, dt_max_u) > 1e-15) & (std::fmod(dt_out_u, dt_max_u) < dt_max_u - 1e-15))
	{
		std::cout << "ERROR: Output time step must be a multiple of the maximum time step!" << std::endl;
		throw std::exception();
	}
	// If adaptivity is not used, check if the total simulation time is a multiple of the maximum time step
	if (!a_u & (std::fmod(tmax_u, dt_max_u) > 1e-14) & (std::fmod(tmax_u, dt_max_u) < dt_max_u - 1e-14))
	{
		std::cout << "ERROR: Total simulation time must be a multiple of the maximum time step!" << std::endl;
		std::cout << "tmax_u: " << tmax_u << std::endl;
		std::cout << "dt_max_u: " << dt_max_u << std::endl;
		std::cout << "fmod: " << std::fmod(tmax_u, dt_max_u) << std::endl;
		throw std::exception();
	}
	// Set the order of the BDF scheme
	N = N_u;
	// Set the adaptive time step flag
	adaptivity = a_u;
	// Set the simulation pointer
	pSim = pIncSim;
	// Set the initial time, maximum time, output time step and state vector
	t = t_u;
	tmax = tmax_u;
	// dt_out = dt_out_u; // TODO: check why this causes convergence problems
	dt_out = dt_max_u;
	y = y_u;
	// Initiallize the time memory vector with the last N times
	tau_i = arma::zeros(N);
	tau_prev = t;
	// Compute the maximum time step
	dt_max = std::min(dt_out, dt_max_u);
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

std::tuple<arma::vec, double> BDFN::set_coefs(int N_u)
{
	// ***Set the coefficients for the BDFN scheme***
	// Set the coefficients for the BDFN scheme depending on the order (1-6)
	arma::vec tmp_y_coefs;
	double tmp_f_coef;
	if (N_u == 1)
	{
		tmp_y_coefs = {-1.0};
		tmp_f_coef = dt;
	}
	else if (N_u == 2)
	{
		tmp_y_coefs = {1.0 / 3.0, -4.0 / 3.0};
		tmp_f_coef = 2.0 * dt / 3.0;
	}
	else if (N_u == 3)
	{
		tmp_y_coefs = {-2.0 / 11.0, 9.0 / 11.0, -18.0 / 11.0};
		tmp_f_coef = 6.0 * dt / 11.0;
	}
	else if (N_u == 4)
	{
		tmp_y_coefs = {3.0 / 25.0, -16.0 / 25.0, 36.0 / 25.0, -48.0 / 25.0};
		tmp_f_coef = 12.0 * dt / 25.0;
	}
	else if (N_u == 5)
	{
		tmp_y_coefs = {-12.0 / 137.0, 75.0 / 137.0, -200.0 / 137.0, 300.0 / 137.0, -300.0 / 137.0};
		tmp_f_coef = 60.0 * dt / 137.0;
	}
	else if (N_u == 6)
	{
		tmp_y_coefs = {10.0 / 147.0, -72.0 / 147.0, 225.0 / 147.0, -400.0 / 147.0, 450.0 / 147.0, -360.0 / 147.0};
		tmp_f_coef = 60.0 * dt / 147.0;
	}
	else
	{
		std::cout << "ERROR: BDFN order must be between 1 and 6!" << std::endl;
		throw std::exception();
	}
	return std::make_tuple(tmp_y_coefs, tmp_f_coef);
}

arma::mat BDFN::fun(double tt, arma::mat yy)
{
	// ***Interface method to the simulation method responsible for the state derivative evaluation***
	return pSim->CalculateSystemDynamics(tt, yy);
}

void BDFN::jac(double tt, arma::mat yy)
{
	// ***Compute the jacobian matrix of the function that returns the state derivative***
	yprime = fun(tt, yy);
	for (int ii = 0; ii < nSystem; ii = ii + 1)
	{
		Jfun.col(ii) = 1.0e14 * (fun(tt, yy + 1.0e-14 * I.col(ii)) - yprime);
	}
}

void BDFN::init()
{
	// ***Initialize the BDFN solver with N iteration of a BDF1 scheme***
	// Set the initial time step size
	if (adaptivity)
	{
		dt = dt_ini;
	}
	else
	{
		dt = dt_max;
	}
	// Create the initialization buffer
	arma::mat y_i_tmp = arma::zeros(nSystem, N + 1);
	y_i_tmp.col(0) = y;

	// // Hardcode the initial values of a exponential decay with C=1.0
	// int i0 = std::min(0, N);
	// for (int ii = 0; ii < N - i0; ii = ii + 1)
	// {
	// 	t = (ii + 1) * dt;
	// 	y = std::exp(-t);
	// 	y_i_tmp.col(ii + 1) = y;
	// 	y_i.col(ii) = y;
	// 	tau_i(ii) = t;
	// }

	// Compute the first N states with increasing order BDF schemes
	for (int ii = 0; ii < N; ii = ii + 1)
	// for (int ii = N - i0; ii < N; ii = ii + 1)
	{
		// Set the order of the BDF scheme for the current iteration
		int order = ii + 1;
		// Get the BDF coefficients for the current iteration
		arma::vec tmp_y_coefs;
		double tmp_f_coef;
		std::tie(tmp_y_coefs, tmp_f_coef) = set_coefs(order);
		// Estimate the new state with an explicit Euler step and compute the function to minimize
		y = y + dt * fun(t, y);
		F = 2 * atol * arma::ones(size(y));
		// Set dy = y
		dy = y;
		// Set the iteration counter to zero
		int k = 0;
		// Nonlinear system iterative solver, checking for convergence
		while (((arma::norm(dy, "inf") > atol) | (arma::norm(F, "inf") > atol)) & (k <= nIterMax))
		{
			// Compute the jacobian matrix
			jac(t + dt, y);
			iJ = iJ + 1;
			// Compute the nonlinear system matrix
			JF = I - tmp_f_coef * Jfun;
			// Solve for the state change vector
			status = arma::solve(dy, JF, -F, arma::solve_opts::fast);
			// If the solver fails, solve the system using the slower method
			if (!status)
			{
				dy = arma::solve(JF, -F);
			}
			// Compute the line-search step
			double alpha = 1e-4;
			double lambda = 0.9;
			double rho_min = 1e-2;
			double rho = 1.0;
			arma::mat F_new = y + rho * dy - tmp_f_coef * fun(t + dt, y + rho * dy);
			for (int jj = 0; jj < order; jj = jj + 1)
			{
				F_new = F_new + arma::as_scalar(tmp_y_coefs(jj)) * y_i_tmp.col(jj);
			}
			while ((arma::norm(F_new) > (1.0 - alpha) * arma::norm(F)) & (rho > rho_min))
			{
				rho = rho * lambda;
				F_new = y + rho * dy - tmp_f_coef * fun(t + dt, y + rho * dy);
				for (int jj = 0; jj < order; jj = jj + 1)
				{
					F_new = F_new + arma::as_scalar(tmp_y_coefs(jj)) * y_i_tmp.col(jj);
				}
			}
			// Update the state vector
			y = y + rho * dy;
			// Recompute the function to minimize
			F = F_new;
			// Update the number of iterations
			k = k + 1;
		}
		// If the number of iterations is larger than the maximum, return an error
		if (k >= nIterMax)
		{
			std::cout << "ERROR: Convergence failed in BDFN during init!" << std::endl;
			throw std::exception();
		}
		// Otherwise, update the time and save the state computed with the BDF scheme to the corresponding state vector
		t = t + dt;
		y_i.col(ii) = y;
		tau_i(ii) = t;
		// Update the initialization buffer
		y_i_tmp.col(ii + 1) = y;
	}
	// Set the constant coefficients for the BDFN scheme if not adaptive
	if (!adaptivity)
	{
		std::tie(y_coefs, f_coef) = set_coefs(N);
	}
	// Compute the jacobian matrix at the current state estimate, so it can be reused in the next iterations
	jac(t, y);
}

arma::mat BDFN::BDFN_fun(double t, arma::mat y)
{
	//***Evaluate the nonlinear BDFN scheme function***
	arma::mat FF;

	if (adaptivity)
	{
		double tmp_coef_y = 0.0;
		for (int jj = 0; jj < N; jj = jj + 1)
		{
			tmp_coef_y = tmp_coef_y + 1.0 / (t - tau_i(jj));
		}
		FF = y - fun(t, y) / tmp_coef_y;
		double tmp_coef;
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
			FF = FF + y_i.col(ii) * tmp_coef / tmp_coef_y;
		}
	}
	else
	{
		FF = y - f_coef * fun(t, y);
		for (int ii = 0; ii < N; ii = ii + 1)
		{
			FF = FF + arma::as_scalar(y_coefs(ii)) * y_i.col(ii);
		}
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
	// If adaptivity is not used, make sure that the current time is a multiple of the time step
	if (!adaptivity)
	{
		t = dt * round(t / dt);
	}
	// Estimate the newstate with an explicit Euler step and compute the function to minimize
	arma::mat y_ini = y_i.col(N - 1) + dt * fun(t, y_i.col(N - 1));
	arma::mat F_ini = BDFN_fun(t + dt, y_ini);
// Loop to try again with a smaller step size if NaN is detected
LOOP:
	if ((!adaptivity) && (dt < dt_max))
	{
		std::cout << std::endl
				  << "ERROR: In BDFN, without adaptivity, it was not possible to try again with a smaller step size!" << std::endl;
		throw std::exception();
	}
	// Compute the maximum number of iterations in terms of the number of times the jacobean matrix is recycled
	if (q <= 2)
	{
		NN = 100;
	}
	else
	{
		NN = nIterMax;
	}
	// Initialize the iteration counter
	k = 0;
	// Set the initial state and function to minimize as the initial estimate
	y = y_ini;
	F = F_ini;
	dy = y;
	// Loop to solve the nonlinear system
	// while ((arma::norm(F, "inf") > atol) & (k <= NN))
	// while (((arma::norm(dy, "inf") > atol + rtol * arma::norm(y, "inf")) | (arma::norm(F, "inf") > atol)) & (k <= NN))
	while (((arma::norm(dy, "inf") > atol) | (arma::norm(F, "inf") > atol)) & (k <= NN))
	{
		// If the jacobean matrix was recycled twice, compute it again
		if (q > 2)
		{
			jac(t + dt, y);
			iJ = iJ + 1;
		}
		// Compute the nonlinear system matrix
		if (adaptivity)
		{
			tmp_coef = 0.0;
			for (int jj = 0; jj < N; jj = jj + 1)
			{
				tmp_coef = tmp_coef + 1.0 / ((t + dt) - tau_i(jj));
			}
			JF = I - Jfun / tmp_coef;
		}
		else
		{
			JF = I - f_coef * Jfun;
		}
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
		if (q <= 2)
		{
			if (adaptivity)
			{
				// If the jacobean matrix was not updated before, reduce the time step and try again
				if (dt > 100 * dt_min)
				{
					dt = std::max(dt / 100, dt_min);
					std::cout << "WARNING: In BDFN, Maximum number of iterations reached! Trying again with smaller step size... " << std::endl;
				}
				else
				{
					if (q < 2)
					{
						jac(t, y_i.col(N - 1));
						iJ = iJ + 1;
						q = 1;
					}
					dt = dt_min;
					std::cout << "WARNING: In BDFN, Maximum number of iterations reached! " << std::endl;
					std::cout << "         Trying again with minimum step size and Jacobean matrix recomputation... " << std::endl;
				}
			}
			else
			{
				// Otherwise, compute the jacobean matrix again and try again, first with the current time step and then with the minimum time step
				if (q == 0)
				{
					jac(t, y_i.col(N - 1));
					iJ = iJ + 1;
					std::cout << "WARNING: In BDFN, Maximum number of iterations reached! Trying again with recomputed Jacobean matrix... " << std::endl;
				}
				else if (q == 1)
				{
					dt = dt_min;
					std::cout << "WARNING: In BDFN, Maximum number of iterations reached! Trying again with minimum time step..." << std::endl;
				}
				else if (q == 2)
				{
					std::cout << "WARNING: In BDFN, Maximum number of iterations reached! Trying again recomputing the Jacobean matrix at every iteration..." << std::endl;
				}
			}
			q = q + 1;
			goto LOOP;
		}
		else
		{
			// Otherwise, return an error
			std::cout << "ERROR: Convergence Failed!" << std::endl;
			throw std::exception();
		}
	}

	// Check for NaN values in the state vector
	if (y.has_nan())
	{
		// If NaN is detected for the first time, try again with a smaller step size, otherwise return an error
		if (flag_nan)
		{
			// update the NaN flag
			flag_nan = false;
			// reduce the time step
			dt = dt_min;
			// recompute the jacobean matrix
			y = y_ini;
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

	// If the solver uses adaptivity, compute the error weighted tolerance and the local truncation error
	// and the adaptive time step multiplier, otherwise set the multiplier to 1.1 so the time step increases
	// slowly from the initial value to the maximum value at the start of the simulation.
	// TODO: sigma_min and sigma_max should be parameters
	double sigma = 1.05;
	// if (adaptivity)
	// {
	// 	// Compute the error weighted tolerance
	// 	EWT = atol * arma::ones(size(y)) + rtol * arma::abs(y);
	// 	// Compute the local truncation error comparing with a BDF(N-1) scheme
	// 	LTE = arma::zeros(size(y));
	// 	arma::mat y_LTE = arma::join_horiz(arma::join_horiz(y_prev, y_i), y);
	// 	arma::vec tau_LTE = arma::join_vert(arma::join_vert(tau_prev * arma::ones(1), tau_i), (t + dt) * arma::ones(1));
	// 	for (int ii = 0; ii < N + 2; ii = ii + 1)
	// 	{
	// 		tmp_coef = 1.0;
	// 		for (int jj = 0; jj < N + 2; jj = jj + 1)
	// 		{
	// 			if (jj != ii)
	// 			{
	// 				tmp_coef = tmp_coef / (tau_LTE(ii) - tau_LTE(jj));
	// 			}
	// 		}
	// 		LTE = LTE + y_LTE.col(ii) * tmp_coef;
	// 	}
	// 	LTE = LTE * pow(dt, N + 1);
	// 	// Compute the adaptive time step multiplier
	// 	sigma = pow(0.5 * arma::norm(EWT) / arma::norm(LTE), 0.25);
	// 	// If the error is too large, reduce the time step and try again
	// 	if ((sigma < 0.9) & (dt > dt_min))
	// 	{
	// 		dt = std::max(dt * sigma, dt_min);
	// 		std::cout << "	WARNING: In BDFN, Error too large! Trying again with smaller step size... " << std::endl;
	// 		goto LOOP;
	// 	}
	// 	// Set an upper bound for sigma
	// 	sigma = std::min(sigma, 2.0);
	// }

	// Update the time, the last three states and time steps, and the next time step size
	if (adaptivity)
	{
		t = t + dt;
	}
	else
	{
		t = dt * round((t + dt) / dt);
	}
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
	if (adaptivity)
	{
		dt = sigma * dt;
		dt = std::max(dt, dt_min);
		dt = std::min(dt, dt_max);
		dt = std::min(dt, dt_out - std::fmod(t, dt_out) + 1e-15);
		dt = std::min(dt, tmax - t + 1e-15);
	}
}

// Methods for the ESDIRK class
ESDIRK::ESDIRK(bool a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u, ISimulation *pIncSim)
{
	// ***ESDIRK constructor***
	// Set the simulation pointer
	pSim = pIncSim;
	// Set the initial time, maximum time, output time step and state vector
	t = t_u;
	tmax = tmax_u;
	dt_out = dt_out_u;
	y = y_u;
	adaptivity = a_u;
	dt_max = dt_max_u;
	// Compute the system size
	nSystem = y_u.n_rows;
	// Initialize identity matrix and jacobian matrix
	I = arma::eye(nSystem, nSystem);
	Jfun = arma::zeros(nSystem, nSystem);
	// Initialize the required vectors and matrices
	ys = arma::zeros(nSystem, s);
	ts = arma::zeros(s);
	// Initialize time step size
	if (adaptivity)
	{
		dt = dt_ini;
	}
	else
	{
		dt = dt_max;
	}
	dt_old = dt;
	// Compute the jacobi matrix at the initial state
	jac(t, y);
}

void ESDIRK::init(void)
{
	std::cout << "WARNING: Initialization of ESDIRK not needed!" << std::endl;
}

arma::mat ESDIRK::fun(double tt, arma::mat yy)
{
	// ***Interface method to the simulation method responsible for the state derivative evaluation***
	return pSim->CalculateSystemDynamics(tt, yy);
}

void ESDIRK::jac(double tt, arma::mat yy)
{
	// ***Compute the jacobian matrix of the function that returns the state derivative***
	arma::mat yprime = fun(tt, yy);
	for (int ii = 0; ii < nSystem; ii = ii + 1)
	{
		Jfun.col(ii) = 1e12 * (fun(tt, yy + 1e-12 * I.col(ii)) - yprime);
	}
}

arma::mat ESDIRK::ESDIRK_fun(double t_i, arma::mat y_i, int ii)
{
	//***Evaluate the nonlinear ESDIRK scheme function***
	arma::mat FF = y_i - y - fun(t_i, y_i) * dt * a(ii, ii);
	for (int jj = 0; jj < ii; jj = jj + 1)
	{
		FF = FF - dt * a(ii, jj) * fun(ts(jj), ys.col(jj));
	}
	return FF;
}

void ESDIRK::step(void)
{
	// ***Perform a time step with the ESDIRK scheme***
	// Check if the jacobi matrix needs to be updated
	bool jac_update = false;
	if (nSteps > nStepsMax)
	{
		jac(t, y);
		iJ = iJ + 1;
		nSteps = 0;
		jac_update = true;
	}
	else
	{
		nSteps = nSteps + 1;
	}
	// Declare the stage state vector
	arma::mat y_i;
LOOP:
	// Loop over the stages
	for (int ii = 0; ii < s; ii = ii + 1)
	{
		// Initialize the stage state vector
		y_i = y;
		// Initialize the iteration counter
		int k = 0;
		// Initialize the stage state change vector
		arma::mat dy = 2.0 * atol * arma::ones(size(y));
		// Define the matrix for the implicit stage
		arma::mat MM = I + dt * a(ii, ii) * Jfun;
		// Iterative scheme to solve the stage
		while ((arma::norm(dy, "inf") > atol) & (k <= nIterMax))
		{
			// Define the vector for the implicit stage for the current iteration
			arma::mat VV = ESDIRK_fun(t + beta(ii) * dt, y_i, ii);
			// Solve for the state change vector
			status = arma::solve(dy, MM, VV, arma::solve_opts::fast);
			// If the fast solver fails, solve the system using the slower method
			if (!status)
			{
				dy = arma::solve(MM, VV);
			}
			// Update the state vector
			y_i = y_i + dy;
			// Update the iteration counter
			k = k + 1;
		}
		// Check if the convergence failed
		if ((k >= nIterMax) || (y_i.has_nan()))
		{
			if (jac_update)
			{
				std::cout << "ERROR: Convergence failed in ESDIRK!" << std::endl;
				throw std::exception();
			}
			else
			{
				jac(t, y);
				iJ = iJ + 1;
				jac_update = true;
				nSteps = 0;
				std::cout << "WARNING: In ESDIRK, Maximum number of iterations reached! Trying again with recomputed Jacobean matrix... " << std::endl;
				goto LOOP;
			}
		}
		// Save the stage state vector
		ys.col(ii) = y_i;
		// Save the stage time vector
		ts(ii) = t + beta(ii) * dt;
	}
	// Update the state vector
	y = y_i;
	// Update the time
	t = t + dt;
	if (adaptivity)
	{
		// Compute the embeded solution
		yhat = arma::zeros(size(y));
		for (int ii = 0; ii < s; ii = ii + 1)
		{
			yhat = yhat + ys.col(ii) * b(ii);
		}
		// Compute the LTE error
		double LTE = arma::norm(yhat - y, "inf");
		// Compute the EWT error
		double EWT = atol + rtol * arma::norm(y, "inf");
		// Compute the adaptive time step multiplier
		error_ratio_old = error_ratio;
		error_ratio = LTE / EWT;
		// Compute rho
		double rho = pow(dt / dt_old, -1.0 / 4.0) * pow(error_ratio * error_ratio_old, -1.0 / 12.0);
		// Update the time step size
		dt_old = dt;
		dt = (1.0 + 2.0 * std::atan((rho - 1.0) / 2.0)) * dt;
		dt = std::max(dt, dt_min);
		dt = std::min(dt, dt_max);
		dt = std::min(dt, dt_out - std::fmod(t, dt_out) + 1e-15);
		dt = std::min(dt, tmax - t + 1e-15);
	}
}