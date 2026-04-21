// SPDX-License-Identifier: GPL-3.0-or-later
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
#include "../os_tools.hpp"
#include "../Logger.hpp"

// Methods for the BDF2 class
BDF2::BDF2(double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u, ISimulation* pIncSim)
{
    // ***BDF2 constructor***
    // Set the simulation pointer
    pSim = pIncSim;
    // Set the initial time, maximum time, output time step and state vector
    t = t_u;
    tmax = tmax_u;
    dt_out = dt_out_u;
    y = y_u;
    // Initialize the variables with the last three time steps
    h_0 = dt_ini;
    h_1 = dt_ini;
    h_2 = dt_ini;
    // Compute the maximum time step
    dt_max = std::min(dt_out, dt_max_u);
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

    // Set the debug flag
    debug_flag = false;
    if (debug_flag)
    {
        Logger::debug("BDF2 debug file opened!");
        // Open the debug file
        char buffer[50];
        int nn = sprintf(buffer, "time_adaptivity_debug.txt");
        std::string file_path = JoinPath(pSim->outputFolderPath, buffer);
        pfile = fopen(file_path.c_str(), "w");
        // Print the header
        fprintf(pfile, "t    dt    EWT    LTE    rho    rejected\n");
    }
}

arma::mat BDF2::fun(double tt, arma::mat yy)
{
    // ***Interface method to the simulation method responsible for the state derivarive evaluation***
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
    // Initialize the function to minimize to a value larger than the tolerance, and the iteration counter
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
        Logger::error("Convergence Failed!");
        throw std::exception();
    }
    // Otherwise, asign the state computed with the BDF1 scheme to the last state and update the solver time
    y_0 = y;
    t = t + h_0;
    // Compute the jacobian matrix at the new state, so it can be reused in the next iterations
    jac(t, y);
}

void BDF2::finalize(void)
{
    if (debug_flag)
    {
        Logger::debug("BDF2 debug file closed!");
        fclose(pfile);
    }
    nNewtonIterAvg = nNewtonIter / nNewton;
}

arma::mat BDF2::BDF2_fun(double t, arma::mat y)
{
    //***Evaluate the nonlinear BDF2 scheme function***
    F = (1.0 + h_0 / (h_1 + h_0)) * y - ((h_1 + h_0) / h_1) * y_0 + ((h_0 * h_0 / h_1) / (h_1 + h_0)) * y_1 -
        h_0 * fun(t + h_0, y);
    return F;
}

void BDF2::step(void)
{
    // ***Perform a time step with the BDF2 scheme***
    int NN;               // Declare local maximum number of iterations
    q = 0;                // Initialize jacobian recycling counter
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
    // Update the number of times the newton method is used
    nNewton = nNewton + 1;
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
    // Update the number of newton iterations
    nNewtonIter = nNewtonIter + k;
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
            // std::cout << "WARNING: In BDF2, Maximum number of iterations reached! Trying again with smaller step
            // size... " << std::endl;
            goto LOOP;
        }
        else
        {
            // Otherwise, return an error
            Logger::error("Convergence Failed!");
            throw std::exception();
        }
    }
    // Set the tolerance values for the error weighted tolerance
    double atol_EWT = 1e-6;
    double rtol_EWT = 1e-3;
    // Compute the error weighted tolerance
    EWT = atol_EWT * arma::ones(size(y)) + rtol_EWT * arma::abs(y);
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
            Logger::warning("In BDF2, NaN detected! Trying again with smaller step size...");
            goto LOOP;
        }
        else
        {
            Logger::error("NaN detected after trying again!");
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
    if (debug_flag)
    {
        // Print the debug information to the file
        fprintf(pfile, "%e    ", t);
        fprintf(pfile, "%e    ", h_0);
        fprintf(pfile, "%e    ", arma::norm(EWT));
        fprintf(pfile, "%e    ", arma::norm(LTE));
        fprintf(pfile, "%e    ", sigma);
    }
    // If the error is too large, reduce the time step and try again
    if ((sigma < 0.9) & (h_0 > dt_min))
    {
        if (debug_flag)
            fprintf(pfile, "%d\n", 1);
        h_0 = h_0 * sigma;
        Logger::warning("In BDF2, Error too large! Trying again with smaller step size...");
        goto LOOP;
    }
    if (debug_flag)
        fprintf(pfile, "%d\n", 0);
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
    if (!debug_flag)
        h_0 = std::min(h_0, dt_out - std::fmod(t, dt_out) + dt_min);
    h_0 = std::min(h_0, tmax - t + h_0);
}

// Methods for the BDFN class
BDFN::BDFN(int N_u, int a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, arma::mat y_u,
           ISimulation* pIncSim)
{
    // ***BDFN constructor***
    // Check for the order of the BDF scheme to be larger than 2
    if (N_u < 1)
    {
        Logger::error("BDFN order must be larger or equal than 1!");
        throw std::exception();
    }
    // Check for the order of the BDF scheme to be smaller than 7
    if (N_u > 6)
    {
        Logger::error("BDFN order must be smaller or equal than 6!");
        throw std::exception();
    }
    // Check for the output time step to be larger than the maximum time step
    if (dt_out_u < dt_max_u)
    {
        Logger::error("Output time step must be larger or equal than the maximum time step!");
        throw std::exception();
    }
    // If adaptivity is not used, check if the output time step is a multiple of the maximum time step
    if (!a_u & (std::fmod(dt_out_u, dt_max_u) > 1e-15) & (std::fmod(dt_out_u, dt_max_u) < dt_max_u - 1e-15))
    {
        Logger::error("Output time step must be a multiple of the maximum time step!");
        throw std::exception();
    }
    // If adaptivity is not used, check if the total simulation time is a multiple of the maximum time step
    if (!a_u & (std::fmod(tmax_u, dt_max_u) > 1e-14) & (std::fmod(tmax_u, dt_max_u) < dt_max_u - 1e-14))
    {
        Logger::error("Total simulation time must be a multiple of the maximum time step!"
                      " tmax_u: " + std::to_string(tmax_u) +
                      " dt_max_u: " + std::to_string(dt_max_u) +
                      " fmod: " + std::to_string(std::fmod(tmax_u, dt_max_u)));
        throw std::exception();
    }
    // Set the order of the BDF scheme
    N = N_u;
    // Set the adaptive time step flag
    if (a_u == 0)
    {
        adaptivity = false;
    }
    else
    {
        adaptivity = true;
    }
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
        Logger::error("BDFN order must be between 1 and 6!");
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
            Logger::error("Convergence failed in BDFN during init!");
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

void BDFN::finalize(void)
{
    // if (debug_flag)
    // {
    // 	std::cout << "DEBUG: BDFN debug file closed!" << std::endl;
    // 	fclose(pfile);
    // }
    nNewtonIterAvg = nNewtonIter / nNewton;
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
    int NN;               // Declare local maximum number of iterations
    double tmp_coef;      // Declare local variable for the temporal coefficient
    q = 0;                // Initialize jacobian recycling counter
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
        Logger::error("In BDFN, without adaptivity, it was not possible to try again with a smaller step size!");
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
    // Update the number of times the newton method is used
    nNewton = nNewton + 1;
    // Loop to solve the nonlinear system
    // while ((arma::norm(F, "inf") > atol) & (k <= NN))
    // while (((arma::norm(dy, "inf") > atol + rtol * arma::norm(y, "inf")) | (arma::norm(F, "inf") > atol)) & (k <=
    // NN))
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
    // Update the number of newton iterations
    nNewtonIter = nNewtonIter + k;
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
                    Logger::warning("In BDFN, Maximum number of iterations reached! Trying again with smaller step size...");
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
                    Logger::warning("In BDFN, Maximum number of iterations reached! Trying again with minimum step size and Jacobean matrix recomputation...");
                }
            }
            else
            {
                // Otherwise, compute the jacobean matrix again and try again, first with the current time step and then
                // with the minimum time step
                if (q == 0)
                {
                    jac(t, y_i.col(N - 1));
                    iJ = iJ + 1;
                    Logger::warning("In BDFN, Maximum number of iterations reached! Trying again with recomputed Jacobean matrix...");
                }
                else if (q == 1)
                {
                    dt = dt_min;
                    Logger::warning("In BDFN, Maximum number of iterations reached! Trying again with minimum time step...");
                }
                else if (q == 2)
                {
                    Logger::warning("In BDFN, Maximum number of iterations reached! Trying again recomputing the Jacobean matrix at every iteration...");
                }
            }
            q = q + 1;
            goto LOOP;
        }
        else
        {
            // Otherwise, return an error
            Logger::error("Convergence Failed!");
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
            Logger::warning("In BDF2, NaN detected! Trying again with smaller step size...");
            goto LOOP;
        }
        else
        {
            Logger::error("NaN detected after trying again!");
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
ESDIRK::ESDIRK(bool a_u, double t_u, double tmax_u, double dt_max_u, double dt_out_u, int nStepsMax_u, arma::mat y_u,
               ISimulation* pIncSim)
{
    // ***ESDIRK constructor***
    // Set the simulation pointer
    pSim = pIncSim;
    // Set the initial time, maximum time, output time step and state vector
    t = t_u;
    tmax = tmax_u;
    dt_out = dt_out_u;
    y = y_u;
    dt_max = dt_max_u;
    nStepsMax = nStepsMax_u;
    adaptivity = a_u;
    // Define the ESDIRK coefficients
    s = 6;
    beta = {0.0, 1.0 / 2.0, 83.0 / 250.0, 31.0 / 50.0, 17.0 / 20.0, 1.0};
    a = arma::zeros<arma::mat>(6, 6);
    a(1, 0) = 1.0 / 4.0;
    a(1, 1) = 1.0 / 4.0;
    a(2, 0) = 8611.0 / 62500.0;
    a(2, 1) = -1743.0 / 31250.0;
    a(2, 2) = 1.0 / 4.0;
    a(3, 0) = 5012029.0 / 34652500.0;
    a(3, 1) = -654441.0 / 2922500.0;
    a(3, 2) = 174375.0 / 388108.0;
    a(3, 3) = 1.0 / 4.0;
    a(4, 0) = 15267082809.0 / 155376265600.0;
    a(4, 1) = -71443401.0 / 120774400.0;
    a(4, 2) = 730878875.0 / 902184768.0;
    a(4, 3) = 2285395.0 / 8070912.0;
    a(4, 4) = 1.0 / 4.0;
    a(5, 0) = 82889.0 / 524892.0;
    a(5, 2) = 15625.0 / 83664.0;
    a(5, 3) = 69875.0 / 102672.0;
    a(5, 4) = -2260.0 / 8211.0;
    a(5, 5) = 1.0 / 4.0;
    b = {4586570599.0 / 29645900160.0, 0.0,
         178811875.0 / 945068544.0,    814220225.0 / 1159782912.0,
         -3700637.0 / 11593932.0,      61727.0 / 225920.0};
    // Compute the system size
    nSystem = y_u.n_rows;
    // Initialize identity matrix and jacobian matrix
    I = arma::eye(nSystem, nSystem);
    Jfun = arma::zeros(nSystem, nSystem);
    // Initialize the required vectors and matrices
    ys = arma::zeros(nSystem, s);
    ys_prime = arma::zeros(nSystem, s);
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
    dt_old2 = dt;
    // Compute the jacobi matrix at the initial state
    jac(t, y);

    // Set the debug flag
    debug_flag = false;
    if (debug_flag)
    {
        Logger::debug("ESDIRK debug file opened!");
        // Open the debug file
        char buffer[50];
        int nn = sprintf(buffer, "time_adaptivity_debug.txt");
        std::string file_path = JoinPath(pSim->outputFolderPath, buffer);
        pfile = fopen(file_path.c_str(), "w");
        // Print the header
        fprintf(pfile, "t    dt    EWT    LTE    rho    rejected\n");
    }

    if (adaptivity)
    {
        char buffer1[50];
        char bufferLine[1000];
        int nn1 = sprintf(buffer1, "ESDIRK_adaptive_time_step_data.dat");
        std::string file_path1 = JoinPath(pSim->inputFolderPath, buffer1);
        FILE* pfile1 = fopen(file_path1.c_str(), "r");
        if (pfile1 == NULL)
        {
            Logger::warning("ESDIRK time step adaptive data file not found! Using default values instead");
            adaptivity_type = 2;
            adaptivity_smooth_flag = 2;
            adaptivity_smooth_param = 2.0;
            LTE_norm_type = 2;
            rho_min = 2.0 / 3.0;
        }
        else
        {
            fscanf(pfile1, "%d %[^\n]\n", &adaptivity_type, bufferLine);
            fscanf(pfile1, "%d %[^\n]\n", &adaptivity_smooth_flag, bufferLine);
            fscanf(pfile1, "%lf %[^\n]\n", &adaptivity_smooth_param, bufferLine);
            fscanf(pfile1, "%d %[^\n]\n", &LTE_norm_type, bufferLine);
            fscanf(pfile1, "%lf %[^\n]\n", &rho_min, bufferLine);
            // close the file
            fclose(pfile1);
        }
    }
}

void ESDIRK::init(void)
{
    Logger::warning("Initialization of ESDIRK not needed!");
}

void ESDIRK::finalize(void)
{
    if (debug_flag)
    {
        Logger::debug("ESDIRK debug file closed!");
        fclose(pfile);
    }
    nNewtonIterAvg = nNewtonIter / nNewton;
}

arma::mat ESDIRK::fun(double tt, arma::mat yy)
{
    // ***Interface method to the simulation method responsible for the state derivative evaluation***
    return pSim->CalculateSystemDynamics(tt, yy);
}

void ESDIRK::jac(double tt, arma::mat yy)
{
    // ***Compute the jacobian matrix of the function that returns the state derivative***
    y_prime = fun(tt, yy);
    for (int ii = 0; ii < nSystem; ii = ii + 1)
    {
        Jfun.col(ii) = 1e12 * (fun(tt, yy + 1e-12 * I.col(ii)) - y_prime);
    }
    iJ = iJ + 1;
}

arma::mat ESDIRK::ESDIRK_fun(double t_i, arma::mat y_i, int ii)
{
    //***Evaluate the nonlinear ESDIRK scheme function***
    y_prime = fun(t_i, y_i);
    arma::mat FF = y_i - y - dt * a(ii, ii) * y_prime;
    for (int jj = 0; jj < ii; jj = jj + 1)
    {
        FF = FF - dt * a(ii, jj) * ys_prime.col(jj);
    }
    return FF;
}

void ESDIRK::step(void)
{
    // ***Perform a time step with the ESDIRK scheme***
    // Check if the jacobi matrix needs to be updated
    bool jac_updated = false;
    bool dt_reduced = false;
    if ((nSteps > nStepsMax) & (nStepsMax >= 0))
    {
        jac(t, y);
        nSteps = 0;
        jac_updated = true;
    }
    nSteps = nSteps + 1;
    // Declare the stage state vector
    arma::mat y_i;
    // Declare the first stage state vector
    ys.col(0) = y;
    // Declare the first stage state derivative vector
    ys_prime.col(0) = y_prime;
    // If adaptivity is not used, make sure that the current time is a multiple of the time step
    if (!adaptivity)
    {
        t = dt * round(t / dt);
    }
    // Declare the first stage time vector
    ts(0) = t;
LOOP:
    // Loop over the stages
    for (int ii = 1; ii < s; ii = ii + 1)
    {
        // Initialize the stage state vector
        y_i = y;
        // Define the time of the implicit stage
        double t_i = t + beta(ii) * dt;
        // Initialize the iteration counter
        int k = 0;
        // Initialize the stage state change vector
        arma::mat dy = 2.0 * atol * arma::ones(size(y));
        // Define the vector for the implicit stage for the current iteration
        arma::mat FF = ESDIRK_fun(t_i, y_i, ii);
        // Define the matrix for the implicit stage
        arma::mat MM = I - dt * a(ii, ii) * Jfun;
        // Update the number of times using the Newton method
        nNewton = nNewton + 1;
        // Iterative scheme to solve the stage
        while (((arma::norm(dy, "inf") > atol) || (arma::norm(FF, "inf") > atol)) & (k <= nIterMax))
        {
            // Solve for the state change vector
            status = arma::solve(dy, MM, -FF, arma::solve_opts::fast);
            // If the fast solver fails, solve the system using the slower method
            if (!status)
            {
                dy = arma::solve(MM, -FF);
            }
            // Compute the line-search step
            double alpha = 1e-4;
            double lambda = 0.9;
            double rho_min = 1e-2;
            double rho = 1.0;
            arma::mat F_new = ESDIRK_fun(t_i, y_i + rho * dy, ii);
            while ((arma::norm(F_new) > (1.0 - alpha) * arma::norm(FF)) & (rho > rho_min))
            {
                rho = rho * lambda;
                F_new = ESDIRK_fun(t_i, y_i + rho * dy, ii);
            }
            // Update the state change vector
            dy = rho * dy;
            // Update the stage state vector
            y_i = y_i + dy;
            // Redefine the function to minimize
            FF = F_new;
            // Update the iteration counter
            k = k + 1;
        }
        // Update the number of iterations
        nNewtonIter = nNewtonIter + k;
        // Check if the convergence failed
        if ((k >= nIterMax) || (y_i.has_nan()))
        {
            if (jac_updated)
            {
                if (adaptivity)
                {
                    if (dt_reduced)
                    {
                        Logger::error("Convergence failed in ESDIRK!");
                        throw std::exception();
                    }
                    else
                    {
                        dt = dt_ini;
                        dt_reduced = true;
                        Logger::warning("In ESDIRK, Maximum number of iterations reached! Trying again with minimum step size...");
                        goto LOOP;
                    }
                }
                else
                {
                    Logger::error("Convergence failed in ESDIRK!");
                    throw std::exception();
                }
            }
            else
            {
                nConvergenceFailed = nConvergenceFailed + 1;
                jac(t, y);
                jac_updated = true;
                nSteps = 0;
                Logger::warning("In ESDIRK, Maximum number of iterations reached! Trying again with recomputed Jacobean matrix...");
                goto LOOP;
            }
        }
        // Save the stage state vector
        ys.col(ii) = y_i;
        // Save the stage state derivative vector
        ys_prime.col(ii) = y_prime;
        // Save the stage time vector
        ts(ii) = t_i;
    }
    if (debug_flag)
    {
        fprintf(pfile, "%e    ", t);
        if (dt_truncated_flag)
        {
            fprintf(pfile, "%e    ", dt_truncated);
        }
        else
        {
            fprintf(pfile, "%e    ", dt);
        }
    }
    // Check if the time step size was too large
    if (adaptivity)
    {
        compute_dt();

        if (debug_flag)
        {
            fprintf(pfile, "%e    ", EWT);
            if (dt_truncated_flag)
            {
                fprintf(pfile, "%e    ", LTE_truncated);
            }
            else
            {
                fprintf(pfile, "%e    ", LTE);
            }
            fprintf(pfile, "%e    ", rho);
        }
        // Check if the error is too large to reject the time step
        if ((rho < rho_min) & (dt > dt_min))
        {
            set_dt();
            if (debug_flag)
            {
                fprintf(pfile, "%d\n", 1);
            }
            goto LOOP;
        }
        else
        {
            if (debug_flag)
            {
                fprintf(pfile, "%d\n", 0);
            }
        }
    }
    else
    {
        if (debug_flag)
        {
            fprintf(pfile, "%e    ", 0.0);
            fprintf(pfile, "%e    ", 0.0);
            fprintf(pfile, "%e    ", 1.0);
            fprintf(pfile, "%d\n", 0);
        }
    }
    // Update the time
    if (adaptivity)
    {
        // Update the time with the current time step size
        t = t + dt;
        // Update the time step size for the next iteration
        if (!dt_truncated_flag)
        {
            dt_old2 = dt_old;
            dt_old = dt;
        }
        set_dt();
    }
    else
    {
        // Update the time with the current time step size and round it
        // to the nearest multiple of the time step size to avoid floating point errors
        t = dt * round((t + dt) / dt);
    }
    // Update the state vector
    y = y_i;
}

void ESDIRK::set_dt(void)
{
    // ***Computes and sets the next time step size***
    if (dt_truncated_flag)
    {
        dt = dt_truncated;
        LTE = LTE_truncated;
        dt_truncated_flag = false;
    }
    else
    {
        // Compute the time step size based on the multiplier
        dt = dt * rho;

        // Set the minimum and maximum time step sizes, making sure that the time step size does not exceed the maximum
        // time
        dt = std::min(dt, tmax - t + 1e-12);
        dt = std::min(dt, dt_max);
        dt = std::max(dt, dt_min);

        // Truncate the time step size to the output time step size
        if ((std::fmod(t, dt_out) > 1e-12) && (!debug_flag))
        {
            double tmp_dt_out = dt_out * std::ceil(t / dt_out) - t;
            if ((dt > tmp_dt_out) && (tmp_dt_out >= dt_min))
            {
                dt_truncated_flag = true;
                dt_truncated = dt;
                LTE_truncated = LTE;
                dt = tmp_dt_out;
            }
        }
    }
}

void ESDIRK::compute_dt(void)
{
    // ***Compute LTE, EWT, error ratio and time step multiplier***

    // Set the tolerance values for the error weighted tolerance
    double atol_EWT = 1e-6;
    double rtol_EWT = 1e-3;

    // Compute the embeded solution
    yhat = arma::zeros(size(y));
    for (int ii = 0; ii < s; ii = ii + 1)
    {
        yhat = yhat + ys.col(ii) * b(ii);
    }

    // Compute the error weighted tolerance
    error_ratio_old2 = error_ratio_old;
    error_ratio_old = error_ratio;
    arma::mat error_LTE = y - yhat;
    arma::mat error_EWT = atol_EWT * arma::ones(size(y)) + rtol_EWT * arma::abs(y);
    if (LTE_norm_type == 1)
    {
        LTE = arma::norm(error_LTE, 1);
        EWT = arma::norm(error_EWT, 1);
        error_ratio = LTE / EWT;
    }
    else if (LTE_norm_type == 2)
    {
        LTE = arma::norm(error_LTE, 2);
        EWT = arma::norm(error_EWT, 2);
        error_ratio = LTE / EWT;
    }
    else if (LTE_norm_type == 3)
    {
        LTE = arma::norm(error_LTE, "inf");
        EWT = arma::norm(error_EWT, "inf");
        error_ratio = LTE / EWT;
    }
    else if (LTE_norm_type == 4)
    {
        error_ratio = 0.0;
        LTE = 0.0;
        EWT = 0.0;
        for (int ii = 0; ii < nSystem; ii = ii + 1)
        {
            double tmp_LTE = arma::as_scalar(error_LTE(ii));
            double tmp_EWT = arma::as_scalar(error_EWT(ii));
            LTE = LTE + pow(tmp_LTE, 2);
            EWT = EWT + pow(tmp_EWT, 2);
            error_ratio = error_ratio + pow(tmp_LTE / tmp_EWT, 2);
        }
        error_ratio = sqrt(error_ratio / nSystem);
        LTE = sqrt(LTE / nSystem);
        EWT = sqrt(EWT / nSystem);
    }
    else
    {
        Logger::error("Invalid LTE norm type!");
        throw std::exception();
    }

    // Define the time step size multiplier limits
    double rho_min = 1.0e-1;
    double rho_max = 1.0e+1;

    // Compute the time step multiplier
    if (abs(error_ratio) < 1e-12)
    {
        rho = rho_max;
    }
    else
    {
        if (adaptivity_type == 1)
        {
            rho = pow(error_ratio, -1.0 / 3.0);
        }
        else if (adaptivity_type == 2)
        {
            rho = pow(dt_old / dt, 1.0 / 4.0) * pow(error_ratio, -1.0 / 12.0) * pow(error_ratio_old, -1.0 / 12.0);
        }
        else if (adaptivity_type == 3)
        {
            rho = pow(dt_old / dt, 3.0 / 8.0) * pow(dt_old2 / dt_old, 1.0 / 8.0) * pow(error_ratio, -1.0 / 24.0) *
                  pow(error_ratio_old, -1.0 / 12.0) * pow(error_ratio_old2, -1.0 / 12.0);
        }
    }

    if (adaptivity_smooth_flag == 1)
    {
        rho = std::max(rho, rho_min);
        rho = std::min(rho, rho_max);
    }
    else if (adaptivity_smooth_flag == 2)
    {
        rho = (1.0 + adaptivity_smooth_param * atan((rho - 1.0) / adaptivity_smooth_param));
    }
    else
    {
        Logger::error("Invalid adaptivity smooth flag!");
        throw std::exception();
    }
}
