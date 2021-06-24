
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

BDF::BDF(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation* pIncSim)
{
	pSim = pIncSim;

	t = t_u;
	tmax = tmax_u;
	dt_out = dt_out_u;
	y = y_u;

	h_0 = dt_ini;
	h_1 = dt_ini;
	h_2 = dt_ini;
	dt_max = std::min(0.5*dt_out,dt_max);

	nSistema = y_u.n_rows;
	y_0 = y;
	y_1 = y;
	y_2 = y;
	yprime = arma::zeros(size(y));
	this->F = arma::zeros(size(y));
	LTE = arma::zeros(size(y));
	EWT = arma::zeros(size(y));
	I = arma::eye(nSistema,nSistema);
	J = arma::zeros(nSistema,nSistema);
	
}


arma::mat BDF::fun(double tt, arma::mat yy)
{
	return pSim->CalculateSystemDynamics(tt, yy);
}


void BDF::jac(double tt, arma::mat yy){

	//std::cout << std::endl << " JACOBEAN : START" << std::endl;

	yprime = fun(tt, yy);
	for(int ii=0;ii<nSistema;ii=ii+1){
		J.col(ii) = 1e12 * (fun(tt, yy + 1e-12 * I.col(ii)) - yprime);
	}

	//std::cout << "              Determinant      : " << arma::det(J) << std::endl;
	//std::cout << "              Condition Number : " << arma::cond(J) << std::endl;
	//std::cout << "          : FINISH" << std::endl << std::endl;

}


void BDF::Initialize()
{
	F(0) = 2*atol;
	do{
		std::cout << "    Computing jac... " << std::endl;
		jac(t + h_0, y);
		std::cout << "    ... done! " << std::endl;
		M = I - h_0 * J;
		F = y - y_0 - h_0 * yprime;
		std::cout << "    size dy: " << arma::size(dy) << std::endl;
		std::cout << "    size M: " << arma::size(M) << std::endl;
		std::cout << "    size F: " << arma::size(F) << std::endl;
		status = arma::solve(dy,M,F,arma::solve_opts::fast);
		std::cout << "    After solver..." << std::endl;
		if (!status){
			dy = arma::solve(M,F);
		}
		y = y - dy;
	} while((arma::norm(F,"inf") > atol));
	y_0 = y;
	t = t + h_0;

	jac(t + h_0, y_0);

}

void BDF::step(void){
	
	int NN;
	q = 0;
	bool flag_nan = true;

	LOOP:

	if (q<2){
		NN = 5*(q+1);
	} else {
		NN = nIterMax;
	}

	k = 0;
	y = y_0 + h_0*(y_0 - y_1)/h_1;
	F = (1.0 + h_0/(h_1+h_0)) * y - ((h_1+h_0)/h_1) * y_0 + ((h_0*h_0/h_1)/(h_1+h_0)) * y_1 - h_0 * fun(t+h_0,y);

	do{
		if (q>=2){
			jac(t + h_0, y);
			iJ = iJ + 1;
			q = q + 1;
		}	
		M = (1.0 + h_0/(h_1+h_0)) * I - h_0 * J;
		status = arma::solve(dy,M,F,arma::solve_opts::fast);
		if (!status){
			dy = arma::solve(M,F);
		}
		y = y - dy;
		F = (1.0 + h_0/(h_1+h_0)) * y - ((h_1+h_0)/h_1) * y_0 + ((h_0*h_0/h_1)/(h_1+h_0)) * y_1 - h_0 * fun(t+h_0,y);
		k = k + 1;
	} while(((arma::norm(dy) > atol + rtol*arma::norm(y)) | (arma::norm(F) > atol)) & (k<NN));

	if(k>=NN){
		if(q<2){
			h_0 = std::max(pow(10.0,-2*q)*h_0, dt_min);
			jac(t + h_0, y_0 + h_0*(y_0 - y_1)/h_1);
			iJ = iJ + 1;
			q = q + 1;
			goto LOOP;
		} else {
			std::cout << std::endl << "ERROR: Convergence Failed!" << std::endl;
			throw std::exception();
		}
	}

	EWT = atol * arma::ones(size(y)) + rtol * arma::abs(y);
	//EWT = atol  + rtol % arma::abs(y);

	if (EWT.has_nan()){		
		if (flag_nan){
			flag_nan = false;
			h_0 = dt_min;
			jac(t + h_0, y_0 + h_0*(y_0 - y_1)/h_1);
			iJ = iJ + 1;
			std::cout << std::endl << " OJO QUE ESTO CASCA ..." << std::endl;
			goto LOOP;
		}else{
			std::cout << std::endl << "ERROR: NaN Detected!" << std::endl;
			throw std::exception();
		}
	}
	
	//LTE = h_0*h_0*(h_0+h_1)*(y/(h_0*(h_0 + h_1)*(h_0 + h_1 + h_2)) - y_2/(h_2*(h_1 + h_2)*(h_0 + h_1 + h_2)) - y_0/(h_0*h_1*(h_1 + h_2)) + y_1/(h_1*h_2*(h_0 + h_1)));
	//sigma = pow(0.5*arma::abs(EWT/LTE).min(),0.25);

	arma::mat M1 = I - h_0 * J;
	arma::mat F1 = y - y_0 - h_0 * fun(t+h_0,y);
	status = arma::solve(LTE,M1,F1,arma::solve_opts::fast);
	if (!status){
		LTE = arma::solve(M1,F1);
	}

	sigma = pow(0.5*arma::norm(EWT)/arma::norm(LTE),0.25);

	if(sigma<0.9){
		h_0 = h_0*sigma;
		goto LOOP;
	}

	//std::cout << "   LTE = " << arma::norm(LTE) << std::endl;
	//std::cout << "   h = " << h_0 << std::endl;

	t = t + h_0;

	h_2 = h_1;
	h_1 = h_0;
	h_0 = sigma * h_0;
	h_0 = std::max(h_0, dt_min);
	h_0 = std::min(h_0, dt_max);
	h_0 = std::min(h_0, dt_out - std::fmod(t,dt_out) + dt_min);
	h_0 = std::min(h_0, tmax - t + h_0);

	y_2 = y_1;
	y_1 = y_0;
	y_0 = y;
	
}