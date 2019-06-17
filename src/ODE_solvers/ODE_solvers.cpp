
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

BDF::BDF(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, arma::mat (*fun_u) (double, arma::mat, solver_data), solver_data SD_u)
{
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
	printf("Before bucle....\n");
	do{
		printf("Before jac\n");
		jac(t + h_0, y);
		M = I - h_0 * J;
		F = y - y_0 - h_0 * yprime;
		printf("Before solve\n");
		std::cout << "size dy: " << arma::size(dy) << std::endl;
		std::cout << "size M: " << arma::size(M) << std::endl;
		std::cout << "size F: " << arma::size(F) << std::endl;
		status = arma::solve(dy,M,F,arma::solve_opts::fast);
		printf("ddddddd\n");
		if (!status){
			dy = arma::solve(M,F);
		}
		y = y - dy;
	} while((arma::norm(F,"inf") > atol));
	y_0 = y;
	t = t + h_0;

	jac(t + h_0, y_0);
}


void BDF::jac(double tt, arma::mat yy){
	yprime = fun(tt, yy, SD);
	printf("After fun...\n");
	for(int ii=0;ii<SD.nSistema;ii=ii+1){
		J.col(ii) = 1e12 * (fun(tt, yy + 1e-12* I.col(ii), SD) - yprime);
	}
}


void BDF::step(void){

	q = 0;

	LOOP:
	k = 0;
	y = y_0 + h_0*(y_0 - y_1)/h_1;
	do{
		F = (1.0 + h_0/(h_1+h_0)) * y - ((h_1+h_0)/h_1) * y_0 + ((h_0*h_0/h_1)/(h_1+h_0)) * y_1 - h_0 * fun(t+h_0,y,SD);
		M = (1.0 + h_0/(h_1+h_0)) * I - h_0 * J;
		status = arma::solve(dy,M,F,arma::solve_opts::fast);
		if (!status){
			dy = arma::solve(M,F);
		}
		y = y - dy;
		k = k + 1;
	} while((arma::norm(dy) > atol + rtol*arma::norm(y))&&(k<=nIterMax+1));

	if(k>=nIterMax-1){
		if((q==0)){
			jac(t + h_0, y_0 + h_0*(y_0 - y_1)/h_1);
			iJ = iJ + 1;
		} else {
			h_0 = std::max(1e-1*h_0, dt_min);
			jac(t + h_0, y_0 + h_0*(y_0 - y_1)/h_1);
			iJ = iJ + 1;
		}
		q = q + 1;
		goto LOOP;
	}

	t = t + h_0;
	
	LTE = h_0*h_0*(h_0+h_1)*(y/(h_0*(h_0 + h_1)*(h_0 + h_1 + h_2)) - y_2/(h_2*(h_1 + h_2)*(h_0 + h_1 + h_2)) - y_0/(h_0*h_1*(h_1 + h_2)) + y_1/(h_1*h_2*(h_0 + h_1)));
	EWT = rtol * arma::abs(y) + atol * arma::ones(size(y));
	sigma = pow(arma::norm(LTE/EWT)/sqrt(nSistema),-1.0/3.0);
	if(isnan(sigma)){
		sigma = 1e-12;
	}

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