
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

void BDF::jac(double tt, arma::mat yy){
	yprime = fun(tt, yy, SD);
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
	h_0 = std::min(h_0, tmax - t);

	y_2 = y_1;
	y_1 = y_0;
	y_0 = y;
}