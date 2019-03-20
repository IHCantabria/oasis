
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

void BDF::step(void){
	if(t > 1e-12){
		y0 = y_prev;
	}
	t_prev = t;
	y_prev = y;
	int p = 0;
	jac(t+dt, y, yprime, J, SD);
	if(t < 1e-12){
		F = y - y_prev - dt * yprime;
		M = arma::eye(nSistema,nSistema) - dt * J;
	}else{
		F = y - (4.0/3.0) * y_prev + (1.0/3.0) * y0 - (2.0/3.0) * dt * yprime;
		M = arma::eye(nSistema,nSistema) - (2.0/3.0) * dt * J;
	}
	status = arma::solve(dy,M,F,arma::solve_opts::fast);
	if (!status){
		dy = arma::solve(M,F);
	}
	do{
		y = y - dy;
		jac(t+dt, y, yprime, J, SD);
		if(t < 1e-12){
			F = y - y_prev - dt * yprime;
			M = arma::eye(nSistema,nSistema) - dt * J;
		}else{
			F = y - (4.0/3.0) * y_prev + (1.0/3.0) * y0 - (2.0/3.0) * dt * yprime;
			M = arma::eye(nSistema,nSistema) - (2.0/3.0) * dt * J;
		}
		status = arma::solve(dy,M,F,arma::solve_opts::fast);
		if (!status){
			dy = arma::solve(M,F);
		}
		p = p + 1;
	} while ((arma::norm(F,2) > atol)&&(p<nIterMax));
	if (y.has_nan()){
		throw std::invalid_argument( "NaN in temporal integration" );
	}
	t = t + dt;
	get_next_dt();

}


void BDF::get_next_dt(void){

	double time;

	eta = arma::norm(y-y_prev) / (arma::norm(y_prev) + atol);

	if ((eta_min<=eta)&&(eta<=eta_max)){
		time = t + dt;
		if (time > tmax){
			dt = tmax - t;
		}
	}else if (eta<eta_min){
		dt = rho * dt;
		dt = std::min(dt,dt_max);
		dt = std::max(dt,dt_min);
		time = t + dt;
		if (time > tmax){
			dt = tmax - t;
		}
	}else if (eta>eta_max){
		dt = sigma * dt;
		dt = std::min(dt,dt_max);
		dt = std::max(dt,dt_min);
		y = y_prev;
		t = t_prev;
		step();
	}
}