
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
	t_prev = t;
	y_prev = y;
	int p = 0;
	F = y - y_prev - dt * fun(t+dt, y, SD);
	M = arma::eye(nSistema,nSistema) - dt * jac(t+dt, y, SD);
	dy = arma::solve(M,F, arma::solve_opts::fast);
	do{
		//y0 = y;
		//F0 = F;
		y = y - dy;
		F = y - y_prev - dt * fun(t+dt, y, SD);
		M = arma::eye(nSistema,nSistema) - dt * jac(t+dt, y, SD);
		//yy = F - F0;
		//ss = y - y0;
		//M = M + ((yy-M*ss) * ss.t()) / arma::as_scalar((ss.t()) * ss);
		dy = arma::solve(M,F, arma::solve_opts::fast);
		p = p + 1;
	} while ((arma::norm(F,2)/nSistema > atol)&&(p<nIterMax));
	if (y.has_nan()){
		throw std::invalid_argument( "NaN in temporal integration" );
	}
	t = t + dt;
	eta = arma::norm(y-y_prev) / (arma::norm(y_prev) + atol);
	get_next_dt();

}


void BDF::get_next_dt(void){

	double time;

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