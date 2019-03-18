#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include "classes.h"

/*
void BDF::step(void){
	arma::mat y0 = y;
	arma::mat M = arma::eye(nSistema,nSistema) - dt * jac(t+dt, y, SD);
	arma::mat F = y - y0 - dt * fun(t+dt, y, SD);
	arma::mat dy = arma::solve(M,F);
	do{
		y = y - dy;
		M = arma::eye(nSistema,nSistema) - dt * jac(t+dt, y, SD);
		F = y - y0 - dt * fun(t+dt, y, SD);
		dy = arma::solve(M,F);
	} while (arma::norm(F,2)/nSistema > atol);
	t = t + dt;
}
*/

/*
void BDF::step(void){
	t_prev = t;
	y_prev = y;
	y = y_prev + dt*fun(t+dt, y, SD);
	arma::mat y_old = y;
	double err = 2*atol;
	int p = 0;
	do{
		std::cout<< "Hola" << std::endl;
		y = y_prev + dt*fun(t+dt, y, SD);
		err = arma::norm(y_old-y,2);
		y_old = y;
		p = p + 1;
	} while ((err>atol)&&(p<nIterMax));
	if (p>=nIterMax){
		throw std::invalid_argument( "Reached nIterMax in predictor-corrector scheme." );
	}
	if (y.has_nan()){
		throw std::invalid_argument( "NaN in temporal integration" );
	}
	t = t + dt;
	eta = arma::norm(y-y_prev) / (arma::norm(y_prev) + atol);
	get_next_dt();

}
*/


void BDF::step(void){
	t_prev = t;
	y_prev = y;
	int p = 0;
	arma::mat y0 = y;
	arma::mat M = arma::eye(nSistema,nSistema) - dt * jac(t+dt, y, SD);
	arma::mat F = y - y0 - dt * fun(t+dt, y, SD);
	arma::mat dy = arma::solve(M,F);
	do{
		y = y - dy;
		M = arma::eye(nSistema,nSistema) - dt * jac(t+dt, y, SD);
		F = y - y0 - dt * fun(t+dt, y, SD);
		dy = arma::solve(M,F);
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