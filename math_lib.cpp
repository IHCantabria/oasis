/*
Libreria para las matematicas que no se encuentran en librerias.
*/

//LIBRERIAS Y OTROS COMANDOS
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <cmath>
#include <armadillo>
#include "math_lib.h"

spline::spline(arma::mat t, arma::mat x, int pp){

	order = pp;

	int nt1 = min(size(t));
	int nt = max(size(t));
	int nx = max(size(x));

	if(size(x)!=size(t)){
		throw std::invalid_argument( "Inputs must have the same size." );
	}
	if(nt1!=1){
		throw std::invalid_argument( "Inputs must be a sigle row or column." );
	}


	t.reshape(nt,1);
	x.reshape(nt,1);

	T = t;

	if(order == 1){

		COEF = arma::zeros(nt,3);
		arma::mat dxdt = arma::zeros(nt,1);
		dxdt(0) = (x(1)-x(0))/(t(1)-t(0));
		dxdt(nt-1) = (x(nt-1)-x(nt-2))/(t(nt-1)-t(nt-2));
		for(int i=1;i<nt-2;i=i+1){
			dxdt(i) = (x(i+1)-x(i-1))/(t(i+1)-t(i-1));
		}

		arma::mat dxdt2 = arma::zeros(nt,1);
		dxdt2(0) = (dxdt(1)-dxdt(0))/(t(1)-t(0));
		dxdt2(nt-1) = (dxdt(nt-1)-dxdt(nt-2))/(t(nt-1)-t(nt-2));
		for(int i=1;i<nt-2;i=i+1){
			dxdt2(i) = (dxdt(i+1)-dxdt(i-1))/(t(i+1)-t(i-1));
		}

		COEF.col(0) = x;
		COEF.col(1) = dxdt;
		COEF.col(2) = dxdt2;

	} else if (order ==3) {

		arma::mat S0 = arma::zeros(nt,4);
		arma::mat S1 = arma::zeros(nt,4);
		arma::mat S2 = arma::zeros(nt,4);
		arma::mat S3 = arma::zeros(nt,4);
		for(int i=0;i<nt;i=i+1){
			for(int j=0;j<4;j=j+1){
				S0(i,j) = pow(t(i),j);
			}
			for(int j=1;j<4;j=j+1){
				S1(i,j) = j * pow(t(i),j-1);
			}
			for(int j=2;j<4;j=j+1){
				S2(i,j) = (j - 1) * j * pow(t(i),j-2);
			}
			for(int j=3;j<4;j=j+1){
				S3(i,j) = (j - 2) * (j - 1) * j * pow(t(i),j-3);
			}
		}
		int nS = 4*nt-4;
		arma::mat SplineM = arma::zeros(nS,nS);
		SplineM.submat(0,0,0,3) = S0.row(0);
		SplineM.submat(1,0,1,3) = S3.row(0);
		SplineM.submat(nS-2,nS-4,nS-2,nS-1) = S0.row(nt-1);
		SplineM.submat(nS-1,nS-4,nS-1,nS-1) = S3.row(nt-1);
		int ni, mi1, mi2;
		for(int i=0;i<nt-2;i=i+1){
			ni  = 4*i + 2;
			mi1 = 4*i;
			mi2 = 4*i + 3;

			SplineM.submat(ni,mi1,ni,mi2) = S0.row(i+1);

			SplineM.submat(ni+1,mi1,ni+1,mi2) = - S1.row(i+1);
			SplineM.submat(ni+1,mi1+4,ni+1,mi2+4) = S1.row(i+1);

			SplineM.submat(ni+2,mi1,ni+2,mi2) = - S2.row(i+1);
			SplineM.submat(ni+2,mi1+4,ni+2,mi2+4) = S2.row(i+1);

			SplineM.submat(ni+3,mi1+4,ni+3,mi2+4) = S0.row(i+1);
		}
		arma::mat M_x = arma::zeros(nS,nt);
		M_x(0,0) = 1.0;
		M_x(nS-2,nt-1) = 1.0;
		for(int i=0;i<nt-2;i=i+1){
			M_x(2+4*i,i+1) = 1.0;
			M_x(2+4*i+3,i+1) = 1.0;
		}
		arma::mat vec_x = M_x * x;
		COEF = arma::solve(SplineM, vec_x);
	} else {
		throw std::invalid_argument( "Order nor available. Only 1 or 3." );
	}

}


arma::mat spline::spl_eval(double t){

	arma::mat sol = arma::zeros(1,3);

	int ni = max(find(T<=t));

	if(order == 1){

		double dt = t-T(ni);
		sol = COEF.row(ni) + dt * (COEF.row(ni+1)-COEF.row(ni))/(T(ni+1)-T(ni));

	} else if (order ==3){

		arma::mat cc = COEF.rows(ni*4,ni*4+3);

		arma::mat S0 = arma::zeros(1,4);
		arma::mat S1 = arma::zeros(1,4);
		arma::mat S2 = arma::zeros(1,4);
		for(int j=0;j<4;j=j+1){
			S0(0,j) = pow(t,j);
		}
		for(int j=1;j<4;j=j+1){
			S1(0,j) = j * pow(t,j-1);
		}
		for(int j=2;j<4;j=j+1){
			S2(0,j) = (j - 1) * j * pow(t,j-2);
		}

		sol(0,0) = arma::as_scalar(S0 * cc);
		sol(1) = arma::as_scalar(S1 * cc);
		sol(2) = arma::as_scalar(S2 * cc);
	}

	return sol;
}

