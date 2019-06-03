
#include <armadillo>
#include <string>
#include "WinchiesController.hpp"

void WinchieController::leer_datosWinchieController(void){

	max_tau = 1000000.0;
	min_tau = 100000.0;

}

void WinchieController::controlWinchies(void){

	double F;

	for(int ii=0; ii<nWinchies; ii=ii+1){

		F = 0.0;

		if (Winchies[ii].LineBCP == 1){
			F = arma::norm(Winchies[ii].LineW->ten_1);
		} else if (Winchies[ii].LineBCP == 2){
			F = arma::norm(Winchies[ii].LineW->ten_N);
		}

		Winchies[ii].tau = F * Winchies[ii].radius;

		if (Winchies[ii].tau>=max_tau) Winchies[ii].tau = max_tau;
		if (Winchies[ii].tau<=min_tau) Winchies[ii].tau = min_tau;
	}
	
}