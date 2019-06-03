
#include <armadillo>
#include <string>
#include "Hydro.hpp"

// Leer inputs
void Hydro::leer_datosHydro(void){

	int ii;

	hydro = arma::zeros(6*nBodies,6*nBodies);
	arma::cube temp_hydro;
	arma::mat temp_mat;
	temp_hydro.load(arma::hdf5_name("input/flotante.h5","hydro"));
	for(ii=0;ii<nBodies;ii=ii+1){
		temp_mat = temp_hydro( arma::span(ii), arma::span::all, arma::span::all);
		hydro( arma::span(6*ii,6*(ii+1)-1), arma::span(6*ii,6*(ii+1)-1) ) = temp_mat;
	}

	nAngles = 2;
	angles = arma::zeros(nAngles,1);

	nFreqs = 2;
	periods = arma::zeros(nFreqs,1);
	frequencies = arma::zeros(nFreqs,1);

	A.set_size(nBodies,nAngles,nFreqs);

	B.set_size(nBodies,nAngles,nFreqs);

	F_mod.set_size(nBodies,nAngles,nFreqs);
	F_phase.set_size(nBodies,nAngles,nFreqs);
	F_real.set_size(nBodies,nAngles,nFreqs);
	F_imag.set_size(nBodies,nAngles,nFreqs);

	QTFd_mod.set_size(nBodies,nAngles,nFreqs);
	QTFd_phase.set_size(nBodies,nAngles,nFreqs);
	QTFd_real.set_size(nBodies,nAngles,nFreqs);
	QTFd_imag.set_size(nBodies,nAngles,nFreqs);

	QTFs_mod.set_size(nBodies,nAngles,nFreqs);
	QTFs_phase.set_size(nBodies,nAngles,nFreqs);
	QTFs_real.set_size(nBodies,nAngles,nFreqs);
	QTFs_imag.set_size(nBodies,nAngles,nFreqs);

	drag = arma::zeros(nBodies,nAngles);
	current_vel = arma::zeros(2,1);

	kl = arma::zeros(nBodies,6);
	knl = arma::zeros(nBodies,6);

	Ainf = arma::zeros(6*nBodies,6*nBodies);

	inertia = arma::zeros(6*nBodies,6*nBodies);
	arma::cube temp_inertia;
	arma::mat temp_mat2;
	temp_inertia.load(arma::hdf5_name("input/flotante.h5","inertia"));
	for(ii=0;ii<nBodies;ii=ii+1){
		temp_mat2 = temp_inertia( arma::span(ii), arma::span::all, arma::span::all);
		inertia( arma::span(6*ii,6*(ii+1)-1), arma::span(6*ii,6*(ii+1)-1) ) = temp_mat2;
	}

	arma::mat MM = Ainf + inertia;
	invM = arma::solve(MM,arma::eye(6*nBodies,6*nBodies));

	t_relax = 60.0;
	dt_IRF = 0.1;
	nt_IRF = 600;
	vel_hist = arma::zeros(6*nBodies,nt_IRF);
	IRF = arma::zeros(6*nBodies,6*nBodies,nt_IRF);

	t_Fe = 100.0; 
	dt_Fe = 0.1;
	nt_Fe = 1000;
	Fe = arma::zeros(6*nBodies,nt_Fe);

	Hs = 1.0;
	Tp = 10.0;

}

// Calcula la impulse response function
void Hydro::computeIRF(void){
	IRF = arma::zeros(6*nBodies,6*nBodies,nt_IRF);
}

// Calcula el espectro
void Hydro::computeWaveSpectrum(void){

	nComp = 10;
	wave_periods = arma::zeros(nComp,1); 
	wave_frequencies = arma::zeros(nComp,1); 
	wave_heights = arma::zeros(nComp,1); 
	wave_phases = arma::zeros(nComp,1); 
}

// Calcula la serie temporal de fuerzas de excitación
void Hydro::computeFe(void){

	Fe = arma::zeros(6*nBodies,nt_Fe);

}

// Obten las fuerzas hidroestaticas e hidrodinamicas en el tiempo deseado
void Hydro::computeHydroForces(double t){

	HydroForces = arma::zeros(6*nBodies,1);

	arma::mat positions = arma::zeros(6*nBodies,1);
	for(int ii=0;ii<nBodies;ii=ii+1){
		positions(arma::span(6*ii,6*(ii+1)-1),arma::span(0)) = Bodies[ii].pos;
	}	

	HydroForces = HydroForces - hydro*positions;

}