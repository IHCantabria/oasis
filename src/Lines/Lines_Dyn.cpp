
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include "Lines.hpp"
#include "../SEM_math/quadrule.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../SeaFloor/SeaFloor.hpp"

int Line::GetId()
{
	return id;
}

Line::Line(int incId, double incG, double incRhoW, double incFondo)
{
	id = incId;
	g = incG;
	rhoW = incRhoW;
	fondo = incFondo;
}

void Line::ReadPropertiesASCII(FILE *pFilePointer)
{
	// Declare variables
	char buffer_line[1000];
	double dtemp;

	// Ignoro las tres primeras lineas, donde pone "New line"
	for (int ii = 0; ii < 3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFilePointer);
	}

	// Leo todo
	fscanf(pFilePointer, "%d %[^\n]\n", &lineType, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &flag_tension, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &nNodos, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &p, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &L, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &rho0, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &d, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &flag_stiffness, buffer_line);
	if (flag_stiffness == 0)
	{
		fscanf(pFilePointer, "%lf %[^\n]\n", &EA, buffer_line);
		fscanf(pFilePointer, "%lf %[^\n]\n", &beta, buffer_line);
	}
	else if (flag_stiffness == 1)
	{
		fscanf(pFilePointer, "%lf %[^\n]\n", &kernel_coef, buffer_line);
		for (int ii = 0; ii < 3; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the elastic polynomial coefficients for line: " << this->GetId() << "\n";
				throw ValueError(ss.str());
			}
			elastic_coef(ii, 0) = dtemp;
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
		for (int ii = 0; ii < 3; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the viscous loading polynomial coefficients for line: " << this->GetId() << "\n";
				throw ValueError(ss.str());
			}
			visc_loading_coef(ii, 0) = dtemp;
			if (visc_loading_coef(ii, 0) != 0.0){
				flag_visc = 1; // There is viscous component	
			}
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
		for (int ii = 0; ii < 3; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1)
			{	
				std::stringstream ss;
				ss << "An error ocurred when trying to read the viscous unloading polynomial coefficients for line: " << this->GetId() << "\n";
				throw ValueError(ss.str());
			}
			visc_unloading_coef(ii, 0) = dtemp;
			if (visc_unloading_coef(ii, 0) != 0.0){
				flag_visc = 1; // There is viscous component	
			}
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
		// Compute EA aproximation for QS initial condition
		if (elastic_coef(0, 0)>0.0)
		{
			EA = elastic_coef(0, 0);
		}
		else
		{
			double strain_EA = 0.1;
			EA = elastic_coef(0, 0) + 
			     2*elastic_coef(1, 0)*strain_EA +
				 3*elastic_coef(2, 0)*strain_EA*strain_EA;
		}
	}
	else if (flag_stiffness > 1)
	{
		strain_data = arma::zeros(flag_stiffness, 1);
		stress_data = arma::zeros(flag_stiffness, 1);
		for (int ii = 0; ii < this->flag_stiffness; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the strain data for line: " << this->GetId() << "\n";
				throw ValueError(ss.str());
			}
			strain_data(ii, 0) = dtemp;
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
		for (int ii = 0; ii < this->flag_stiffness; ii++)
		{
			if (fscanf(pFilePointer, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error ocurred when trying to read the stress data for line: " << this->GetId() << "\n";
				throw ValueError(ss.str());
			}
			stress_data(ii, 0) = dtemp;
		}
		fscanf(pFilePointer, "%[^\n]\n", buffer_line);
		fscanf(pFilePointer, "%lf %[^\n]\n", &beta, buffer_line);
		// Compute EA aproximation for QS initial condition
		EA = 0;
		int i_EA = 0;
		while (EA<=0)
		{
			EA = (stress_data(i_EA+1,0) - stress_data(i_EA,0)) / 
			     (strain_data(i_EA+1,0) - strain_data(i_EA,0));
		    i_EA++;
			if (i_EA>flag_stiffness-1)
			{
				std::stringstream ss;
				ss << "Strain-Stress curve should have positive slope for line: " << this->GetId() << "\n";
				throw ValueError(ss.str());
			}
		}
	}

	fscanf(pFilePointer, "%lf %[^\n]\n", &CB, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Cmn, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Cdn, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &Cdt, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &GK, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &GC, buffer_line);
	fscanf(pFilePointer, "%d %[^\n]\n", &indexSeaFloor, buffer_line);
	indexSeaFloor--;
	fscanf(pFilePointer, "%d %[^\n]\n", &BCP_N, buffer_line);
	BCP_N -= 1;
	indexBcps[1] = BCP_N;
	fscanf(pFilePointer, "%d %[^\n]\n", &BCP_1, buffer_line);
	BCP_1 -= 1;
	indexBcps[0] = BCP_1;

	if (lineType == 1)
	{
		floor_flag = 1;
	}
	fscanf(pFilePointer, "%d %[^\n]\n", &smoothstep, buffer_line); // UNUSED?!!!!!!!!!
	fscanf(pFilePointer, "%d %[^\n]\n", &frictionModel, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &vth, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &ust, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &usn, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &ud, buffer_line);
	fscanf(pFilePointer, "%lf %[^\n]\n", &deltamax, buffer_line);
	A = arma::datum::pi * d * d * 0.25;
	dL = L / (nNodos - 1);
	dL0 = dL;
	N = p * (nNodos - 1) + 1;
	Kn = Cmn * A * rhoW;

	dampCoef = 2.0 * sqrt(rho0 * GK * d);
	VR = 0.01 * (d * d * GK) / (dampCoef);

	pos = arma::zeros(3 * N);
	vel = arma::zeros(N, 3);
	acc = arma::zeros(N, 3);
	F = arma::zeros(N, 3);
	s = arma::zeros(N);
	xc = arma::zeros(N);
	zc = arma::zeros(N);
	dxcds = arma::zeros(N);
	dzcds = arma::zeros(N);
	Te = arma::zeros(N);
	roots = arma::zeros(p + 1);
	weights = arma::zeros(p + 1);
	FF = arma::zeros(N, 3);
	ff = arma::zeros(N, 3);
	t = arma::zeros(N, 3);
	e_z = arma::zeros(1, 3);
	e_z(0, 2) = 1.0;
	T = arma::zeros(this->N, 1);
	T_visc = arma::zeros(this->N, 1);
	T_elast = arma::zeros(this->N, 1);
	if (flag_stiffness == 1) {
		// Initiallize viscoelasticity buffer vectors
		// TODO: review if the size is well chosen with Y(t0)=1e-3 -> nbuffer = int(t0/dt)
		num_buffer = std::min(std::max(static_cast<int>(ceil(- log(0.001)/(kernel_coef*dt))), 1000), 10000);
		//min 1000 max 10000
		time_vector = arma::zeros(num_buffer, 1);
		strain_vector = arma::zeros(N, num_buffer);
		strain_rate_vector = arma::zeros(N, num_buffer);
		visc_resp_tmp_vector = arma::zeros(N, num_buffer);
		visc_resp_0 = arma::zeros(N, 1); 
		zc_times = arma::zeros(N, num_buffer);
		zc_num = arma::zeros<arma::Mat<int>>(N, 1);
		zc_visc_resp_load = arma::zeros(N, num_buffer);
		zc_visc_resp_unload = arma::zeros(N, num_buffer);
	}
	projectionDirection_1 = arma::zeros(1, 3);
	projectionDirection_1(0, 2) = 1.0;
	projectionDirection_N = arma::zeros(1, 3);
	projectionDirection_N(0, 2) = 1.0;
	posFriccion = arma::zeros(N, 3);
	isSlip = arma::zeros(N, 1);

	if (frictionModel == 1)
	{
		// voy a resolver un solo sistema que dara los coeficientes por si hay friccion de velocidad al principio del problemna.
		// es el step fuerte
		// Los coeficientes van a ser siempre esos asi que solo hace falta hacerlo una vez y asi se obtiene el polinomio
		arma::mat tmpA = arma::zeros(4, 4);
		arma::mat tmpB = arma::zeros(4, 1);
		double x1 = 0;
		double x2 = 1e-04;
		tmpA = {{pow(x1, 3), pow(x1, 2), 1 * x1, 1},
				{3 * pow(x1, 2), 2 * (x1), 1, 0},
				{pow(x2, 3), pow(x2, 2), 1 * x2, 1},
				{3 * pow(x2, 2), 2 * x2, 1, 0}};
		tmpB(0, 0) = 1e-07;
		tmpB(1, 0) = 0.0;
		tmpB(2, 0) = x2;
		tmpB(3, 0) = 1;
		a_1 = arma::solve(tmpA, tmpB);
	}
	double *roots_temp = new double[p + 1];
	double *weights_temp = new double[p + 1];

	lobatto_set(p + 1, roots_temp, weights_temp);

	for (int ii = 0; ii < p + 1; ii++)
	{
		roots(ii) = roots_temp[ii];
		weights(ii) = weights_temp[ii];
	}

	int kk;
	for (int ii = 0; ii < N; ii++)
	{
		kk = ii % p;
		s(ii, 0) = dL * ((ii - kk) / p + (roots(kk) + 1.0) * 0.5);
	}
}

void Line::SEM_getBaseFunctions(void)
{
	pos.reshape(3, N);
	pos = pos.t();

	this->SEM_coefficients();

	arma::mat D_local = SEM_get_D_local();
	D = arma::zeros(N, N);
	for (int ii = 0; ii < nNodos - 1; ii = ii + 1)
	{
		D.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) = D.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) + D_local;
	}
	for (int ii = 1; ii < nNodos - 1; ii = ii + 1)
	{
		D.row(ii * p) = D.row(ii * p) * 0.5;
	}

	int nIntegrate = p + 2;
	double *roots_temp = new double[nIntegrate + 1];
	double *weights_temp = new double[nIntegrate + 1];
	lobatto_set(nIntegrate + 1, roots_temp, weights_temp);

	arma::mat MassMatrix_local = arma::zeros(p + 1, p + 1);
	arma::mat StiffMatrix_local = arma::zeros(p + 1, p + 1);
	arma::mat MSMatrix_local = arma::zeros(p + 1, p + 1);

	for (int i = 0; i < p + 1; i = i + 1)
	{
		for (int j = 0; j < p + 1; j = j + 1)
		{
			for (int k = 0; k < nIntegrate + 1; k = k + 1)
			{
				MassMatrix_local(i, j) = MassMatrix_local(i, j) + weights_temp[k] * this->SEM_poly(roots_temp[k], i) * this->SEM_poly(roots_temp[k], j);
				StiffMatrix_local(i, j) = StiffMatrix_local(i, j) + weights_temp[k] * this->SEM_poly_first_derivative(roots_temp[k], i) * this->SEM_poly_first_derivative(roots_temp[k], j);
				MSMatrix_local(i, j) = MSMatrix_local(i, j) + weights_temp[k] * this->SEM_poly_first_derivative(roots_temp[k], i) * this->SEM_poly(roots_temp[k], j);
			}
		}
	}

	MassMatrix = arma::zeros(N, N);
	StiffMatrix = arma::zeros(N, N);
	MSMatrix = arma::zeros(N, N);

	for (int ii = 0; ii < nNodos - 1; ii = ii + 1)
	{
		MassMatrix.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) = MassMatrix.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) + MassMatrix_local;
		StiffMatrix.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) = StiffMatrix.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) + StiffMatrix_local;
		MSMatrix.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) = MSMatrix.submat(ii * p, ii * p, (ii + 1) * p, (ii + 1) * p) + MSMatrix_local;
	}

	MM = 0.5 * (rho0 + Kn) * MassMatrix;
	inv_MM = arma::solve(MM, arma::eye(N, N));
	MM.row(0) = arma::zeros(1, N);
	MM(0, 0) = 1.0;
	inv_MM_1 = arma::solve(MM, arma::eye(N, N));
	MM.row(N - 1) = arma::zeros(1, N);
	MM(N - 1, N - 1) = 1.0;
	inv_MM_1N = arma::solve(MM, arma::eye(N, N));

	MM = 0.5 * (rho0 + Kn) * MassMatrix;
	MM.row(N - 1) = arma::zeros(1, N);
	MM(N - 1, N - 1) = 1.0;
	inv_MM_N = arma::solve(MM, arma::eye(N, N));

	MM = 0.5 * (rho0 + Kn) * MassMatrix;

	D_sp = arma::sp_mat(D);
	MassMatrix_sp = arma::sp_mat(MassMatrix);
	StiffMatrix_sp = arma::sp_mat(StiffMatrix);
	MSMatrix_sp = arma::sp_mat(MSMatrix);
	MM_sp = arma::sp_mat(MM);
}

void Line::SEM_coefficients(void)
{
	C = arma::zeros(p + 1, p + 1);
	arma::mat M = arma::zeros(p + 1, p + 1);
	for (int ii = 0; ii < p + 1; ii = ii + 1)
	{
		M.insert_cols(ii, arma::pow(roots, ii));
	}
	arma::mat I = arma::eye(p + 1, p + 1);
	C = arma::solve(M, I);
	C = C.t();
}

double Line::SEM_poly(double x, int i)
{
	double y = 0.0;
	arma::mat cc = C.row(i);
	for (int ii = 0; ii < p + 1; ii = ii + 1)
	{
		y = y + cc(ii) * pow(x, ii);
	}
	return y;
}

double Line::SEM_poly_first_derivative(double x, int i)
{
	double y = 0.0;
	arma::mat cc = C.row(i);
	for (int ii = 1; ii < p + 1; ii = ii + 1)
	{
		y = y + ii * cc(ii) * pow(x, ii - 1);
	}
	return y;
}

arma::mat Line::SEM_get_D_local(void)
{
	arma::mat D = arma::zeros(p + 1, p + 1);
	for (int ii = 0; ii < p + 1; ii = ii + 1)
	{
		for (int jj = 0; jj < p + 1; jj = jj + 1)
		{
			D(ii, jj) = this->SEM_poly_first_derivative(roots(ii), jj);
		}
	}
	return D;
}

double Line::kernel(double time)
{
	/**
	 * Proposed Kernel for the memory capability of the viscoelastic material
	 */
	double kernel = exp(-kernel_coef*time); // Usamos exp de <cmath>
	return kernel;
}

double Line::kernel_rate(double time)
{
	double kernel_rate = -kernel_coef * kernel(time);
	return kernel_rate;
}

void Line::smooth_tension(void)
{
	T = 0.5 * (T + arma::abs(T));
		double T0 = 0.1; // size of the smoothing region in newtons
		double a2 = 2.0 / T0, a3 = -1.0 / (T0 * T0);
		arma::mat temp_T = a3 * arma::pow(T, 3) + a2 * arma::pow(T, 3);
		arma::umat ind = arma::find(((T > 0.0) && (T < T0)));
		T.elem(ind) = temp_T.elem(ind);
}

void Line::compute_tension(double time)
{	
	// Compute strain
	arma::mat temp_strain;
	temp_strain = norm_drds - dL / dL0;

	if (flag_stiffness == 0)
	{
		T = EA * temp_strain + beta * dedt;
		
		if ((flag_tension == 2))
		{
			smooth_tension();
		}
	}
	else if (flag_stiffness == 1)
	{	
		T = arma::zeros(N, 1);
		for (int k = 0; k < N; k++)
		{	
			// Extract current node strain and strain rate
			double strain_k = temp_strain(k,0);
			// Compute the elastic tension
			double elastic_tension_k = 0.0;
			for (int i = 0; i < 3; i++){
				elastic_tension_k = elastic_tension_k + elastic_coef(i,0)*pow(strain_k, i+1);
			}
			T_elast(k,0) = elastic_tension_k;
		}
		
		// Check if viscous tension should be computed/updated
		if ((flag_visc == 1) && (time>last_time+dt) && (time>0.0))
		{ 	
			// Update last time
			// std::cout << "  --> Line::compute_tension: Update last time " << std::endl;
			last_time = time;
			num_time_steps += 1;
			// shift time values to the left
			// std::cout << "  --> Line::compute_tension: shift time values to the left time_vector " << std::endl;
			for (int i = 0; i < num_buffer - 1; i++) {
				time_vector[i] = time_vector[i + 1];
			}
			// assign new value to last position
			time_vector[num_buffer - 1] = time;
			// shift strain columns to the left
			// std::cout << "  --> Line::compute_tension: shift time values to the left strain_vector " << std::endl;
			for (int i = 0; i < num_buffer - 1; i++) {
				strain_vector.col(i) = strain_vector.col(i + 1);
			}
			// assign new value to last position
			strain_vector.col(num_buffer - 1) = temp_strain;
			// shift strain rate columns to the left
			// std::cout << "  --> Line::compute_tension: shift time values to the left strain_rate_vector " << std::endl;
			for (int i = 0; i < strain_rate_vector.n_cols - 1; i++) {
				strain_rate_vector.col(i) = strain_rate_vector.col(i + 1);
			}
			// assign new value to last position
			strain_rate_vector.col(strain_rate_vector.n_cols - 1) = dedt;

			// Extract from time vector buffer the assigned values
			// std::cout << "  --> Line::compute_tension: Extract from time vector buffer the assigned values " << std::endl;
			arma::mat tau;
			if (num_time_steps<num_buffer){
				tau = time_vector.rows(num_buffer - num_time_steps, num_buffer - 1);
			} 
			else
			{
				tau = time_vector;
			}
			// Recompute kernel for integral [dY/dtau(t-tau)]
			// std::cout << "  --> Line::compute_tension: Recompute kernel for integral [dY/dtau(t-tau)] " << std::endl;
			arma::mat kernel_rate_tmtau = arma::zeros(size(tau));
			int num_tau = std::min(num_buffer,num_time_steps);
			for (int i = 0; i < num_tau; i++)
			{	
				kernel_rate_tmtau(i,0) = kernel_rate(time-arma::as_scalar(tau(i,0)));
			}

			// Loop over all nodes
			// std::cout << "  --> Line::compute_tension: Loop over all nodes " << std::endl;
			for (int k = 0; k < N; k++)
			{
				// Extract current node strain and strain rate
				// std::cout << "  --> Line::compute_tension: Extract current node strain and strain rate " << std::endl;
				double strain_k = temp_strain(k,0);
				double strain_rate_k = dedt(k,0);

				// Check if current state is loading or unloading and choose polynomial coefficients
				// TODO: filtering may be necessary
				// std::cout << "  --> Line::compute_tension: Check if current state is loading or unloading and choose polynomial coefficients " << std::endl;
				arma::mat strain_poly_coef = arma::zeros(3, 1);
				if (strain_rate_k < -tol_zero){
					strain_poly_coef = visc_unloading_coef;
				}
				else {
					strain_poly_coef = visc_loading_coef;
				}
				
				// Compute instant viscous response
				// std::cout << "  --> Line::compute_tension: Compute instant viscous response " << std::endl;
				double visc_resp_tmp_k = 0.0; // gv(strain(t))
				for (int i = 0; i < 3; i++){
					visc_resp_tmp_k += strain_poly_coef(i,0)*pow(strain_k,i+1);
				}

				// shift strain columns to the left
				// std::cout << "  --> Line::compute_tension: shift strain columns to the left visc_resp_tmp_vector " << std::endl;
				for (int i = 0; i < num_buffer - 1; i++) {
					visc_resp_tmp_vector(k,i) = visc_resp_tmp_vector(k,i+1);
				}
				// assign new value to last position
				visc_resp_tmp_vector(k, num_buffer - 1) = visc_resp_tmp_k;
				
				// Extract from viscous response vector buffer the assigned values
				// std::cout << "  --> Line::compute_tension: Extract from viscous response vector buffer the assigned values " << std::endl;
				arma::mat visc_resp;
				if (num_time_steps<num_buffer){
					visc_resp = visc_resp_tmp_vector(k,arma::span(num_buffer - num_time_steps, num_buffer - 1));
				} 
				else
				{
					visc_resp = visc_resp_tmp_vector.row(k);
				}

				// Integral
				// std::cout << "  --> Line::compute_tension: Compute Integral " << std::endl;
				arma::mat integrand = kernel_rate_tmtau.t() % visc_resp;
				double visc_tension_integral_k = trapzi(tau.t(), integrand);
				// if (time>100){
				// std::cout << "  --> Line::compute_tension: tau.t() = "<<tau.t()<<" " << std::endl;
				// std::cout << "  --> Line::compute_tension: integrand = "<<integrand<<" " << std::endl;
				// throw std::exception();
				// }
				// Initial tension term from the integral
				// std::cout << "  --> Line::compute_tension: Compute Initial tension term from the integral " << std::endl;
				double final_tension_integral_k = visc_resp_tmp_k;
				// Final tension term from the integral
				// std::cout << "  --> Line::compute_tension: Compute Final tension term from the integral " << std::endl;
				double initial_tension_integral_k = kernel(time)*visc_resp_0(k,0);


				// Check for new strain rate zero crossing events
				// std::cout << "  --> Line::compute_tension: Check for new strain rate zero crossing events " << std::endl;
				// TODO: Filtering may be necesary
				double time_old = time_vector[num_buffer-2];
				double strain_k_old = strain_vector(k,num_buffer-2);
				double strain_rate_k_old = strain_rate_vector(k,num_buffer-2);
				double strain_k_3_old = strain_vector(k,num_buffer-3);
				double strain_rate_k_3_old = strain_rate_vector(k,num_buffer-3);
				if ((arma::sign(strain_rate_k) != arma::sign(strain_rate_k_old)) &&
					(std::abs(strain_rate_k)>tol_zero) && 
					(std::abs(strain_rate_k_old)>tol_zero))
				{
					// Update the number of zero crossings for the current node
					// std::cout << "  --> Line::compute_tension: Update the number of zero crossings for the current node " << std::endl;
					zc_num(k,0) += 1;

					// Compute the zero crossing time
					// std::cout << "  --> Line::compute_tension: Compute the zero crossing time " << std::endl;
					double dt_temp = std::abs(time-time_old);
					double zc_time_tmp = time_old 
										- strain_rate_k_old * dt_temp
										/ (strain_rate_k - strain_rate_k_old);
					
					// shift strain columns to the left
					// std::cout << "  --> Line::compute_tension: shift strain columns to the left zc_times " << std::endl;
					for (int i = 0; i < num_buffer - 1; i++) {
						zc_times(k,i) = zc_times(k,i+1);
					}
					// assign new value to last position
					zc_times(k, num_buffer - 1) = zc_time_tmp;
					
					// Estimate strain at strain rate zero crossing time
					// std::cout << "  --> Line::compute_tension: Estimate strain at strain rate zero crossing time " << std::endl;

					// Cubic
					// double zc_strain_tmp = cubic_interp(time_old,
					// 									time, 
					// 									strain_k_old, 
					// 									strain_k, 
					// 									strain_rate_k_old, 
					// 									strain_rate_k, 
					// 									zc_time_tmp);

					// Linear
					double zc_strain_tmp = strain_k_old + (strain_k-strain_k_old)/(time-time_old)*(zc_time_tmp-time_old);

					// Compute viscous loading and unloading response at new zero crossing
					// std::cout << "  --> Line::compute_tension: Compute viscous loading and unloading response at new zero crossing " << std::endl;
					double visc_resp_load_tmp_k = 0.0;
					double visc_resp_unload_tmp_k = 0.0;
					for (int j = 0; j < 3; j++) {
						visc_resp_load_tmp_k += (pow(zc_strain_tmp, j + 1) * visc_loading_coef(j, 0));
						visc_resp_unload_tmp_k += (pow(zc_strain_tmp, j + 1) * visc_unloading_coef(j, 0));
					}

					// shift strain columns to the left
					// std::cout << "  --> Line::compute_tension: shift strain columns to the left zc_visc_resp_load" << std::endl;
					for (int i = 0; i < num_buffer - 1; i++) {
						zc_visc_resp_load(k,i) = zc_visc_resp_load(k,i+1);
						zc_visc_resp_unload(k,i) = zc_visc_resp_unload(k,i+1);
					}
					// assign new value to last position
					zc_visc_resp_load(k, num_buffer - 1) = visc_resp_load_tmp_k;
					zc_visc_resp_unload(k, num_buffer - 1) = visc_resp_unload_tmp_k;
				}

				double sum_stress_integral = 0.0;
				if (zc_num(k,0)>0)
				
				{
					// Extract from zero crossing vector buffers the assigned values
					// std::cout << "  --> Line::compute_tension: Extract from zero crossing vector buffers the assigned values" << std::endl;
					arma::mat zc_times_k, zc_visc_resp_load_k, zc_visc_resp_unload_k;
					if (zc_num(k,0)<num_buffer){
						zc_times_k = zc_times(k,arma::span(num_buffer - zc_num(k,0), num_buffer - 1));
						zc_visc_resp_load_k = zc_visc_resp_load(k,arma::span(num_buffer - zc_num(k,0), num_buffer - 1));
						zc_visc_resp_unload_k = zc_visc_resp_unload(k,arma::span(num_buffer - zc_num(k,0), num_buffer - 1));
					} 
					else
					{
						zc_times_k = zc_times.row(k);
						zc_visc_resp_load_k = zc_visc_resp_load.row(k);
						zc_visc_resp_unload_k = zc_visc_resp_unload.row(k);
					}

					// Sum of initial/final terms from the integral between zero crossing events
					// std::cout << "  --> Line::compute_tension: Sum of initial/final terms from the integral between zero crossing events" << std::endl;
					int num_zc_tmp_k = std::min(num_buffer,zc_num(k,0));
					for (int i = 0; i < num_zc_tmp_k; i++) {
						sum_stress_integral += kernel(time - zc_times_k[i])
											* pow(-1,i)
											* (zc_visc_resp_load_k[i] - zc_visc_resp_unload_k[i]);							
					}
				}

				// Save the computed tension in the class atribute vector
				// std::cout << "  --> Line::compute_tension: Save the computed tension in the class atribute vector" << std::endl;
				T_visc(k,0) = visc_tension_integral_k
							+ final_tension_integral_k
							- initial_tension_integral_k
							+ sum_stress_integral;

				if (k == 0){
					fprintf(pfile_debug, "%f    ", time);
					fprintf(pfile_debug, "%f    ", strain_k);
					fprintf(pfile_debug, "%f    ", strain_rate_k);
					fprintf(pfile_debug, "%f    ", T_elast(k,0));
					fprintf(pfile_debug, "%f    ", visc_tension_integral_k);
					fprintf(pfile_debug, "%f    ", final_tension_integral_k);
					fprintf(pfile_debug, "%f    ", initial_tension_integral_k);
					fprintf(pfile_debug, "%f    ", sum_stress_integral);
					fprintf(pfile_debug, "\n");
				}
			}
		}

		// Add elastic and viscous tension
		T = T_elast + T_visc;
		// T = T_elast;

		if ((flag_tension == 2))
		{
			smooth_tension();
		}
	}
	else if (flag_stiffness > 1)
	{
		// TODO: investigate how to use qubic interpolation
		T = interp1(strain_data, stress_data, temp_strain) + beta * dedt;
	}
	
}

void Line::SEM_compute_derivarives(void)
{
	// Compute the variables needed to get strain and strain rate
	drds = (D * pos) * (2.0 / dL0);
	drdsdt = (D * vel) * (2.0 / dL0);
	norm_drds = sqrt(pow(drds.col(0), 2) + pow(drds.col(1), 2) + pow(drds.col(2), 2));
	dedt = drds.col(0) % drdsdt.col(0) + drds.col(1) % drdsdt.col(1) + drds.col(2) % drdsdt.col(2);
}

void Line::SEM_computeF(double time)
{
	SEM_compute_derivarives();
	// Compute tension
	compute_tension(time);	

	if (floor_flag == 1)
	{
		arma::field<arma::mat> temp = pLineSeaFloor->projectPoints(pos);
		// ahora toman valor
		projectedPoints = temp(0, 0);
		projectionDirection = temp(0, 2);
		zCoordinates = temp(0, 1);
	}

	for (int k = 0; k < N; k++)
	{
		t.row(k) = drds.row(k) / norm_drds(k);
		FF.row(k) = T(k) * t.row(k);

		fg = (rho0 - rhoW * A) * g / norm_drds(k);
		ff.row(k) = -fg * e_z;

		v = vel.row(k);
		vt = (v * t.row(k).t()) * t.row(k);
		vn = v - vt;

		ff.row(k) = ff.row(k) - 0.5 * Cdt * d * rhoW * arma::norm(vt, 2) * vt;
		ff.row(k) = ff.row(k) - 0.5 * Cdn * d * rhoW * arma::norm(vn, 2) * vn;

		if (floor_flag == 1)
		{

			// GROUND NORMAL FORCES -- SMOOTHED PALM'S MODEL
			double ultimaCoordVel = arma::as_scalar(v * arma::strans(projectionDirection.row(k)));
			double zCoordinate = -zCoordinates(k); // porque al ir las normales hacia arriba es el criterio contrario

			double paramNormal, parammuelle1, parammuelle2, paramVel;
			paramNormal = step(zCoordinate, -d / 2, 0, 0, 1);
			paramVel = step(ultimaCoordVel, -VR, 1, 0, 0);
			parammuelle1 = step(zCoordinate, 0, 0, d / 2, 1);
			parammuelle2 = GK * d * (zCoordinate);

			ff.row(k) = ff.row(k) + (fg * paramNormal + parammuelle1 * parammuelle2 - GC * paramNormal * dampCoef * paramVel * ultimaCoordVel) * projectionDirection.row(k) / arma::norm(projectionDirection.row(k), 2);

			if (k == 0)
			{
				paramNormal_1 = paramNormal;
				parammuelle1_1 = parammuelle1;
				parammuelle2_1 = parammuelle2;
				paramVel_1 = paramVel;
				ultimaCoordVel_1 = ultimaCoordVel;
				projectionDirection_1 = projectionDirection.row(k);
			}
			if (k == N - 1)
			{
				paramNormal_N = paramNormal;
				parammuelle1_N = parammuelle1;
				parammuelle2_N = parammuelle2;
				paramVel_N = paramVel;
				ultimaCoordVel_N = ultimaCoordVel;
				projectionDirection_N = projectionDirection.row(k);
			}

			// FRICTION FORCES
			// STICK-SLIP MODEL
			if (frictionModel == 2 || frictionModel == 1)
			{

				fn = fg * step(zCoordinate, -d / 2, 0, 0, 1);

				if (fn > 0)
				{ // en otro caso se considera que no hay friccion

					// parametros a utilizar
					arma::mat vpi = arma::zeros(1, 3); // velocidad en el plano del suelo
					vpi = v - arma::as_scalar(v * arma::strans(projectionDirection.row(k))) * projectionDirection.row(k) / arma::norm(projectionDirection.row(k), 2);
					double normVel = arma::norm(vpi, 2);

					// isSlip vale 0 si en la iter.anterior es slip
					// isSlip vale 1 si en la iter.anterior es stick
					double betaStick, usstep, udstep;

					if (frictionModel == 1)
					{
						if (isSlip(k, 0) == 1)
						{
							if (normVel > vth)
							{
								// en este caso cambia a slip
								isSlip(k, 0) = 0;
								// no guarda posicion, ya que solo la guarda si pasa de slip a stick
							}
						}
						if (isSlip(k, 0) == 0)
						{
							if (normVel <= vth)
							{
								// en este caso cambia a stick
								isSlip(k, 0) = 1;
								// se deben guardar las posiciones cuando ocurre esto
								posFriccion.row(k) = projectedPoints.row(k);
							}
						}

						arma::mat deformationpi = projectedPoints.row(k) - posFriccion.row(k); // deformacion en el plano del suelo
						double delta = arma::norm(deformationpi, 2);

						// parametros distintos si es stick o slip
						double usstep;
						if (isSlip(k, 0) == 0)
						{
							betaStick = 1;
							usstep = 0;
							udstep = ud;
						}
						else
						{
							// betaStick=step(normVel, -vth, -1.0, vth, 1.0);
							betaStick = step(normVel, 0.0, 0.0, vth, 1.0);
							usstep = step(delta, 0.0, 0.0, deltamax, ust);
							// usstep=step(delta, -deltamax, -ust, deltamax, ust);
							// udstep=step(normVel, -vth, -ud, vth, ud);
							udstep = step(normVel, 0.0, 0.0, vth, ud);
						}
						// hay que calcular las fuerzas
						double tol_max;
						tol_max = vth * 0.01;
						double velocidadSuave = std::max(normVel, tol_max);
						arma::mat ffslid = -vpi * udstep * fn * step(normVel, tol_max, 0, 2 * tol_max, 1) / velocidadSuave;
						ff.row(k) = ff.row(k) + ffslid;
						// std::cout << "ffslid del nodo " << k << "= " << ffslid << std::endl;
						if ((isSlip(k, 0) == 1))
						{
							tol_max = deltamax * 0.01;
							double deltaSuave = std::max(delta, tol_max);
							arma::mat ffstick = -deformationpi * (1 - betaStick) * usstep * fn * step(delta, tol_max, 0, 2 * tol_max, 1) / deltaSuave;
							ff.row(k) = ff.row(k) + ffstick;
							// std::cout << "delta del nodo " << k << "= " << delta << std::endl;
							// std::cout << "ffstick del nodo " << k << "= " << ffstick << std::endl;
						}
					}
					if (frictionModel == 2)
					{
						arma::mat ffsticktan, ffslidtan, ffslidnorm, ffsticknorm;
						arma::mat tpi = arma::zeros(1, 3);
						tpi = t.row(k) - arma::as_scalar(t.row(k) * arma::strans(projectionDirection.row(k))) * projectionDirection.row(k) / arma::norm(projectionDirection.row(k));
						tpi = tpi / arma::norm(tpi, 2);

						// arma::mat velocityTan = vt-  arma::as_scalar(vt* arma::strans(projectionDirection.row(k)))*projectionDirection.row(k)/ arma::norm(projectionDirection.row(k));
						// arma::mat velocityNorm = vn-  arma::as_scalar(vn* arma::strans(projectionDirection.row(k)))*projectionDirection.row(k)/ arma::norm(projectionDirection.row(k));
						arma::mat velocityTan = arma::as_scalar(vpi * arma::strans(tpi)) * tpi;
						arma::mat velocityNorm = vpi - velocityTan;

						double normVelTan = arma::norm(velocityTan, 2);
						double normVelNorm = arma::norm(velocityNorm, 2);

						if (isSlip(k, 0) == 1)
						{
							if ((normVelTan > vth) || normVelNorm > vth)
							{
								// en este caso cambia a slip
								isSlip(k, 0) = 0;
								// no guarda posicion, ya que solo la guarda si pasa de slip a stick
							}
						}
						if (isSlip(k, 0) == 0)
						{
							if (normVelTan <= vth && normVelNorm <= vth)
							{
								// en este caso cambia a stick
								isSlip(k, 0) = 1;
								// se deben guardar las posiciones cuando ocurre esto
								posFriccion.row(k) = projectedPoints.row(k);
							}
						}

						arma::mat deformationpi = projectedPoints.row(k) - posFriccion.row(k); // deformacion en el plano del suelo
						double delta = arma::norm(deformationpi, 2);
						// arma::mat deformationTan = arma::as_scalar(tpi * arma::strans(deformationpi)) *tpi /arma::norm(tpi, 2);
						arma::mat deformationTan = arma::as_scalar(tpi * arma::strans(deformationpi)) * tpi;
						arma::mat deformationNor = deformationpi - deformationTan;
						double deltaTan = arma::norm(deformationTan, 2);
						double deltaNorm = arma::norm(deformationNor, 2);

						double uststep, usnstep, udstep;
						// parametros distintos si es stick o slip
						if (isSlip(k, 0) == 0)
						{
							betaStick = 1;
							betaStick = 1;
							uststep = 0;
							usnstep = 0;
							udstep = ud;
						}
						else
						{
							betaStick = step(normVel, 0.0, 0.0, vth, 1.0);
							uststep = step(deltaTan, 0.0, 0.0, deltamax, ust);
							usnstep = step(deltaNorm, 0.0, 0.0, deltamax, usn);
							udstep = step(normVel, 0.0, 0.0, vth, ud);
						}
						double tol_max;
						tol_max = vth * 0.01;
						double velocidadSuaveTan = std::max(normVelTan, tol_max);
						double velocidadSuaveNor = std::max(normVelNorm, tol_max);
						// hay que calcular las fuerzas
						ffslidtan = -velocityTan * udstep * fn * step(normVelTan, tol_max, 0, 2 * tol_max, 1) / velocidadSuaveTan;
						ffslidnorm = -velocityNorm * udstep * fn * step(normVelNorm, tol_max, 0, 2 * tol_max, 1) / velocidadSuaveNor;

						ff.row(k) = ff.row(k) + ffslidtan + ffslidnorm;

						if ((isSlip(k, 0) == 1))
						{
							tol_max = deltamax * 0.01;
							double deltaSuaveTan = std::max(deltaTan, tol_max);
							double deltaSuaveNor = std::max(deltaNorm, tol_max);
							ffsticktan = -deformationTan * (1 - betaStick) * uststep * fn * step(deltaTan, tol_max, 0, 2 * tol_max, 1) / deltaSuaveTan;
							ffsticknorm = -deformationNor * (1 - betaStick) * usnstep * fn * step(deltaNorm, tol_max, 0, 2 * tol_max, 1) / deltaSuaveNor;
							ff.row(k) = ff.row(k) + ffsticktan + ffsticknorm;
						}
					}
				}
			}
		}
	}

	// F = 0.5 * dL * (MassMatrix_sp * ff) - (MSMatrix_sp * FF);
	F = 0.5 * dL * (MassMatrix * ff) - (MSMatrix * FF);
	// F.row(0) = F.row(0) + FF.row(0);
	// F.row(N-1) = F.row(N-1) - FF.row(N-1);

	if (F.has_nan())
	{
		std::cout << std::endl
				  << "ERROR: NaN Detected on line with id = " << id << std::endl;
		std::cout << std::endl
				  << "  F = " << F << std::endl;
		std::cout << std::endl
				  << "  ff = " << ff << std::endl;
		std::cout << std::endl
				  << "  FF = " << FF << std::endl;
		std::cout << std::endl
				  << "  pos = " << pos << std::endl;
		std::cout << std::endl
				  << "  vel = " << vel << std::endl;

		throw std::exception();
	}

	ten_1 = FF.row(0).t();
	ten_N = FF.row(N - 1).t();
	F_1 = F.row(0).t();
	F_N = F.row(N - 1).t();
	// F_1 = d*ff.row(0).t() + ten_1;
	// F_N = d*ff.row(N-1).t() - ten_N;

	pLineBcps[0]->forceBcp.rows(0, 2) = pLineBcps[0]->forceBcp.rows(0, 2) + F_1;
	pLineBcps[1]->forceBcp.rows(0, 2) = pLineBcps[1]->forceBcp.rows(0, 2) + F_N;
	pLineBcps[0]->temp = pLineBcps[0]->forceBcp;
	pLineBcps[1]->temp = pLineBcps[1]->forceBcp;
}

void Line::initiallize_strain_memory(void)
{
	SEM_compute_derivarives();
	std::cout << "  --> Working on initiallization of strain memory, derivatives computed succesfully ..." << std::endl;
	arma::mat temp_strain = norm_drds - dL / dL0;
	std::cout << "  --> Working on initiallization of strain memory,norm_drds = "<<norm_drds<<" dL / dL0 ="<<dL / dL0<<" and  "<<temp_strain<<" = temp_strain ..." << std::endl;
	for (int i = 0; i < strain_vector.n_cols; i++) {
		strain_vector.col(i) = temp_strain;
	}
	for (int i = 0; i < 3; i++){
		visc_resp_0 = visc_resp_0 + visc_loading_coef(i,0)*pow(temp_strain,i+1);
	}
	
}

void Line::print_out(void)
{

	std::cout << "Para la linea " << this->id << " , se ha leido:" << std::endl
			  << std::endl;
	std::cout << "nNodos   " << this->nNodos << std::endl;
	std::cout << "p        " << this->p << std::endl;
	std::cout << "L        " << this->L << std::endl;
	std::cout << "rho0     " << this->rho0 << std::endl;
	std::cout << "d        " << this->d << std::endl;
	std::cout << "EA       " << this->EA << std::endl;
	std::cout << "beta     " << this->beta << std::endl;
	std::cout << "CB       " << this->CB << std::endl;
	std::cout << "Cmn      " << this->Cmn << std::endl;
	std::cout << "Cdn      " << this->Cdn << std::endl;
	std::cout << "Cdt      " << this->Cdt << std::endl;
	std::cout << "GK       " << this->GK << std::endl;
	std::cout << "GC       " << this->GC << std::endl;
	std::cout << "pos_N  " << this->pos_N(0, 0) << " " << this->pos_N(1, 0) << " " << this->pos_N(2, 0) << std::endl;
	std::cout << "pos_1  " << this->pos_1(0, 0) << " " << this->pos_1(1, 0) << " " << this->pos_1(2, 0) << std::endl;
	std::cout << "smoothstep      " << this->smoothstep << std::endl;
	std::cout << "frictionModel      " << this->frictionModel << std::endl;
	std::cout << "vth      " << this->vth << std::endl;
	std::cout << "us    tangential  " << this->ust << std::endl;
	std::cout << "us    normal  " << this->usn << std::endl;
	std::cout << "ud      " << this->ud << std::endl;
	std::cout << "deltamax      " << this->deltamax << std::endl
			  << std::endl;
}

void Line::initLine(void)
{
	// Load BCP data
	pos_1 = pLineBcps[0]->posG_BCP;
	pos_N = pLineBcps[1]->posG_BCP;

	// Init Line
	if (pos_N(2, 0) < fondo || pos_1(2, 0) < fondo)
		throw 0;

	if (lineType == 1)
	{
		int flag = 0;
		int flagTense = 0;
		xF = sqrt(pow((pos_N(0, 0) - pos_1(0, 0)), 2) + pow((pos_N(1, 0) - pos_1(1, 0)), 2));
		zF = (pos_N(2, 0) - pos_1(2, 0));
		if (xF == 0)
		{
			cosa = 0;
			sina = 0;
		}
		else
		{
			cosa = (pos_N(0, 0) - pos_1(0, 0)) / xF;
			sina = (pos_N(1, 0) - pos_1(1, 0)) / xF;
		}
		double Li = sqrt(pow(xF, 2) + pow(zF, 2));
		if (Li < L)
		{
			if (pos_1(2, 0) == fondo && pos_N(2, 0) == fondo)
				throw 2;
			if (xF < 1e-5)
				throw 3;
			floor_flag = -1;
			this->qs_GetTen();
			this->qs_Solution();
			if (isnan(xc(0, 0)))
				flagTense = 1;
			for (int ii = 0; ii < N; ii = ii + 1)
			{
				pos(3 * ii, 0) = pos_1(0, 0) + cosa * xc(ii, 0);
				pos(3 * ii + 1, 0) = pos_1(1, 0) + sina * xc(ii, 0);
				pos(3 * ii + 2, 0) = pos_1(2, 0) + zc(ii, 0);
				if (pos(3 * ii + 2, 0) <= fondo && ii >= 1)
					flag = 1;
			}
			if ((flag == 1) && pos_1(2, 0) > fondo)
				throw 6;
			if ((flag == 1) && ((pos_1(2, 0) == fondo) || (pos_N(2, 0) == fondo)))
			{
				floor_flag = 1;
				this->qs_GetTen();
				this->qs_Solution();
				if (isnan(xc(0, 0)))
					throw 5;
				for (int ii = 0; ii < N; ii = ii + 1)
				{
					pos(3 * ii, 0) = pos_1(0, 0) + cosa * xc(ii, 0);
					pos(3 * ii + 1, 0) = pos_1(1, 0) + sina * xc(ii, 0);
					pos(3 * ii + 2, 0) = pos_1(2, 0) + zc(ii, 0);
				}
			}
			arma::mat vec1(pos.rows(3, 5) - pos.rows(0, 2));
			double nvec1 = norm(vec1);
			ten_1 = (Te(0) / nvec1) * vec1;
			arma::mat vecN(pos.rows(3 * N - 6, 3 * N - 4) - pos.rows(3 * N - 3, 3 * N - 1));
			double nvecN = norm(vecN);
			ten_N = (Te(N - 1) / nvecN) * vecN;
		}
		if ((Li >= L) || (flagTense == 1))
		{
			double cost = xF / Li;
			double sint = zF / Li;
			for (int ii = 0; ii < N; ii = ii + 1)
			{
				xc(ii, 0) = cost * (Li / L) * s(ii);
				zc(ii, 0) = sint * (Li / L) * s(ii);
			}
			for (int ii = 0; ii < N; ii = ii + 1)
			{
				pos(3 * ii, 0) = pos_1(0, 0) + cosa * xc(ii, 0);
				pos(3 * ii + 1, 0) = pos_1(1, 0) + sina * xc(ii, 0);
				pos(3 * ii + 2, 0) = pos_1(2, 0) + zc(ii, 0);
			}
			arma::mat vec1(pos.rows(3, 5) - pos.rows(0, 2));
			double nvec1 = norm(vec1);
			ten_1 = (EA * (nvec1 - 1.0) / (0.5 * (roots(1) + 1) * dL * nvec1)) * vec1;
			arma::mat vecN(pos.rows(3 * N - 6, 3 * N - 4) - pos.rows(3 * N - 3, 3 * N - 1));
			double nvecN = norm(vecN);
			ten_N = (EA * (nvecN - 1.0) / (0.5 * (roots(1) + 1) * dL * nvecN)) * vecN;
			if (Li == L && pos_1(2, 0) == fondo && pos_N(2, 0) == fondo)
			{
				std::cout << "     WARNING: Mooring line " << id << " is laying on the floor " << std::endl;
			}
			else if (xF < 1e-5)
			{
				std::cout << "     WARNING: Mooring line " << id << " is vertical. " << std::endl;
			}
			else
			{
				std::cout << "     WARNING: Mooring line " << id << " tension is high. " << std::endl;
			}
		}
		floor_flag = 1;
		std::cout << "    Tension at anchor for line: " << id << " ; is: (" << ten_1(0) << " , " << ten_1(1) << " , " << ten_1(2) << " ) N" << std::endl;
		std::cout << "    Tension at fairlead for line: " << id << " ; is: (" << ten_N(0) << " , " << ten_N(1) << " , " << ten_N(2) << " ) N" << std::endl
				  << std::endl;
	}
	else if (lineType == 2)
	{
		floor_flag = -1;
		int flagTense = 0;
		xF = sqrt(pow((pos_N(0, 0) - pos_1(0, 0)), 2) + pow((pos_N(1, 0) - pos_1(1, 0)), 2));
		zF = (pos_N(2, 0) - pos_1(2, 0));
		if (xF == 0)
		{ // Vertical case
			cosa = 0;
			sina = 0;
		}
		else
		{ // other
			cosa = (pos_N(0, 0) - pos_1(0, 0)) / xF;
			sina = (pos_N(1, 0) - pos_1(1, 0)) / xF;
		}
		double Li = sqrt(pow(xF, 2) + pow(zF, 2));
		if (Li < L)
		{
			this->qs_GetTen();
			this->qs_Solution();
			if (isnan(xc(0, 0)))
				flagTense = 1;
			for (int ii = 0; ii < N; ii = ii + 1)
			{
				pos(3 * ii, 0) = pos_1(0, 0) + cosa * xc(ii, 0);
				pos(3 * ii + 1, 0) = pos_1(1, 0) + sina * xc(ii, 0);
				pos(3 * ii + 2, 0) = pos_1(2, 0) + zc(ii, 0);
				if (pos(3 * ii + 2, 0) <= fondo)
					throw 1;
			}
			arma::mat vec1(pos.rows(3, 5) - pos.rows(0, 2));
			double nvec1 = norm(vec1);
			ten_1 = (Te(0) / nvec1) * vec1;
			arma::mat vecN(pos.rows(3 * N - 6, 3 * N - 4) - pos.rows(3 * N - 3, 3 * N - 1));
			double nvecN = norm(vecN);
			ten_N = (Te(N - 1) / nvecN) * vecN;
		}
		if ((Li >= L) || (flagTense == 1))
		{
			double cost = xF / Li;
			double sint = zF / Li;
			for (int ii = 0; ii < N; ii = ii + 1)
			{
				xc(ii, 0) = cost * (Li / L) * s(ii);
				zc(ii, 0) = sint * (Li / L) * s(ii);
			}
			for (int ii = 0; ii < N; ii = ii + 1)
			{
				pos(3 * ii, 0) = pos_1(0, 0) + cosa * xc(ii, 0);
				pos(3 * ii + 1, 0) = pos_1(1, 0) + sina * xc(ii, 0);
				pos(3 * ii + 2, 0) = pos_1(2, 0) + zc(ii, 0);
			}
			arma::mat vec1(pos.rows(3, 5) - pos.rows(0, 2));
			double nvec1 = norm(vec1);
			ten_1 = (EA * (nvec1 - 1.0) / (0.5 * (roots(1) + 1) * dL * nvec1)) * vec1;
			arma::mat vecN(pos.rows(3 * N - 6, 3 * N - 4) - pos.rows(3 * N - 3, 3 * N - 1));
			double nvecN = norm(vecN);
			ten_N = (EA * (nvecN - 1.0) / (0.5 * (roots(1) + 1) * dL * nvecN)) * vecN;
			std::cout << "    WARNING: Towing line " << id << " tension is high. " << std::endl;
		}
		std::cout << "    Tension at anchor for line: " << id << " ; is: (" << ten_1(0) << " , " << ten_1(1) << " , " << ten_1(2) << " ) N" << std::endl;
		std::cout << "    Tension at fairlead for line: " << id << " ; is: (" << ten_N(0) << " , " << ten_N(1) << " , " << ten_N(2) << " ) N" << std::endl
				  << std::endl;
	}
	else if (lineType == 3)
	{
		xF = sqrt(pow((pos_N(0, 0) - pos_1(0, 0)), 2) + pow((pos_N(1, 0) - pos_1(1, 0)), 2));
		zF = (pos_N(2, 0) - pos_1(2, 0));
		double Li = sqrt(pow(xF, 2) + pow(zF, 2));
		if (Li < L)
			throw 4;
		if (xF == 0)
		{
			cosa = 0;
			sina = 0;
		}
		else
		{
			cosa = (pos_N(0, 0) - pos_1(0, 0)) / xF;
			sina = (pos_N(1, 0) - pos_1(1, 0)) / xF;
		}
		double cost = xF / Li;
		double sint = zF / Li;
		for (int ii = 0; ii < N; ii = ii + 1)
		{
			xc(ii, 0) = cost * (Li / L) * s(ii);
			zc(ii, 0) = sint * (Li / L) * s(ii);
			pos(3 * ii, 0) = pos_1(0, 0) + cosa * xc(ii, 0);
			pos(3 * ii + 1, 0) = pos_1(1, 0) + sina * xc(ii, 0);
			pos(3 * ii + 2, 0) = pos_1(2, 0) + zc(ii, 0);
		}
		arma::mat vec1(pos.rows(3, 5) - pos.rows(0, 2));
		double nvec1 = norm(vec1);
		ten_1 = (EA * (nvec1 - 1.0) / (0.5 * (this->roots(1) + 1) * dL * nvec1)) * vec1;
		arma::mat vecN(pos.rows(3 * N - 6, 3 * N - 4) - pos.rows(3 * N - 3, 3 * N - 1));
		double nvecN = norm(vecN);
		ten_N = (EA * (nvecN - 1.0) / (0.5 * (this->roots(1) + 1) * dL * nvecN)) * vecN;
		std::cout << "     Tension at anchor for line: " << id << " ; is: (" << ten_1(0) << " , " << ten_1(1) << " , " << ten_1(2) << " ) N" << std::endl;
		std::cout << "     Tension at fairlead for line: " << id << " ; is: (" << ten_N(0) << " , " << ten_N(1) << " , " << ten_N(2) << " ) N" << std::endl
				  << std::endl;
		floor_flag = -1;
	}
	else
	{
		throw std::invalid_argument("Type of line not available.");
	}
}

void Line::OpenOutputFilesASCII(std::string path)
{
	char buffer1[50], buffer2[50], buffer3[50], buffer4[50], buffer5[50], buffer6[50];

	int nn1 = sprintf(buffer1, "NodePosX_%d.txt", GetId());
	int nn2 = sprintf(buffer2, "NodePosY_%d.txt", GetId());
	int nn3 = sprintf(buffer3, "NodePosZ_%d.txt", GetId());
	int nn4 = sprintf(buffer4, "EndsTen_%d.txt", GetId());
	int nn5 = sprintf(buffer5, "LineTen_%d.txt", GetId());
	int nn6 = sprintf(buffer6, "debug_%d.txt", GetId());

	std::string file_path1 = JoinPath(path, buffer1);
	std::string file_path2 = JoinPath(path, buffer2);
	std::string file_path3 = JoinPath(path, buffer3);
	std::string file_path4 = JoinPath(path, buffer4);
	std::string file_path5 = JoinPath(path, buffer5);
	std::string file_path6 = JoinPath(path, buffer6);

	pfile_xpos = fopen(file_path1.c_str(), "w");
	pfile_ypos = fopen(file_path2.c_str(), "w");
	pfile_zpos = fopen(file_path3.c_str(), "w");
	pfile_ten = fopen(file_path4.c_str(), "w");
	pfile_ten_line = fopen(file_path5.c_str(), "w");
	pfile_debug = fopen(file_path6.c_str(), "w");
}

void Line::CloseOutputFilesASCII(void)
{
	fclose(pfile_xpos);
	fclose(pfile_ypos);
	fclose(pfile_zpos);
	fclose(pfile_ten);
	fclose(pfile_ten_line);
	fclose(pfile_debug);
}

void Line::WriteOut(double t)
{
	int ii;

	fprintf(pfile_xpos, "%f    ", t);
	for (ii = 0; ii < this->N; ii = ii + 1)
		fprintf(pfile_xpos, "%f    ", this->pos(ii, 0));
	fprintf(pfile_xpos, "\n");

	fprintf(pfile_ypos, "%f    ", t);
	for (ii = 0; ii < this->N; ii = ii + 1)
		fprintf(pfile_ypos, "%f    ", this->pos(ii, 1));
	fprintf(pfile_ypos, "\n");

	fprintf(pfile_zpos, "%f    ", t);
	for (ii = 0; ii < this->N; ii = ii + 1)
		fprintf(pfile_zpos, "%f    ", this->pos(ii, 2));
	fprintf(pfile_zpos, "\n");

	fprintf(pfile_ten, "%f    %f    %f    %f    %f    %f    %f \n", t, ten_1(0, 0), ten_1(1, 0), ten_1(2, 0), ten_N(0, 0), ten_N(1, 0), ten_N(2, 0));

	if (t == 0.0)
	{
		fprintf(pfile_ten_line, "%f    ", t);
		for (ii = 0; ii < this->N; ii = ii + 1)
			fprintf(pfile_ten_line, "%f    ", this->s(ii, 0));
		fprintf(pfile_ten_line, "\n");
	}
	fprintf(pfile_ten_line, "%f    ", t);
	for (ii = 0; ii < this->N; ii = ii + 1)
		fprintf(pfile_ten_line, "%f    ", this->T(ii, 0));
	fprintf(pfile_ten_line, "\n");
}
