
#include <armadillo>
#include <string>
#include "WinchiesController.hpp"
#include "../Simulations/Simulation.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"

WinchieController::WinchieController(void){
	nWinchies=0; 
}

WinchieController::WinchieController(int n, Winchie** Ws, Simulation* pIncSim){
	nWinchies=n; 
	Winchies = Ws; 
	pSim = pIncSim;
}

void WinchieController::ReadPropertiesASCII(FILE* pFile){

	// Declare variables
	char buffer_line [1000];
	double dtemp;
	int itemp, itemp2;

	//Ignoro las tres primeras lineas, donde pone "Controller parameters"
	for(int ii=0; ii<3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile);
	}

	fscanf(pFile, "%d %[^\n]\n", &indBody, buffer_line); indBody = indBody - 1; 

	//Leo todo
	if (fscanf(pFile, "%lf", &Ac) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Bc) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Cc) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Dc) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	for(int ii=0; ii<3; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the lead controller gains \n";
			throw ValueError(ss.str());
		}
		Kc(ii,0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	for(int ii=0; ii<3; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the integral time gains \n";
			throw ValueError(ss.str());
		}
		Ki(ii,0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	if (fscanf(pFile, "%lf", &Af) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Bf) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Cf) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	if (fscanf(pFile, "%lf", &Ar) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the reference state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Br) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the reference state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Cr) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the reference state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Dr) != 1)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the reference state-space model parameters \n";
		throw ValueError(ss.str());
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	for(int ii=0; ii<3; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the vector of reference positions \n";
			throw ValueError(ss.str());
		}
		ur(ii,0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%lf %[^\n]\n", &Kw, buffer_line);



	//Ignoro las tres primeras lineas, donde pone "Inversor block inputs"
	for(int ii=0; ii<3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile);
	}

	fscanf(pFile, "%lf %[^\n]\n", &T_max, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &T_min, buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_x_pos = arma::zeros<arma::uvec>(itemp,1);
	coef_x_pos = arma::zeros(itemp,1);
	std::cout << "itemp = " << itemp << std::endl;
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::cout << "itemp2 = " << itemp2 << std::endl;
			std::stringstream ss;
			ss << "An error ocurred when trying to read the indices of lines for positive X force \n";
			throw ValueError(ss.str());
		}
		ind_x_pos(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the coefficients of lines for positive X force \n";
			throw ValueError(ss.str());
		}
		coef_x_pos(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_x_neg = arma::zeros<arma::uvec>(itemp,1);
	coef_x_neg = arma::zeros(itemp,1);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the indices of lines for negative X force \n";
			throw ValueError(ss.str());
		}
		ind_x_neg(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the coefficients of lines for negative X force \n";
			throw ValueError(ss.str());
		}
		coef_x_neg(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_y_pos = arma::zeros<arma::uvec>(itemp,1);
	coef_y_pos = arma::zeros(itemp,1);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the indices of lines for positive Y force \n";
			throw ValueError(ss.str());
		}
		ind_y_pos(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the coefficients of lines for positive Y force \n";
			throw ValueError(ss.str());
		}
		coef_y_pos(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_y_neg = arma::zeros<arma::uvec>(itemp,1);
	coef_y_neg = arma::zeros(itemp,1);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the indices of lines for negative Y force \n";
			throw ValueError(ss.str());
		}
		ind_y_neg(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the coefficients of lines for negative Y force \n";
			throw ValueError(ss.str());
		}
		coef_y_neg(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_g_pos = arma::zeros<arma::uvec>(itemp,1);
	coef_g_pos = arma::zeros(itemp,1);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the indices of lines for positive YAW force \n";
			throw ValueError(ss.str());
		}
		ind_g_pos(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the coefficients of lines for positive YAW force \n";
			throw ValueError(ss.str());
		}
		coef_g_pos(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_g_neg = arma::zeros<arma::uvec>(itemp,1);
	coef_g_neg = arma::zeros(itemp,1);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the indices of lines for negative YAW force \n";
			throw ValueError(ss.str());
		}
		ind_g_neg(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for(int ii=0; ii<itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the coefficients of lines for negative YAW force \n";
			throw ValueError(ss.str());
		}
		coef_g_neg(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	int ne = ceil(pSim->simulationTime/pSim->maxTimeStep);
	error = arma::zeros(ne,3);

	T_min = T_min*pSim->gravity*1000; T_max = T_max*pSim->gravity*1000;
	std::cout << "Winches controller T_min = " << T_min << " N" << std::endl;
	std::cout << "Winches controller T_max = " << T_max << " N" << std::endl;

	coef_x_pos = arma::pow(coef_x_pos,-1); coef_x_pos = coef_x_pos/arma::accu(coef_x_pos);
	coef_x_neg = arma::pow(coef_x_neg,-1); coef_x_neg = coef_x_neg/arma::accu(coef_x_neg);
	coef_y_pos = arma::pow(coef_y_pos,-1); coef_y_pos = coef_y_pos/arma::accu(coef_y_pos);
	coef_y_neg = arma::pow(coef_y_neg,-1); coef_y_neg = coef_y_neg/arma::accu(coef_y_neg);
	coef_g_pos = arma::pow(coef_g_pos,-1); coef_g_pos = coef_g_pos/arma::accu(coef_g_pos);
	coef_g_neg = arma::pow(coef_g_neg,-1); coef_g_neg = coef_g_neg/arma::accu(coef_g_neg);

	for(int ii=0; ii<nWinchies; ii=ii+1){
		Winchies[ii]->tau = T_min * Winchies[ii]->radius;
		std::cout << "  Tension required on winch " << ii+1 << " is T = " << arma::as_scalar(T_min/9800) << " tons" << std::endl;
	}

}

void WinchieController::controlWinchies(void){

	arma::mat e1, pos; 

	xr = Ar*xr + Br*ur; yr = Cr*xr + Dr*ur;

	pos = pSim->pBodies[indBody]->pos; 
	yb(0,0) = pos(0,0);  yb(1,0) = pos(1,0); yb(2,0) = pos(5,0);

	xf = Af*xf + Bf*yb; yf = Cf*xf;

	error.row(k) = (yr-yf).t();
	e1 = Ki%arma::trapz(error.rows(0,k)).t() - yb;

	xc = Ac*xc + Bc*e1; yc = Kc%(Cc*xc + Dc*e1);

	double Fx, Fy, Fg;
	Fx = arma::as_scalar(yc(0,0)); Fy = arma::as_scalar(yc(1,0)); Fg = arma::as_scalar(yc(2,0));
	arma::mat f = arma::zeros(nWinchies,3);
	arma::uvec ind_0 = arma::zeros<arma::uvec>(1);
	if (Fx>0) f.submat(ind_x_pos,ind_0) = coef_x_pos;
	if (Fx<0) f.submat(ind_x_neg,ind_0) = coef_x_neg;
	if (Fy>0) f.submat(ind_y_pos,ind_0+1) = coef_y_pos;
	if (Fy<0) f.submat(ind_y_neg,ind_0+1) = coef_y_neg;
	if (Fg>0) f.submat(ind_g_pos,ind_0+2) = coef_g_pos;
	if (Fg<0) f.submat(ind_g_neg,ind_0+2) = coef_g_neg;

	arma::mat T  = arma::clamp(T_min+Kw*f*arma::abs(yc), T_min, T_max);

	for(int ii=0; ii<nWinchies; ii=ii+1){

		Winchies[ii]->tau = T(ii,0) * Winchies[ii]->radius;
		std::cout << "  Tension required on winch " << ii+1 << " is T = " << arma::as_scalar(T(ii,0)/9800) << " tons" << std::endl;
	}

	k = k + 1;
	
}
