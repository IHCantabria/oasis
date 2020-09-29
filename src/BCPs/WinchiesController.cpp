
#include <armadillo>
#include <string>
#include "WinchiesController.hpp"
#include "../Simulations/Simulation.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"

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

	fscanf(pFile, "%d %[^\n]\n", &inversor_flag);
	fscanf(pFile, "%lf %[^\n]\n", &T_max, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &T_min, buffer_line);

	//Ignoro la linea donde pone "Straigt lines inversor inputs"
	fgets(buffer_line, sizeof(buffer_line), pFile);

	fscanf(pFile, "%d %[^\n]\n", &num_sol);

	//Ignoro la linea donde pone "Coefficients inversor inputs"
	fgets(buffer_line, sizeof(buffer_line), pFile);

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

	T = arma::ones(nWinchies,1)*T_min;

	coef_x_pos = arma::pow(coef_x_pos,-1); coef_x_pos = coef_x_pos/arma::accu(coef_x_pos);
	coef_x_neg = arma::pow(coef_x_neg,-1); coef_x_neg = coef_x_neg/arma::accu(coef_x_neg);
	coef_y_pos = arma::pow(coef_y_pos,-1); coef_y_pos = coef_y_pos/arma::accu(coef_y_pos);
	coef_y_neg = arma::pow(coef_y_neg,-1); coef_y_neg = coef_y_neg/arma::accu(coef_y_neg);
	coef_g_pos = arma::pow(coef_g_pos,-1); coef_g_pos = coef_g_pos/arma::accu(coef_g_pos);
	coef_g_neg = arma::pow(coef_g_neg,-1); coef_g_neg = coef_g_neg/arma::accu(coef_g_neg);

	posAnchG = arma::zeros(3,nWinchies);
	posFairL = arma::zeros(3,nWinchies);
	int indAnch, indFair;

	for(int ii=0; ii<nWinchies; ii=ii+1){

		Winchies[ii]->tau = T_min * Winchies[ii]->radius;

		indFair = Winchies[ii]->LineBCP - 1; indAnch = 0;
		if (indFair == 0) indAnch = 1;
		posAnchG.col(ii) = Winchies[ii]->LineW->pLineBcps[indAnch]->pos;
		posFairL.col(ii) = Winchies[ii]->LineW->pLineBcps[indFair]->pos;

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

	k = k + 1;

	inversorBlock();
	
}


void WinchieController::inversorBlock(void){

	if (inversor_flag == 1) {

		arma::mat rotMat = pSim->pBodies[indBody]->rotMat;
		arma::mat posBody = pSim->pBodies[indBody]->pos;
		arma::mat posFairG_temp, rG;
		double alpha, beta, dx, dy, dz, rx, ry;
		for(int ii=0; ii<nWinchies; ii=ii+1){
			rG = rotMat*posFairL.col(ii);
			posFairG_temp = posBody + rG;
			dx = arma::as_scalar(posAnchG(0,ii)-posFairG_temp(0,0));
			dy = arma::as_scalar(posAnchG(1,ii)-posFairG_temp(1,0));
			dz = arma::as_scalar(posAnchG(2,ii)-posFairG_temp(2,0));
			rx = arma::as_scalar(rG(0,0));
			ry = arma::as_scalar(rG(1,0));
			alpha = atan2(dy,dx); beta = atan2(-dz,sqrt(dx*dx+dy*dy));
			Aeq(0,ii) = cos(alpha)*cos(beta); Aeq(1,ii) = sin(alpha)*cos(beta);
			Aeq(2,ii) = (rx*sin(alpha)-ry*cos(alpha))*cos(beta);
		}
		beq = Kw*yc;

		status_flag = true; int nIterMax = 5; int k = 1;
		arma::mat x0 = T; arma::mat x;
		while (status_flag && (k<=nIterMax)) {
			x = find_tension(x0); x0 = x; k = k + 1;
		}
		if (status_flag) {
			k = 1;
			while (status_flag && (k<=nIterMax)) {
				x0 = T_min + (T_max-T_min)*arma::randu(arma::size(T));
				x = find_tension(x0); k = k + 1;
			}
		}
		if (status_flag) {
			x = find_tension(T);
		}
		T = x;

	} else if (inversor_flag == 2) {

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
		T  = arma::clamp(T_min+Kw*f*arma::abs(yc), T_min, T_max);

	} else {
		std::stringstream ss;
		ss << "Inversor flag not available. \n";
		throw ValueError(ss.str());
	}

	for(int ii=0; ii<nWinchies; ii=ii+1){
		Winchies[ii]->tau = T(ii,0) * Winchies[ii]->radius;
		//std::cout << "  Tension required on winch " << ii+1 << " is T = " << arma::as_scalar(T(ii,0)/9800) << " tons" << std::endl;
	}

}

arma::mat WinchieController::find_tension(arma::mat x0){
	// Initiallize the status flag to success status
	status_flag = false;
	// Initiallize output
	arma::mat x;
	// Check if the problem has solutions, if it does not, return an error
	if (nWinchies<3) {
		std::stringstream ss;
		ss << "ERROR: No solution is possible. \n";
		throw ValueError(ss.str());
    }
    // If the number of equations is equal to the number of variables, solve
    // the square system of equations
    if (nWinchies==3) {
		bool status = solve(x, Aeq, beq);
		// If something goes wrong simply use the initial guess and return a failure status flag
		if (status){
			x = x0;
			status_flag = true;
		}
		// Clamp the solution if needed and return a failure status flag
		arma::mat x_old = x;
		x = arma::clamp(x,T_min,T_max);
		if (arma::norm(x-x0)>1e-6) status_flag = true;
    }
    if (nWinchies>3) {
    	arma::umat ind = comb_n_k(nWinchies,3); int nC = ind.n_rows;
    	arma::mat xx = arma::zeros(nWinchies,nC);
    	arma::mat err = arma::zeros(nC,1);
    	arma::uvec ind_i, indC_i;
    	arma::mat A, Ac, temp, x_i;
    	for(int i=0; i<nC; i=i+1){
    		ind_i = ind.row(i);
    		indC_i = comp_ind(nWinchies,ind_i);
    		A = Aeq.cols(ind_i);
    		Ac = Aeq.cols(indC_i);
    		x_i = x0;
    		bool status = solve(temp, A, beq-Ac*x0.rows(indC_i));
    		if (!status){
    			temp = arma::clamp(temp,T_min,T_max);
    			x_i.rows(ind_i) = temp;
    		}
    		xx.col(i) = x_i;
    		err(i,0) = arma::as_scalar(arma::max(arma::abs(Aeq*x_i-beq)));
    	}
    	arma::uvec ind_great = arma::find(err<1e-12);
    	if (ind_great.is_empty()) {
    		double err0 = arma::as_scalar(arma::max(arma::abs(Aeq*x0-beq)));
    		double err_min = err.min();
    		arma::uword ii = err.index_min();
    		if (err0<err_min){
    			x = x0; status_flag = true;
    		} else {
				x = xx.col(ii); status_flag = true;
    		}
    	} else {
    		xx = xx.cols(ind_great);
    		x = arma::mean(xx,1);
    		int nx = ind_great.n_elem;
    		if (nx>1) {
    			arma::mat coefs = 2*arma::randu(nx,num_sol)-1;
    			coefs = coefs/arma::sum(coefs);	xx = xx*coefs;
    			arma::uvec ind_max = arma::find(arma::max(xx)>T_max); 
    			arma::uvec ind_min = arma::find(arma::min(xx)<T_min);
    			arma::uvec ind_bounds = comp_ind(num_sol,arma::unique(arma::join_vert(ind_max,ind_min)));
    			xx = xx.cols(ind_bounds); xx = arma::join_horiz(xx,x);
    			arma::mat max_dif = arma::max(arma::abs(xx-T)); 
    			double minimax_dif = arma::as_scalar(arma::min(max_dif));
    			arma::uvec ind_lowTchange = arma::find(max_dif<minimax_dif*1.05);
    			if (ind_lowTchange.n_elem>1) {
    				arma::mat temp = arma::sum(xx.cols(ind_lowTchange));
    				arma::uword ind_best = temp.index_min();
    				x = xx.col(ind_best);
    			} else {
    				x = xx.cols(ind_lowTchange);
    			}
    		}
    	}
    }
    return x;
}
