
#include <armadillo>
#include <string>
#include <math.h>
#include "WinchiesController.hpp"
#include "../Simulations/Simulation.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"

WinchieController::WinchieController(void)
{
	nWinchies = 0;
}

WinchieController::WinchieController(int n, Winchie **Ws, Simulation *pIncSim)
{
	nWinchies = n;
	Winchies = Ws;
	pSim = pIncSim;
}

void WinchieController::ReadPropertiesASCII(FILE *pFile)
{

	// Declare variables
	char buffer_line[1000];
	double dtemp;
	int itemp, itemp2;

	// Ignoro las tres primeras lineas, donde pone "Controller parameters"
	for (int ii = 0; ii < 3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile);
	}

	fscanf(pFile, "%d %[^\n]\n", &indBody, buffer_line);
	indBody = indBody - 1;

	// Leo todo
	if (fscanf(pFile, "%lf", &Ac) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Bc) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Cc) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Dc) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the controller state-space model parameters \n";
		throw ValueError(ss.str());
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	for (int ii = 0; ii < 3; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the lead controller gains \n";
			throw ValueError(ss.str());
		}
		Kc(ii, 0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	for (int ii = 0; ii < 3; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the integral time gains \n";
			throw ValueError(ss.str());
		}
		Ki(ii, 0) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	if (fscanf(pFile, "%lf", &Af) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Bf) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Cf) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	if (fscanf(pFile, "%lf", &Df) != 1)
	{
		std::stringstream ss;
		ss << "An error occurred when trying to read the filter state-space model parameters \n";
		throw ValueError(ss.str());
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%lf %[^\n]\n", &Kw, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &time_ini, buffer_line);
	fscanf(pFile, "%ld %[^\n]\n", &reference_flag, buffer_line);

	if (reference_flag < 1)
	{
		// Ignoro la linea, donde pone "State space reference position"
		fgets(buffer_line, sizeof(buffer_line), pFile);

		if (fscanf(pFile, "%lf", &Ar) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the reference state-space model parameters \n";
			throw ValueError(ss.str());
		}
		if (fscanf(pFile, "%lf", &Br) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the reference state-space model parameters \n";
			throw ValueError(ss.str());
		}
		if (fscanf(pFile, "%lf", &Cr) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the reference state-space model parameters \n";
			throw ValueError(ss.str());
		}
		if (fscanf(pFile, "%lf", &Dr) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the reference state-space model parameters \n";
			throw ValueError(ss.str());
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);

		for (int ii = 0; ii < 3; ii++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the vector of reference positions \n";
				throw ValueError(ss.str());
			}
			ur(ii, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);

		// Ignoro las lineas de "Discrete points reference position"
		for (int ii = 0; ii < 5; ii++)
		{
			fgets(buffer_line, sizeof(buffer_line), pFile);
		}
	}
	else
	{
		// Ignoro las lineas de "State space reference position"
		for (int ii = 0; ii < 3; ii++)
		{
			fgets(buffer_line, sizeof(buffer_line), pFile);
		}

		// Ignoro la linea, donde pone "Discrete points reference position"
		fgets(buffer_line, sizeof(buffer_line), pFile);

		// Inicializo los datos de la señal de referencia
		t_ref = arma::zeros(reference_flag, 1);
		x_ref = arma::zeros(reference_flag, 1);
		y_ref = arma::zeros(reference_flag, 1);
		yaw_ref = arma::zeros(reference_flag, 1);

		for (int ii = 0; ii < reference_flag; ii++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the times for reference positions \n";
				throw ValueError(ss.str());
			}
			t_ref(ii, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		for (int ii = 0; ii < reference_flag; ii++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the x for reference positions \n";
				throw ValueError(ss.str());
			}
			x_ref(ii, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		for (int ii = 0; ii < reference_flag; ii++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the y for reference positions \n";
				throw ValueError(ss.str());
			}
			y_ref(ii, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
		for (int ii = 0; ii < reference_flag; ii++)
		{
			if (fscanf(pFile, "%lf", &dtemp) != 1)
			{
				std::stringstream ss;
				ss << "An error occurred when trying to read the yaw for reference positions \n";
				throw ValueError(ss.str());
			}
			yaw_ref(ii, 0) = dtemp;
		}
		fscanf(pFile, "%[^\n]\n", buffer_line);
	}

	// Ignoro las tres primeras lineas, donde pone "Inversor block inputs"
	for (int ii = 0; ii < 3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile);
	}

	fscanf(pFile, "%d %[^\n]\n", &inversor_flag, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &T_max, buffer_line, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &T_min, buffer_line, buffer_line);

	// Ignoro la linea donde pone "Straigt lines inversor inputs"
	fgets(buffer_line, sizeof(buffer_line), pFile);

	fscanf(pFile, "%d %[^\n]\n", &nIterMax, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &atol, buffer_line);

	// Ignoro la linea donde pone "Coefficients inversor inputs"
	fgets(buffer_line, sizeof(buffer_line), pFile);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_x_pos = arma::zeros<arma::uvec>(itemp, 1);
	coef_x_pos = arma::zeros(itemp, 1);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::cout << "itemp2 = " << itemp2 << std::endl;
			std::stringstream ss;
			ss << "An error occurred when trying to read the indices of lines for positive X force \n";
			throw ValueError(ss.str());
		}
		ind_x_pos(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the coefficients of lines for positive X force \n";
			throw ValueError(ss.str());
		}
		coef_x_pos(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_x_neg = arma::zeros<arma::uvec>(itemp, 1);
	coef_x_neg = arma::zeros(itemp, 1);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the indices of lines for negative X force \n";
			throw ValueError(ss.str());
		}
		ind_x_neg(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the coefficients of lines for negative X force \n";
			throw ValueError(ss.str());
		}
		coef_x_neg(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_y_pos = arma::zeros<arma::uvec>(itemp, 1);
	coef_y_pos = arma::zeros(itemp, 1);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the indices of lines for positive Y force \n";
			throw ValueError(ss.str());
		}
		ind_y_pos(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the coefficients of lines for positive Y force \n";
			throw ValueError(ss.str());
		}
		coef_y_pos(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_y_neg = arma::zeros<arma::uvec>(itemp, 1);
	coef_y_neg = arma::zeros(itemp, 1);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the indices of lines for negative Y force \n";
			throw ValueError(ss.str());
		}
		ind_y_neg(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the coefficients of lines for negative Y force \n";
			throw ValueError(ss.str());
		}
		coef_y_neg(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_g_pos = arma::zeros<arma::uvec>(itemp, 1);
	coef_g_pos = arma::zeros(itemp, 1);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the indices of lines for positive YAW force \n";
			throw ValueError(ss.str());
		}
		ind_g_pos(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the coefficients of lines for positive YAW force \n";
			throw ValueError(ss.str());
		}
		coef_g_pos(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
	ind_g_neg = arma::zeros<arma::uvec>(itemp, 1);
	coef_g_neg = arma::zeros(itemp, 1);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%d", &itemp2) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the indices of lines for negative YAW force \n";
			throw ValueError(ss.str());
		}
		ind_g_neg(ii) = itemp2 - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
	for (int ii = 0; ii < itemp; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error occurred when trying to read the coefficients of lines for negative YAW force \n";
			throw ValueError(ss.str());
		}
		coef_g_neg(ii) = dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);
}

void WinchieController::SetUpWinchiesController(void)
{

	int ne = ceil(pSim->simulationTime / pSim->maxTimeStep);
	error = arma::zeros(ne, 3);

	T_min = T_min * pSim->gravity * 1000;
	T_max = T_max * pSim->gravity * 1000;
	std::cout << "Winches controller T_min = " << T_min << " N" << std::endl;
	std::cout << "Winches controller T_max = " << T_max << " N" << std::endl;

	T = arma::ones(nWinchies, 1) * T_min;

	coef_x_pos = arma::pow(coef_x_pos, -1);
	coef_x_pos = coef_x_pos / arma::accu(coef_x_pos);
	coef_x_neg = arma::pow(coef_x_neg, -1);
	coef_x_neg = coef_x_neg / arma::accu(coef_x_neg);
	coef_y_pos = arma::pow(coef_y_pos, -1);
	coef_y_pos = coef_y_pos / arma::accu(coef_y_pos);
	coef_y_neg = arma::pow(coef_y_neg, -1);
	coef_y_neg = coef_y_neg / arma::accu(coef_y_neg);
	coef_g_pos = arma::pow(coef_g_pos, -1);
	coef_g_pos = coef_g_pos / arma::accu(coef_g_pos);
	coef_g_neg = arma::pow(coef_g_neg, -1);
	coef_g_neg = coef_g_neg / arma::accu(coef_g_neg);

	posAnchG = arma::zeros(3, nWinchies);
	posFairL = arma::zeros(3, nWinchies);
	int indAnch, indFair;

	for (int ii = 0; ii < nWinchies; ii = ii + 1)
	{

		Winchies[ii]->tau = T_min * Winchies[ii]->radius;

		indFair = Winchies[ii]->LineBCP - 1;
		indAnch = 0;
		if (indFair == 0)
			indAnch = 1;
		posAnchG.col(ii) = Winchies[ii]->LineW->pLineBcps[indAnch]->pos;
		posFairL.col(ii) = Winchies[ii]->LineW->pLineBcps[indFair]->pos;
	}
}

void WinchieController::controlWinchies(double time)
{

	arma::mat pos;

	if (time >= time_ini)
	{
		arma::mat e1;
		// Reference signal
		if (reference_flag < 1)
		{
			xr = Ar * xr + Br * ur;
			yr = Cr * xr + Dr * ur;
		}
		else
		{
			arma::mat temp_input, temp_output;
			temp_input = arma::ones(1) * time;
			arma::interp1(t_ref, x_ref, temp_input, temp_output, "linear", arma::as_scalar(x_ref(reference_flag - 1, 0)));
			yr(0, 0) = arma::as_scalar(temp_output);
			arma::interp1(t_ref, y_ref, temp_input, temp_output, "linear", arma::as_scalar(y_ref(reference_flag - 1, 0)));
			yr(1, 0) = arma::as_scalar(temp_output);
			arma::interp1(t_ref, yaw_ref, temp_input, temp_output, "linear", arma::as_scalar(yaw_ref(reference_flag - 1, 0)));
			yr(2, 0) = arma::as_scalar(temp_output);
		}
		// Extract position from body
		// Aqui habria que meter ruido gausiano para el ruido de los sensores
		pos = pSim->pBodies[indBody]->pos;
		yb(0, 0) = pos(0, 0);
		yb(1, 0) = pos(1, 0);
		yb(2, 0) = pos(5, 0);
		// First order filter
		xf = Af * xf + Bf * yb;
		yf = Cf * xf + Df * yb;
		// Compute error
		error.row(k) = (yr - yf).t();
		// Integrate error
		e1 = Ki % arma::trapz(error.rows(0, k)).t() - yf;
		xc = Ac * xc + Bc * e1;
		yc = Kc % (Cc * xc + Dc * e1);
		k = k + 1;
		inversorBlock();
	}

	for (int ii = 0; ii < nWinchies; ii = ii + 1)
	{
		Winchies[ii]->tau = T(ii, 0) * Winchies[ii]->radius;
	}
}

void WinchieController::inversorBlock(void)
{

	if (inversor_flag == 1)
	{

		arma::mat rotMat = pSim->pBodies[indBody]->rotMat;
		arma::mat posBody = pSim->pBodies[indBody]->pos.rows(0, 2);
		arma::mat posFairG_temp, rG;
		arma::mat Aeq = arma::zeros(3, nWinchies);
		double alpha, beta, dx, dy, dz, rx, ry, rz;

		for (int ii = 0; ii < nWinchies; ii = ii + 1)
		{
			rG = rotMat * posFairL.col(ii);
			posFairG_temp = posBody + rG;
			dx = arma::as_scalar(posAnchG(0, ii) - posFairG_temp(0, 0));
			dy = arma::as_scalar(posAnchG(1, ii) - posFairG_temp(1, 0));
			dz = arma::as_scalar(posAnchG(2, ii) - posFairG_temp(2, 0));
			rx = arma::as_scalar(rG(0, 0));
			ry = arma::as_scalar(rG(1, 0));
			rz = arma::as_scalar(rG(2, 0));
			alpha = atan2(dy, dx);
			beta = atan2(-dz, sqrt(dx * dx + dy * dy));
			Aeq(0, ii) = cos(alpha) * cos(beta);
			Aeq(1, ii) = sin(alpha) * cos(beta);
			Aeq(2, ii) = (ry * sin(beta) - rz * sin(alpha) * cos(beta)) * rotMat(0, 2) +
						 (rz * cos(alpha) * cos(beta) - rx * sin(beta)) * rotMat(1, 2) +
						 (rx * sin(alpha) * cos(beta) - ry * cos(alpha) * cos(beta)) * rotMat(2, 2);
		}

		arma::mat beq = Kw * yc;

		arma::mat A = arma::join_vert(arma::eye(nWinchies, nWinchies), arma::ones(1, nWinchies) / nWinchies);
		arma::mat b = arma::join_vert(T, arma::zeros(1, 1));
		arma::mat C = Aeq;
		arma::mat d = beq;
		arma::mat AA = A.t() * A;
		arma::mat invAA = arma::inv(AA);

		arma::mat x = invAA * (A.t() * b - C.t() * arma::inv(C * invAA * C.t()) * (C * invAA * A.t() * b - d));

		int kk = 1;
		arma::mat CC = C.t() * arma::inv(C * C.t());
		while ((arma::any(arma::any(x > T_max + atol)) || arma::any(arma::any(x < T_min - atol))) && (kk <= nIterMax))
		{
			b = arma::clamp(x, T_min, T_max);
			x = b - CC * (C * b - d);
			kk = kk + 1;
		}

		if (kk > nIterMax)
		{
			std::cout << " WARNING: Tensions clamped on winches controller! " << std::endl;
		}

		T = arma::clamp(x, T_min, T_max);
	}
	else if (inversor_flag == 2)
	{

		double Fx, Fy, Fg;
		Fx = arma::as_scalar(yc(0, 0));
		Fy = arma::as_scalar(yc(1, 0));
		Fg = arma::as_scalar(yc(2, 0));
		arma::mat f = arma::zeros(nWinchies, 3);
		arma::uvec ind_0 = arma::zeros<arma::uvec>(1);
		if (Fx > 0)
			f.submat(ind_x_pos, ind_0) = coef_x_pos;
		if (Fx < 0)
			f.submat(ind_x_neg, ind_0) = coef_x_neg;
		if (Fy > 0)
			f.submat(ind_y_pos, ind_0 + 1) = coef_y_pos;
		if (Fy < 0)
			f.submat(ind_y_neg, ind_0 + 1) = coef_y_neg;
		if (Fg > 0)
			f.submat(ind_g_pos, ind_0 + 2) = coef_g_pos;
		if (Fg < 0)
			f.submat(ind_g_neg, ind_0 + 2) = coef_g_neg;
		T = arma::clamp(T_min + Kw * f * arma::abs(yc), T_min, T_max);
	}
	else
	{
		std::stringstream ss;
		ss << "Inversor flag not available. \n";
		throw ValueError(ss.str());
	}
}

void WinchieController::OpenOutputFilesASCII(std::string path)
{

	char buffer1[50];
	int nn1 = sprintf(buffer1, "WinchesTensions.txt");
	std::string file_path1 = JoinPath(path, buffer1);
	pfile_TW = fopen(file_path1.c_str(), "w");
	if (pfile_TW == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn1 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}

	char buffer2[50];
	int nn2 = sprintf(buffer2, "ControlForce.txt");
	std::string file_path2 = JoinPath(path, buffer2);
	pfile_FC = fopen(file_path2.c_str(), "w");
	if (pfile_FC == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn2 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}

	char buffer3[50];
	int nn3 = sprintf(buffer3, "ReferencePosition.txt");
	std::string file_path3 = JoinPath(path, buffer3);
	pfile_RP = fopen(file_path3.c_str(), "w");
	if (pfile_RP == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn3 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}

	char buffer4[50];
	int nn4 = sprintf(buffer4, "WinchedLinesLengths.txt");
	std::string file_path4 = JoinPath(path, buffer4);
	pfile_LL = fopen(file_path4.c_str(), "w");
	if (pfile_LL == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn4 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}
}

void WinchieController::CloseOutputFilesASCII(void)
{

	fclose(pfile_TW);
	fclose(pfile_FC);
	fclose(pfile_RP);
	fclose(pfile_LL);
}

void WinchieController::WriteOut(double t)
{

	fprintf(pfile_TW, "%f    ", t);
	for (int ii = 0; ii < nWinchies; ii = ii + 1)
		fprintf(pfile_TW, "%f    ", T(ii, 0));
	fprintf(pfile_TW, "\n");

	fprintf(pfile_FC, "%f    ", t);
	for (int ii = 0; ii < 3; ii = ii + 1)
		fprintf(pfile_FC, "%f    ", Kw * yc(ii, 0));
	fprintf(pfile_FC, "\n");

	fprintf(pfile_RP, "%f    ", t);
	for (int ii = 0; ii < 3; ii = ii + 1)
		fprintf(pfile_RP, "%f    ", yr(ii, 0));
	fprintf(pfile_RP, "\n");

	fprintf(pfile_LL, "%f    ", t);
	for (int ii = 0; ii < nWinchies; ii = ii + 1)
		fprintf(pfile_LL, "%f    ", Winchies[ii]->LineW->L * Winchies[ii]->LineW->dL / Winchies[ii]->LineW->dL0);
	fprintf(pfile_LL, "\n");
}
