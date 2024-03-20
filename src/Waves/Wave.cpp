
#include <armadillo>
#include <tuple>
#include "Wave.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"
#include <math.h>

Wave::Wave(Simulation *pSimInc, double H, double T, double D)
{
	height = H;
	period = T;
	heading = D * arma::datum::pi / 180.0;
	pSim = pSimInc;
	simulationTime = pSim->simulationTime;
	gravity = pSim->gravity;
	waterDepth = abs(pSim->waterDepth);
	if (period == 0)
	{
		std::stringstream ss;
		ss << "The wave period cannot be zero. \n";
		throw ValueError(ss.str());
	}
	lambda_peak = solve_lambda(period);
}

void Wave::CheckBreakingWave(void)
{
	std::cout << "--> Checking wave breaking limits." << std::endl;

	double lambda = solve_lambda(period);
	double hL = waterDepth / lambda;

	double hmax;

	if (hL > 0.108)
	{
		hmax = 0.142 * tanh(2.0 * arma::datum::pi * hL) * lambda;
	}
	else if (hL <= 0.108)
	{
		hmax = 0.78 * waterDepth;
	}

	if (height > hmax)
	{
		std::stringstream ss;
		ss << "The water is outside of the breaking limits. \n";
		throw ValueError(ss.str());
	}
	else
	{
		std::cout << "--> Wave within wave breaking limits." << std::endl;
	}
}

void Wave::GetWaveLengths(void)
{
	std::cout << "--> Getting wave lengths..." << std::endl;

	std::tie(this->lambdas,
			 this->k,
			 this->kx,
			 this->ky,
			 this->kx_1D,
			 this->ky_1D) = this->GetWaveLengths(this->periods, this->headings);

	std::cout << "--> ... wave lengths computed!" << std::endl;
}

std::tuple<arma::vec, // lambdas
		   arma::vec, // k
		   arma::mat, // kx
		   arma::mat, // ky
		   arma::vec, // kx_1D
		   arma::vec> // ky_1D
Wave::GetWaveLengths(arma::vec periods, arma::vec headings)
{
	double T;
	arma::vec lambdas = arma::zeros(size(periods));
	arma::vec k = lambdas;
	for (int ii = 0; ii < num_comps; ii++)
	{
		T = arma::as_scalar(periods(ii));
		if (T == 0)
		{
			lambdas(ii) = 0.0;
			k(ii) = 0.0;
		}
		else if (std::isinf(T))
		{
			lambdas(ii) = 0.0;
			k(ii) = 0.0;
		}
		else
		{
			lambdas(ii) = solve_lambda(T);
			k(ii) = 2.0 * arma::datum::pi / lambdas(ii);
		}
	}

	int num_comps = periods.n_elem;
	int num_headings = headings.n_elem;
	arma::mat kx = arma::zeros(num_comps, num_headings);
	arma::mat ky = arma::zeros(num_comps, num_headings);
	arma::vec kx_1D = arma::zeros(num_comps);
	arma::vec ky_1D = arma::zeros(num_comps);
	for (int ii = 0; ii < num_comps; ii++)
	{
		// The 1D wave number is computed with the main wave heading
		kx_1D(ii) = k(ii) * cos(this->heading);
		ky_1D(ii) = k(ii) * sin(this->heading);
		for (int jj = 0; jj < num_headings; jj++)
		{
			kx(ii, jj) = k(ii) * cos(headings(jj));
			ky(ii, jj) = k(ii) * sin(headings(jj));
		}
	}

	return std::make_tuple(lambdas, k, kx, ky, kx_1D, ky_1D);
}

double Wave::solve_lambda(double T)
{
	double lambda = gravity * T * T / (2.0 * arma::datum::pi);
	double f = f_lambda(lambda, T);
	double df = 0;
	int ii = 0;
	int nIterMax = 100;
	double atol = 1e-6;
	while ((abs(f) > atol) && (ii < nIterMax))
	{
		df = df_lambda(lambda, T);
		lambda = lambda - f / df;
		f = f_lambda(lambda, T);
		ii++;
	}

	if (abs(f) > atol)
	{
		arma::uvec ind = arma::find(T == periods);
		double H = arma::accu(amplitudes.rows(ind));
		std::stringstream ss;
		ss << "Convergence failed for computing wave length for period T = " << T << " s. \n";
		throw ValueError(ss.str());
	}
	return lambda;
}

double Wave::f_lambda(double lambda, double T)
{
	return gravity * T * T / (2.0 * arma::datum::pi) * tanh(2.0 * arma::datum::pi * waterDepth / lambda) - lambda;
}

double Wave::df_lambda(double lambda, double T)
{
	return -gravity * T * T * waterDepth / (lambda * lambda * cosh(2.0 * arma::datum::pi * waterDepth / lambda)) - 1.0;
}

void Wave::GetFreeSurface(void)
{
	std::cout << "----> Computing Free Surface" << std::endl;
	this->t_FS = arma::linspace(0.0, this->simulationTime, this->num_points);
	this->eta_FS = GetFreeSurface(this->amplitudes, this->phases, this->num_points);
	std::cout << "----> Free Surface Computed" << std::endl;
}

arma::vec Wave::GetFreeSurface(arma::mat amplitudes, arma::mat phases, int num_points)
{
	arma::mat Cm = arma::sum(amplitudes % arma::cos(phases), 1);
	arma::mat Sm = arma::sum(amplitudes % arma::sin(phases), 1);
	arma::mat Am = arma::sqrt(arma::pow(Cm, 2) + arma::pow(Sm, 2));
	arma::mat PHIm = arma::atan2(Sm, Cm);
	arma::cx_mat iPHIm(arma::zeros(size(PHIm)), PHIm);
	arma::cx_mat Y1 = (0.5 * num_points * Am) % arma::exp(-iPHIm);
	Y1(0, 0) = 0.0;
	arma::cx_mat Y;
	if (num_points % 2 == 0)
	{
		Y = arma::join_vert(arma::conj(Y1), arma::flipud(Y1.rows(1, num_comps - 2)));
	}
	else
	{
		Y = arma::join_vert(arma::conj(Y1), arma::flipud(Y1.rows(1, num_comps - 1)));
	}
	arma::cx_mat eta_cx = arma::ifft(Y);
	double imag = arma::as_scalar(arma::sum(arma::abs(arma::imag(eta_cx)), 0));
	if (imag > 1e-9 * num_points)
	{
		std::stringstream ss;
		ss << "Something went wrong with the ifft. imag = " << imag << " \n";
		throw ValueError(ss.str());
	}
	arma::vec eta = arma::real(eta_cx);
	return eta;
}

arma::vec Wave::GetFreeSurface(double time, arma::vec x, arma::vec y)
{

	arma::vec eta = arma::sum(amplitudes * ones(size(amplitudes.t())) * arma::cos(phases * ones(size(x.t())) + kx * x.t() + ky * y.t() - time * ang_freqs * ones(size(x.t())))).t();

	return eta;
}

arma::vec Wave::GetPressure(double time, arma::vec x, arma::vec y, arma::vec z)
{
	// std::cout << "--> Getting Pressure " << std::endl;

	arma::vec dPhidt, dPhidx, dPhidy, dPhidz, pressure;
	arma::vec A = amplitudes, W = ang_freqs;
	double H = waterDepth;
	arma::mat onesNp = arma::ones(size(x.t()));
	arma::mat aux = phases * onesNp + kx * x.t() + ky * y.t() - time * ang_freqs * onesNp;

	dPhidt = arma::sum(-((A % arma::pow(W, 2)) * onesNp) % arma::cos(aux) % arma::cosh(k * (H * onesNp + z.t())) /
					   ((k % arma::sinh(H * k)) * onesNp))
				 .t();
	dPhidx = arma::sum(((A % kx % W) * onesNp) % arma::cos(aux) % arma::cosh(k * (H * onesNp + z.t())) /
					   ((k % arma::sinh(H * k)) * onesNp))
				 .t();
	dPhidy = arma::sum(((A % ky % W) * onesNp) % arma::cos(aux) % arma::cosh(k * (H * onesNp + z.t())) /
					   ((k % arma::sinh(H * k)) * onesNp))
				 .t();
	dPhidz = arma::sum(((A % W) * onesNp % arma::sin(aux)) % arma::sinh(k * (H * onesNp + z.t())) /
					   (arma::sinh(H * k) * onesNp))
				 .t();

	pressure = -pSim->waterDensity * (pSim->gravity * z + dPhidt + 0.5 * (dPhidx % dPhidx + dPhidy % dPhidy + dPhidz % dPhidz));

	return pressure;
}

void Wave::SetSinglePiece(void)
{
	num_pieces = 1;

	time_ref = arma::vec(1);
	time_ini = arma::vec(1);
	time_end = arma::vec(1);
	amplitudes_piece = arma::field<arma::mat>(1);
	phases_piece = arma::field<arma::mat>(1);

	time_ref(0) = 0.0;
	time_ini(0) = 0.0;
	time_end(0) = simulationTime;

	num_points_piece = num_points;
	num_comps_piece = num_comps;
	df_piece = df;
	freqs_piece = freqs;
	ang_freqs_piece = ang_freqs;
	periods_piece = periods;
	dw_piece = dw;
	amplitudes_piece(0) = amplitudes;
	phases_piece(0) = phases;
	lambdas_piece = lambdas;
	k_piece = k;
	kx_piece = kx;
	ky_piece = ky;
	kx_1D_piece = kx_1D;
	ky_1D_piece = ky_1D;
	headings_piece = headings;
	num_headings_piece = num_headings;
}

void Wave::WriteOut(std::string path)
{
	char buffer1[50];
	int nn1 = sprintf(buffer1, "WaveSpectrum.txt");
	std::string file_path1 = JoinPath(path, buffer1);
	pfile_SPEC = fopen(file_path1.c_str(), "ang_freqs");
	if (pfile_SPEC == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn1 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}
	for (int ii = 0; ii < num_comps; ii = ii + 1)
		fprintf(pfile_SPEC, "%f    %f  \n", freqs(ii, 0), S_w(ii, 0));
	fclose(pfile_SPEC);

	char buffer2[50];
	int nn2 = sprintf(buffer2, "WaveTimeSeries.txt");
	std::string file_path2 = JoinPath(path, buffer2);
	pfile_TIME = fopen(file_path2.c_str(), "ang_freqs");
	if (pfile_TIME == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: " << nn2 << "\n    ->Dir: " << path << std::endl;
		throw IOError(ss.str());
	}

	for (int ii = 0; ii < num_points; ii = ii + 1)
		fprintf(pfile_TIME, "%f    %f  \n", t_FS(ii, 0), eta_FS(ii, 0));

	fclose(pfile_TIME);
}

void RegularWave::GetWaveSpectrum(void)
{
	num_comps = 1;
	num_headings = 1;
	df = 1.0 / simulationTime;
	dw = 2.0 * arma::datum::pi * df;
	periods = arma::ones<arma::vec>(1) * period;
	freqs = 1.0 / periods;
	ang_freqs = arma::ones<arma::vec>(1) * (2.0 * arma::datum::pi / period);
	headings = arma::ones<arma::vec>(1) * heading;
	headings = mod(headings, 2.0 * arma::datum::pi);
	amplitudes = arma::ones(1, 1) * height * 0.5;
	phases = arma::zeros(1, 1);
	headings_1D = headings;
	amplitudes_1D = amplitudes;
	phases_1D = phases;
	GetWaveLengths();
}

void IrregularWave::GetWaveSpectrum(void)
{
	// Convert heading step to radians
	dtheta = dtheta * arma::datum::pi / 180;

	if (specType_flag == 1)
	{

		// Find the spectrum width, the maximum frequency step, and the minimum wave time
		// TODO: Review the methods used to find these variables
		FindSpectrumWidth();
		this->df_max = spectrum_width / 11.0;
		double time_sim_min = std::ceil(1 / this->df_max);

		// Update the wave time if necessary
		simulationTime = std::max(simulationTime, time_sim_min);

		// Compute the number of pieces and the time per piece
		GetPiecesNumber();
		this->time_piece = (this->simulationTime + (this->num_pieces - 1) * this->time_gap) / this->num_pieces;

		// Compute the wave for a single piece
		// Compute the time vector
		if (this->num_pieces > 1)
		{
			t_FS = arma::regspace(0, this->dt, simulationTime + 2 * this->time_piece) - this->time_piece;
			std::cout << "----> Computing Wave Spectrum (JONSWAP) in pieces" << std::endl;
			std::cout << "    ----> Number of pieces: " << this->num_pieces << std::endl;
		}
		else
		{
			t_FS = arma::regspace(0, this->dt, simulationTime);
			std::cout << "--> Computing Wave Spectrum (JONSWAP) in a single piece" << std::endl;
		}
		// Compute the frequency spectrum
		num_points = t_FS.n_elem;
		if (num_points % 2 == 0)
		{
			num_points = num_points++;
		}
		num_comps = (num_points + 1) / 2;
		df = 1 / (num_points * dt);
		freqs = arma::linspace(0, (num_comps - 1) * df, num_comps);
		ang_freqs = 2.0 * arma::datum::pi * freqs;
		periods = 1.0 / freqs;
		dw = 2.0 * arma::datum::pi * df;
		GetJonswapSpectrum();
		// Compute the directional spectrum
		if (s > 1.0)
		{
			headings = arma::regspace(heading - 0.5 * arma::datum::pi, dtheta, heading + 0.5 * arma::datum::pi);
			num_headings = headings.n_elem;
			GetSpreadingFunction();
		}
		else
		{
			headings = arma::ones<arma::vec>(1) * heading;
			num_headings = 1;
			g_theta = arma::ones<arma::vec>(1);
			dtheta = 1.0;
		}
		headings = mod(headings, 2.0 * arma::datum::pi);
		// Get the wave amplitudes combining the frequency and directional spectra
		this->amplitudes = arma::zeros(num_comps, num_headings);
		for (int ii = 0; ii < num_comps; ii++)
		{
			for (int jj = 0; jj < num_headings; jj++)
			{
				this->amplitudes(ii, jj) = sqrt(2.0 * S_w(ii) * df * g_theta(jj) * dtheta);
			}
		}
		// Get the wave checked free surface and phases
		this->eta_FS = GetCheckedFreeSurface(this->amplitudes, this->num_points, this->simulationTime);
		// Get the wave lengths
		GetWaveLengths();

		// Compute the time intervals for each piece
		this->time_ref = arma::vec(this->num_pieces);
		this->time_ini = arma::vec(this->num_pieces);
		this->time_end = arma::vec(this->num_pieces);
		for (int ii = 0; ii < this->num_pieces; ii++)
		{
			this->time_ref(ii) = (this->time_piece - this->time_gap) * ii - this->time_piece;
		}
		this->time_ini(0) = 0.0;
		for (int ii = 1; ii < this->num_pieces; ii++)
		{
			this->time_ini(ii) = (ii + 1) * this->time_piece - ii * this->time_gap;
		}
		for (int ii = 0; ii < this->num_pieces - 1; ii++)
		{
			this->time_end(ii) = this->time_ini(ii + 1) - this->time_gap;
		}
		this->time_end(this->num_pieces - 1) = this->simulationTime;

		// Compute the wave frequency spectrum for each piece
		this->amplitudes_piece = arma::field<arma::mat>(this->num_pieces);
		this->phases_piece = arma::field<arma::mat>(this->num_pieces);
		if (this->num_pieces > 1)
		{
			arma::vec time_piece_total = arma::regspace(0, this->dt, 3 * this->time_piece);
			num_points_piece = time_piece_total.n_elem;
			if (num_points_piece % 2 == 0)
			{
				num_points_piece = num_points_piece++;
			}
			num_comps_piece = (num_points_piece + 1) / 2;
			df_piece = 1.0 / (num_points_piece * dt);
			freqs_piece = arma::linspace(0, (num_comps_piece - 1) * df_piece, num_comps_piece);
			ang_freqs_piece = 2.0 * arma::datum::pi * freqs_piece;
			periods_piece = 1.0 / freqs_piece;
			dw_piece = 2.0 * arma::datum::pi * df_piece;
			GetPiecesWaveLengths();
			for (int ii = 0; ii < this->num_pieces; ii++)
			{
				this->amplitudes_piece(ii) = arma::zeros(num_comps_piece, num_headings);
				this->phases_piece(ii) = arma::zeros(num_comps_piece, num_headings);
			}
			for (int jj = 0; jj < num_headings; jj++)
			{
				// Get the wave free surface for this heading
				arma::vec eta_heading = GetFreeSurface(amplitudes.col(jj), phases.col(jj), num_points);
				for (int ii = 0; ii < this->num_pieces; ii++)
				{
					// Get the wave free surface for the times corresponding to this piece
					arma::vec tmp_t = time_piece_total + arma::as_scalar(this->time_ref(ii));
					arma::vec tmp_eta;
					arma::interp1(t_FS, eta_heading, tmp_t, tmp_eta, "*linear");
					// Get the wave frequency amplitude and phases for this piece
					arma::cx_vec yf = arma::fft(tmp_eta) / num_points_piece;
					arma::vec psd = arma::pow(arma::abs(yf.rows(0, num_comps_piece - 1)), 2) / df_piece;
					psd.rows(1, num_comps_piece - 1) = 2.0 * psd.rows(1, num_comps_piece - 1);
					psd(0, 0) = 0.0;
					this->amplitudes_piece(ii).col(jj) = arma::sqrt(2.0 * psd * df_piece);
					this->phases_piece(ii).col(jj) = arma::atan2(arma::imag(yf.rows(0, num_comps_piece - 1)),
																 arma::real(yf.rows(0, num_comps_piece - 1))) -
													 ang_freqs_piece * this->time_ref(ii); // TODO: Check this
				}
			}
			headings_piece = headings;
			num_headings_piece = num_headings;
		}
		else
		{
			SetSinglePiece();
		}

		CutPiecesSpectrumZeros();
		CutSpectrumZeros();

		std::cout << "----> Wave Spectrum Computed" << std::endl;
	}
	else if (specType_flag == 2)
	{
		ReadWaveSpectrumHDF5();
		// TODO: Review the wave input types and implement properly to use with pieces
	}
	else if (specType_flag == 3)
	{
		ReadWaveSpectrumASCII();
		// TODO: Review the wave input types and implement properly to use with pieces
	}
	else
	{
		std::stringstream ss;
		ss << "Error while parsing file: dataWaves.dat; Unexpected spectrum type. \n";
		throw ValueError(ss.str());
	}
}

void IrregularWave::GetJonswapSpectrum(void)
{
	this->S_w = GetJonswapSpectrum(this->freqs, this->height, this->period, this->gamma);
}

arma::vec IrregularWave::GetJonswapSpectrum(arma::vec freqs, double height, double period, double gamma)
{

	// Function jonswap: It returns the jonswap frequency spectra according
	// to the IEC 61400-3 standard. It does it for certain frequencies.
	//
	//
	//  Alvaro Rodriguez Luis
	//  Feb. 2020
	//  IH Cantabria

	double c1 = 0.3125;
	double c2 = 0.287;

	double fp = 1.0 / period;

	arma::vec sigma = arma::zeros(size(freqs));
	for (int ii = 0; ii < num_comps; ii++)
	{
		if (freqs(ii, 0) < fp)
		{
			sigma(ii, 0) = 0.07;
		}
		else
		{
			sigma(ii, 0) = 0.09;
		}
	}

	double cte = c1 * height * height * period * (1.0 - c2 * log(gamma));
	arma::vec alpha = arma::exp(-0.5 * (arma::pow((freqs - fp) / (fp * sigma), 2.0)));
	arma::vec temp1 = arma::pow(freqs / fp, -5.0);
	arma::vec temp2 = arma::exp(-1.25 * arma::pow(freqs / fp, -4.0));
	arma::vec temp3 = arma::exp(log(gamma) * alpha);
	arma::vec S = cte * temp1 % temp2 % temp3;
	S(0) = 0.0;
	return S;
}

void IrregularWave::GetSpreadingFunction(void)
{

	// Function spreading_function: It returns the spreading function
	// cos(dir)^s. Where s ... depends on tp.
	// It does it for certain directions (theta) between -pi/2 and pi/2.
	//
	//
	//  Alvaro Rodriguez Luis
	//  Feb. 2020
	//  IH Cantabria

	arma::vec theta = arma::regspace(-90, dtheta, 90) * arma::datum::pi / 180.0;
	arma::vec G = arma::pow(arma::cos(theta), s);
	double integral = arma::as_scalar(arma::trapz(theta, G, 0));
	g_theta = G / integral;
	integral = arma::as_scalar(arma::trapz(theta, g_theta, 0));
	if (abs(integral - 1) > 1e-6)
	{
		g_theta = g_theta / integral;
	}
}

arma::vec IrregularWave::GetCheckedFreeSurface(arma::mat amplitudes, int num_points, double time)
{
	arma::vec t = arma::linspace(0, time, num_points);
	arma::vec eta;

	int flag = 0; // 0: something is wrong; 1: everything is ok

	int ii = 1;
	int nIterMax = 200;
	while ((flag == 0) && (ii <= nIterMax))
	{
		this->phases = arma::randn(size(amplitudes)) * arma::datum::pi;
		eta = GetFreeSurface(amplitudes, this->phases, num_points);
		flag = CheckPhases(t, eta);
		ii++;
	}
	if (flag == 0)
	{
		std::stringstream ss;
		ss << "It was not possible to find a realistic wave. \n";
		throw ValueError(ss.str());
	}

	return eta;
}

int IrregularWave::CheckPhases(arma::vec t, arma::vec eta)
{
	int flag = 0; // 0: something is wrong; 1: everything is ok

	arma::mat T, H;
	std::tie(T, H) = upcrossing(t, eta);

	double Tmean = arma::as_scalar(arma::mean(T));
	double Hmax = H.max();
	arma::mat Haux = arma::sort(H, "descend");
	int nH = std::round(H.n_elem / 3);
	double Hsig = arma::as_scalar(arma::mean(Haux.rows(0, nH)));

	double TmeanC = period / (sqrt(2.0) * pow(gamma, -0.082));
	double HsigC = height;
	double HmaxC = HsigC * 1.8;

	if ((Tmean <= TmeanC * (1 + rel_tol)) && (Tmean >= TmeanC * (1 - rel_tol)) &&
		(Hsig <= HsigC * (1 + rel_tol)) && (Hsig >= HsigC * (1 - rel_tol)) &&
		(Hmax <= HmaxC * (1 + rel_tol)) && (Hmax >= HmaxC * (1 - rel_tol)))
	{
		flag = 1;
	}

	return flag;
}

void IrregularWave::CutSpectrumZeros(void)
{
	// Crop the spectrum to avoid zeros
	arma::uvec ind_rows = arma::find(S_w >= factor * arma::as_scalar(arma::mean(S_w)));
	arma::uvec ind_cols = arma::find(g_theta >= factor * arma::as_scalar(arma::mean(g_theta)));
	arma::uvec ind_0 = arma::zeros<arma::uvec>(1);

	num_comps = ind_rows.n_elem;
	num_headings = ind_cols.n_elem;

	periods = periods.submat(ind_rows, ind_0);
	freqs = freqs.submat(ind_rows, ind_0);
	ang_freqs = ang_freqs.submat(ind_rows, ind_0);
	headings = headings.submat(ind_cols, ind_0);
	amplitudes = amplitudes.submat(ind_rows, ind_cols);
	phases = phases.submat(ind_rows, ind_cols);
	lambdas = lambdas.submat(ind_rows, ind_0);
	S_w = S_w.submat(ind_rows, ind_0);
	g_theta = g_theta.submat(ind_cols, ind_0);
	k = k.submat(ind_rows, ind_0);
	kx = kx.submat(ind_rows, ind_cols);
	ky = ky.submat(ind_rows, ind_cols);

	amplitudes_1D = amplitudes_1D.submat(ind_rows, ind_0);
	phases_1D = phases_1D.submat(ind_rows, ind_0);
	kx_1D = kx_1D.submat(ind_rows, ind_0);
	ky_1D = ky_1D.submat(ind_rows, ind_0);
}

void IrregularWave::ReadWaveSpectrumHDF5(void)
{
	std::cout << "--> Reading Wave Spectrum (HDF5 format)" << std::endl;
	std::stringstream ss;
	ss << "Method ReadWaveSpectrumHDF5 in class IrregularWave not implemented yet. \n";
	throw NotImplementedError(ss.str());
	std::cout << "----> Wave Spectrum Read" << std::endl;
}

void IrregularWave::ReadWaveSpectrumASCII(void)
{
	std::cout << "--> Reading Wave Spectrum (ASCII format)" << std::endl;

	int nn;
	arma::mat time, eta;
	double dtemp;

	FILE *file_pointer = fopen(file_path.c_str(), "r");

	if (file_pointer == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: dataWaves.dat\n    ->Dir: " << file_path << std::endl;
		throw IOError(ss.str());
	}

	char bufferLine[1000];

	fscanf(file_pointer, "%i", &nn, bufferLine);
	fscanf(file_pointer, "%[^\n]\n", bufferLine);
	time = arma::zeros(nn);
	eta = arma::zeros(nn);

	for (int ii = 0; ii < nn; ii++)
	{
		fscanf(file_pointer, "%lf", &dtemp);
		time(ii, 0) = dtemp;
		fscanf(file_pointer, "%lf", &dtemp);
		eta(ii, 0) = dtemp;
		fscanf(file_pointer, "%[^\n]\n", bufferLine);
	}

	// Close file
	fclose(file_pointer);

	simulationTime = time.max() - time(0, 0);
	dt = arma::as_scalar(time(1, 0) - time(0, 0));
	num_points = nn;
	num_comps = floor(num_points / 2.0) + 1;
	freqs = arma::linspace(0.0, 1.0 / (dt * 2.0), num_comps);
	df = arma::as_scalar(freqs(1, 0) - freqs(0, 0));
	ang_freqs = 2.0 * arma::datum::pi * freqs;
	dw = 2.0 * arma::datum::pi * df;
	periods = 1.0 / freqs;
	dtheta = 1.0;

	arma::cx_vec yf = arma::fft(eta) / num_points;
	phases = arma::atan2(arma::imag(yf.rows(0, num_comps - 1)), arma::real(yf.rows(0, num_comps - 1)));
	arma::vec psd = arma::pow(arma::abs(yf.rows(0, num_comps - 1)), 2) / df;
	psd.rows(1, num_comps - 1) = 2.0 * psd.rows(1, num_comps - 1);
	psd(0, 0) = 0.0;
	amplitudes = arma::sqrt(2.0 * psd * df);

	num_headings = 1;
	g_theta = arma::ones<arma::vec>(1);
	dtheta = 1.0;
	headings = arma::ones<arma::vec>(1) * heading;
	headings_1D = headings;
	amplitudes_1D = amplitudes;
	phases_1D = phases;
	S_w = psd;
	std::cout << "    ----> Compute wave lengths" << std::endl;
	GetWaveLengths();
	std::cout << "    ----> Build free surface" << std::endl;
	GetFreeSurface();
	std::cout << "    ----> Remove unnecesary frequencies" << std::endl;
	CutSpectrumZeros();
	std::cout << "    ----> Write the wave into a single piece. Multiple pieces not implemented yet!" << std::endl;
	SetSinglePiece();
	// TODO: Implement multiple pieces

	std::cout << "----> Wave Spectrum Read" << std::endl;
}

void IrregularWave::FindSpectrumWidth(void)
{
	// Search of FWHM Jonswap spectrum
	if (specType_flag == 1)
	{
		arma::vec fm = arma::linspace(0, 1, 100001);
		arma::mat S = GetJonswapSpectrum(fm, this->height, this->period, this->gamma);
		S(0, 0) = 0.0;

		double Smax = S.max();

		arma::uvec i1 = arma::find(S > Smax / 2, 1, "first");
		double f1 = arma::as_scalar(fm(i1));

		arma::uvec i2 = arma::find(S > Smax / 2, 1, "last");
		double f2 = arma::as_scalar(fm(i2));

		this->spectrum_width = f2 - f1;
	}
	else
	{
		std::stringstream ss;
		ss << "FindSpectrumWidth is not implemented for custom spectra f \n";
		throw ValueError(ss.str());
	}
}

void IrregularWave::GetPiecesNumber(void)
{
	std::cout << "--> Computing number of pieces to divide wave time series..." << std::endl;
	double time_sim = this->simulationTime;
	int tmp_int;

	std::cout << "    --> Computing the gap time between two pieces..." << std::endl;
	this->time_gap = std::min(1.2 * this->period, time_sim * (0.5 - std::sqrt(6.0) / 6.0));
	std::cout << "        --> Time gap: " << this->time_gap << " s" << std::endl;

	std::cout << "    --> Imposing basic constraints..." << std::endl;
	// Constraint given because each piece uses 3*time_piece for fft spectrum computation
	int n_min = 4;
	std::cout << "        --> n_min: 4" << std::endl;
	// Constraint to avoid gap times overlapping
	int n_max = std::floor((time_sim - this->time_gap) / this->time_gap);
	std::cout << "        --> n_max: " << n_max << std::endl;

	// Check the maximum frequency step constraint
	if (this->time_gap < 1.0 / (3.0 * this->df_max))
	{
		std::cout << "    --> Imposing maximum frequency step constraint..." << std::endl;
		tmp_int = std::floor(3.0 * this->df_max * (time_sim - this->time_gap) / (1.0 - 3.0 * this->df_max * this->time_gap));
		n_max = std::min(n_max, tmp_int);
		std::cout << "        --> n_max: " << n_max << std::endl;
	}

	// Check the saving memory constraints
	tmp_int = std::ceil(-(time_sim * std::sqrt(36 * std::pow(this->time_gap, 2) -
											   36 * this->time_gap * time_sim +
											   std::pow(time_sim, 2)) +
						  (18 * this->time_gap * time_sim -
						   18 * std::pow(this->time_gap, 2) -
						   std::pow(time_sim, 2))) /
						(18 * std::pow(this->time_gap, 2)));
	if (n_min < tmp_int)
	{
		std::cout << "    --> Imposing memory saving constraints (n_min)..." << std::endl;
		n_min = tmp_int;
		std::cout << "        --> n_min: " << n_min << std::endl;
	}
	tmp_int = std::ceil(-(time_sim * std::sqrt(36 * std::pow(this->time_gap, 2) -
											   36 * this->time_gap * time_sim +
											   std::pow(time_sim, 2)) -
						  (18 * this->time_gap * time_sim -
						   18 * std::pow(this->time_gap, 2) -
						   std::pow(time_sim, 2))) /
						(18 * std::pow(this->time_gap, 2)));
	if (n_max > tmp_int)
	{
		std::cout << "    --> Imposing memory saving constraints (n_max)..." << std::endl;
		n_max = tmp_int;
		std::cout << "        --> n_max: " << n_max << std::endl;
	}

	// If it is possible, get the optimum number of pieces
	std::cout << "    --> Finding optimum number of pieces..." << std::endl;
	if (n_min <= n_max)
	{
		double cost_min = std::pow(time_sim + (n_min - 1) * this->time_gap, 2) / n_min;
		double cost_max = std::pow(time_sim + (n_max - 1) * this->time_gap, 2) / n_max;
		int n_opt = std::round((time_sim - this->time_gap) / this->time_gap);
		if (n_opt >= n_min && n_opt <= n_max)
		{
			double cost_opt = std::pow(time_sim + (n_opt - 1) * this->time_gap, 2) / n_opt;
			if (cost_opt < cost_min && cost_opt < cost_max)
			{
				this->num_pieces = n_opt;
			}
			else
			{
				if (cost_min < cost_max)
				{
					this->num_pieces = n_min;
				}
				else
				{
					this->num_pieces = n_max;
				}
			}
		}
		else
		{
			if (cost_min < cost_max)
			{
				this->num_pieces = n_min;
			}
			else
			{
				this->num_pieces = n_max;
			}
		}
		std::cout << "        --> Optimum number of pieces: " << this->num_pieces << std::endl;
	}
	else
	{
		std::cout << "        --> It is not worthed to split the wave in pieces!" << std::endl;
		this->num_pieces = 1;
	}
	std::cout << "... done!" << std::endl;
}

void IrregularWave::GetPiecesWaveLengths(void)
{
	std::cout << "--> Getting pieces wave lengths..." << std::endl;

	std::tie(this->lambdas_piece,
			 this->k_piece,
			 this->kx_piece,
			 this->ky_piece,
			 this->kx_1D_piece,
			 this->ky_1D_piece) = this->GetWaveLengths(this->periods_piece, this->headings);

	std::cout << "--> ... pieces wave lengths computed!" << std::endl;
}

void IrregularWave::CutPiecesSpectrumZeros(void)
{
	std::cout << "--> Cropping pieces spectrum to avoid zeros..." << std::endl;

	// Find the maximum amplitude value among the pieces
	arma::mat amplitudes_max = arma::zeros(num_comps_piece, num_headings);
	for (int ii = 0; ii < this->num_pieces; ii++)
	{
		for (int jj = 0; jj < num_comps_piece; jj++)
		{
			for (int kk = 0; kk < num_headings; kk++)
			{
				amplitudes_max(jj, kk) = std::max(amplitudes_max(jj, kk), this->amplitudes_piece(ii)(jj, kk));
			}
		}
	}
	// Compute the sum of the maximum amplitudes
	double amplitudes_max_mean = arma::accu(amplitudes_max) / (num_comps_piece * num_headings_piece);
	// Find the indices of the frequencies and headings with amplitudes greater than factor * amplitudes_max_mean
	arma::vec amplitudes_max_freqs = arma::max(amplitudes_max, 1);
	arma::uvec ind_rows = arma::find(amplitudes_max_freqs >= factor * amplitudes_max_mean);
	arma::vec amplitudes_max_headings = arma::max(amplitudes_max, 0);
	arma::uvec ind_cols = arma::find(amplitudes_max_headings >= factor * amplitudes_max_mean);
	// Crop the spectrum to avoid zeros
	this->num_comps_piece = ind_rows.n_elem;
	this->num_headings_piece = ind_cols.n_elem;
	this->periods_piece = this->periods_piece.elem(ind_rows);
	this->freqs_piece = this->freqs_piece.elem(ind_rows);
	this->ang_freqs_piece = this->ang_freqs_piece.elem(ind_rows);
	this->headings_piece = this->headings.elem(ind_cols);
	this->k_piece = this->k_piece.elem(ind_rows);
	this->kx_piece = this->kx_piece.elem(ind_rows, ind_cols);
	this->ky_piece = this->ky_piece.elem(ind_rows, ind_cols);
	this->kx_1D_piece = this->kx_1D_piece.elem(ind_rows);
	this->ky_1D_piece = this->ky_1D_piece.elem(ind_rows);
	for (int ii = 0; ii < this->num_pieces; ii++)
	{
		this->amplitudes_piece(ii) = this->amplitudes_piece(ii).submat(ind_rows, ind_cols);
		this->phases_piece(ii) = this->phases_piece(ii).submat(ind_rows, ind_cols);
	}

	std::cout << "--> ... pieces spectrum cropped!" << std::endl;
}