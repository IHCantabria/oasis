
#include <armadillo>
#include <tuple>
#include <math.h>

#include "Wave.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"

Wave::Wave(Simulation *pSimInc, double H, double T, double D)
{
	height = H;
	period = T;
	heading = D * arma::datum::pi / 180.0;
	heading = std::fmod(heading + arma::datum::pi, 2.0 * arma::datum::pi);
	if (heading < 0)
	{
		heading += 2.0 * arma::datum::pi;
	}
	heading -= arma::datum::pi;
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

void Wave::SetZeroHeight(void)
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
	amplitudes = arma::zeros(1, 1);
	phases = arma::zeros(1, 1);
	headings_1D = headings;
	amplitudes_1D = amplitudes;
	phases_1D = phases;
	GetWaveLengths();
	SetSinglePiece();
}

void Wave::CheckBreakingWave(void)
{
	std::cout << "    --> Checking wave breaking limits..." << std::endl;

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
		std::cout << "    --> ... wave within wave breaking limits!" << std::endl;
	}
}

void Wave::GetWaveLengths(void)
{
	std::cout << "    --> Computing wave lengths..." << std::endl;

	std::tie(this->lambdas,
			 this->k,
			 this->kx,
			 this->ky,
			 this->kx_1D,
			 this->ky_1D) = this->GetWaveLengths(this->periods, this->headings);

	std::cout << "    --> ... done!" << std::endl;
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
	int num_comps = periods.n_elem;
	int num_headings = headings.n_elem;
	arma::vec lambdas = arma::zeros(size(periods));
	arma::vec k = arma::zeros(size(periods));
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
	std::cout << "    --> Computing Free Surface" << std::endl;
	this->t_FS = arma::linspace(0.0, this->simulationTime, this->num_points);
	this->eta_FS = GetFreeSurface(this->amplitudes, this->phases, this->num_points);
	std::cout << "    --> Free Surface Computed" << std::endl;
}

arma::vec Wave::GetFreeSurface(arma::mat amplitudes, arma::mat phases, int num_points)
{
	arma::vec Cm = arma::sum(amplitudes % arma::cos(phases), 1);
	arma::vec Sm = arma::sum(amplitudes % arma::sin(phases), 1);
	arma::vec Am = arma::sqrt(arma::pow(Cm, 2.0) + arma::pow(Sm, 2.0));
	arma::vec PHIm = arma::atan2(Sm, Cm);
	arma::cx_vec iPHIm(arma::zeros(size(PHIm)), PHIm);
	arma::cx_vec Y1 = (0.5 * num_points * Am) % arma::exp(-iPHIm);
	Y1(0) = 0.0;
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
	// std::cout << "    --> Computing Pressure..." << std::endl;
	// std::cout << "      --> num_points = " << x.n_elem << std::endl;
	// std::cout << "      --> num_freqs = " << num_comps << std::endl;
	// std::cout << "      --> num_headings = " << num_headings << std::endl;

	// arma::mat onesNp = arma::ones(size(x.t()));
	// arma::mat eta_term = (A * onesNp) % arma::cos(KX * x.t() + KY * y.t() - time * W * onesNp + P * onesNp);
	// arma::mat coef = arma::cosh(K * z.t()) + arma::sinh(K * z.t()) % arma::tanh(K * waterDepth * onesNp);
	// arma::vec dynamic_term = ((eta_term % coef).t() * arma::ones(size(A)));
	// arma::vec eta = (eta_term.t() * arma::ones(size(A)));
	arma::vec time_arg = P - time * W;
	arma::mat eta_term = KX * x.t();
	eta_term += KY * y.t();
	eta_term.each_col() += time_arg;
	eta_term = arma::cos(eta_term);
	arma::vec eta = (A.t() * eta_term).t();
	// std::cout << "      --> Computed free surface!" << std::endl;

	arma::vec tkw = arma::tanh(K * waterDepth);
	arma::mat tmp = K * z.t();
	arma::mat coef = arma::cosh(tmp);
	tmp = arma::sinh(tmp);
	tmp.each_col() %= tkw;
	coef += tmp;
	tmp.reset();
	eta_term.each_col() %= A;
	arma::vec dynamic_term = arma::sum(eta_term % coef, 0).t();
	// std::cout << "      --> Computed dynamic term!" << std::endl;

	arma::uvec ind_z_pos = arma::find(z > 0);
	if (ind_z_pos.n_elem > 0)
	{
		dynamic_term(ind_z_pos) = eta(ind_z_pos);
	}
	// std::cout << "      --> Corrected points over MWL!" << std::endl;

	// arma::vec dynamic_term = eta;

	arma::uvec ind_z_eme = arma::find(z > eta);
	if (ind_z_eme.n_elem > 0)
	{
		dynamic_term(ind_z_eme) = z(ind_z_eme);
	}
	// std::cout << "      --> Corrected points over free surface!" << std::endl;

	arma::vec pressure = pSim->waterDensity * pSim->gravity * (dynamic_term - z);
	// std::cout << "      --> Pressure computed!" << std::endl;

	// std::cout << "    --> ... done!" << std::endl;

	return pressure;
}

void Wave::VectoriseComponentsMatrices(void)
{
	A = arma::vectorise(amplitudes);
	P = arma::vectorise(phases);
	KX = arma::vectorise(kx);
	KY = arma::vectorise(ky);
	K = arma::repmat(k, num_headings, 1);
	W = arma::repmat(ang_freqs, num_headings, 1);
	// std::cout << "      --> Vectorized wave terms!" << std::endl;
}

void Wave::SetSinglePiece(void)
{
	num_pieces = 1;

	time_ref = arma::vec(1);
	time_ini = arma::vec(1);
	time_end = arma::vec(1);
	amplitudes_piece = arma::field<arma::mat>(1);
	phases_piece = arma::field<arma::mat>(1);
	amplitudes_1D_piece = arma::field<arma::vec>(1);
	phases_1D_piece = arma::field<arma::vec>(1);

	time_ref(0) = 0.0;
	time_ini(0) = 0.0;
	time_end(0) = simulationTime;
	time_gap = 1.0;

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
	headings_1D_piece = headings_1D;
	amplitudes_1D_piece(0) = amplitudes_1D;
	phases_1D_piece(0) = phases_1D;
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
		fprintf(pfile_SPEC, "%f    %f  \n", freqs(ii, 0), spectral_density(ii, 0));
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
		fprintf(pfile_TIME, "%f    %f  \n", t_FS(ii), eta_FS(ii));

	fclose(pfile_TIME);

	std::string filename = JoinPath(path, "WavePhases.txt");
	this->original_phases.save(filename, arma::arma_ascii);
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
	SetSinglePiece();
	VectoriseComponentsMatrices();
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
		this->df_max = spectrum_width / 10.0;
		// Update the wave time if necessary
		this->simulationTime = std::max(simulationTime, std::ceil(1 / this->df_max));

		// Compute the number of pieces and the time per piece
		if (piecewise_flag == 1)
		{
			GetPiecesNumber();
		}
		else
		{
			num_pieces = 1;
			time_piece = simulationTime;
		}

		// Compute the wave for a single piece
		// Compute the time vector
		if (this->num_pieces > 1)
		{
			t_FS = arma::regspace(0, this->dt, simulationTime + 2 * this->time_piece) - this->time_piece;
			std::cout << "    --> Computing Wave Spectrum (JONSWAP) in pieces" << std::endl;
			std::cout << "    --> Number of pieces: " << this->num_pieces << std::endl;
		}
		else
		{
			t_FS = arma::regspace(0, this->dt, simulationTime);
			std::cout << "    --> Computing Wave Spectrum (JONSWAP) in a single piece" << std::endl;
		}
		// Compute the frequency spectrum
		num_points = t_FS.n_elem;
		if (num_points % 2 == 0)
		{
			t_FS = arma::join_vert(t_FS, arma::ones<arma::vec>(1) * (t_FS(num_points - 1) + this->dt));
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
			directional_spreading = arma::ones<arma::vec>(1);
			dtheta = 1.0;
		}
		headings = mod(headings, 2.0 * arma::datum::pi);
		// Get the wave amplitudes combining the frequency and directional spectra
		this->amplitudes = arma::zeros(num_comps, num_headings);
		for (int ii = 0; ii < num_comps; ii++)
		{
			for (int jj = 0; jj < num_headings; jj++)
			{
				this->amplitudes(ii, jj) = std::sqrt(2.0 * spectral_density(ii) * df * directional_spreading(jj) * dtheta);
			}
		}
		// Get the wave checked free surface and phases
		GetCheckedFreeSurface();
		// TODO: Review how to define the 1D wave
		// Get the 1D amplitudes and headings
		this->headings_1D = arma::ones<arma::vec>(1) * this->heading;
		this->amplitudes_1D = arma::sqrt(2.0 * spectral_density * df);
		// Get the 1D phases
		arma::cx_vec spc = arma::fft(this->eta_FS);
		this->phases_1D = arma::atan2(arma::imag(spc), arma::real(spc));
		// Get the wave lengths
		GetWaveLengths();

		std::cout << "    --> Wave Spectrum Computed" << std::endl;
	}
	else if (specType_flag == 2)
	{
		// Read the wave time series and preprocess it
		ReadWaveTimeSeries();
		// Compute the maximum frequency step
		this->df_max = arma::min(arma::diff(this->freqs));
		// Update the wave time if necessary
		this->simulationTime = std::max(simulationTime, std::ceil(1 / this->df_max));
		// Compute the number of pieces and the time per piece
		if (piecewise_flag == 1)
		{
			GetPiecesNumber();
		}
		else
		{
			num_pieces = 1;
			time_piece = simulationTime;
		}
		// CRITICAL-TODO: Review the computation of the time vector for pieces
		if (this->num_pieces > 1)
		{
			std::stringstream ss;
			ss << "Custom wave not implemented with several pieces. \n";
			throw NotImplementedError(ss.str());
		}
	}
	else if (specType_flag == 3)
	{
		ReadWavePSD();
		// TODO: Implement the reading of the wave PSD
	}
	else
	{
		// TODO: Implement the reading of the wave PSD and time series with multiple headings
		std::stringstream ss;
		ss << "Error while parsing file: dataWaves.dat; Unexpected spectrum type. \n";
		throw ValueError(ss.str());
	}

	// Compute the time intervals for each piece
	GetPiecesTimeIntervals();

	// Compute the wave frequency spectrum for each piece
	if (this->num_pieces > 1)
	{
		GetPiecesSpectra();
	}
	else
	{
		SetSinglePiece();
	}

	// Cut the zeros from the wave spectra
	if (piecewise_flag == 1)
	{
		CutPiecesSpectrumZeros();
	}
	CutSpectrumZeros();
	VectoriseComponentsMatrices();
}

void IrregularWave::GetJonswapSpectrum(void)
{
	this->spectral_density = GetJonswapSpectrum(this->freqs, this->height, this->period, this->gamma);
	zero_order_moment = arma::as_scalar(arma::trapz(this->freqs, this->spectral_density));
	first_order_moment = arma::as_scalar(arma::trapz(this->freqs, this->freqs % this->spectral_density));
	second_order_moment = arma::as_scalar(arma::trapz(this->freqs, arma::pow(this->freqs, 2.0) % this->spectral_density));
}

arma::vec IrregularWave::GetJonswapSpectrum(arma::vec freqs, double height, double period, double gamma)
{
	// GetJonswapSpectrum returns the jonswap frequency spectra according
	// to the IEC 61400-3 standard. It does it for certain frequencies.

	double fp = 1.0 / period;
	arma::vec S_PM = 0.3125 * pow(height, 2.0) * pow(fp, 4.0) * arma::pow(freqs, -5.0) % arma::exp(-1.25 * arma::pow(fp / freqs, 4.0));
	if (freqs(0) == 0)
	{
		S_PM(0) = 0.0;
	}

	arma::vec gamma_alpha = arma::zeros(size(freqs));
	double sigma, alpha, tmp_freq;
	for (int ii = 0; ii < freqs.n_elem; ii++)
	{
		tmp_freq = arma::as_scalar(freqs(ii));
		if (tmp_freq <= fp)
		{
			sigma = 0.07;
		}
		else
		{
			sigma = 0.09;
		}
		alpha = exp(-0.5 * pow((tmp_freq - fp) / (sigma * fp), 2));
		gamma_alpha(ii) = pow(gamma, alpha);
	}

	double numerator = arma::as_scalar(arma::trapz(freqs, S_PM));
	double denominator = arma::as_scalar(arma::trapz(freqs, S_PM % gamma_alpha));
	double normalization_factor = numerator / denominator;

	arma::vec S_JS = normalization_factor * S_PM % gamma_alpha;
	if (freqs(0) == 0)
	{
		S_JS(0) = 0.0;
	}

	return S_JS;
}

void IrregularWave::GetSpreadingFunction(void)
{
	// GetSpreadingFunction returns the spreading function cos(dir/2)^2s.
	// Based in Goda (2000), p. 32.
	arma::vec theta = arma::regspace(-90, dtheta, 90) * arma::datum::pi / 180.0;
	arma::vec dir_spr = arma::pow(arma::cos(0.5 * theta), 2.0 * s);
	double normalization_coeff = arma::as_scalar(arma::trapz(theta, dir_spr, 0));
	arma::vec normalized_dir_spr = dir_spr / normalization_coeff;
	// Interpolate the directional spreading to the wave headings
	arma::interp1(wrapToPi(theta - heading), normalized_dir_spr, wrapToPi(headings), directional_spreading, "linear", 0.0);
}

void IrregularWave::GetCheckedFreeSurface(void)
{
	int flag = 0; // Status flag [0: something is wrong; 1: everything is ok]

	// Randomly generate or read the wave phases
	if (readPhases_flag == 0)
	{
		int ii = 1;
		int nIterMax = 200;
		std::cout << "    --> Generating random phases..." << std::endl;
		while ((flag == 0) && (ii <= nIterMax))
		{
			this->phases = arma::randn(size(this->amplitudes)) * arma::datum::pi;
			this->eta_FS = GetFreeSurface(this->amplitudes, this->phases, this->num_points);
			flag = CheckPhases(this->t_FS, this->eta_FS);
			ii++;
			if (flag == 0)
			{
				std::cout << "    --> Iteration: " << ii << " / " << nIterMax << std::endl;
			}
		}
		if (flag == 0)
		{
			std::stringstream ss;
			ss << "It was not possible to find a realistic wave. \n";
			throw ValueError(ss.str());
		}
		else
		{
			std::cout << "    --> ... random phases generated!" << std::endl;
			// Save the original phases to write them out without zero-cutting
			original_phases = phases;
		}
	}
	else if (readPhases_flag == 1)
	{
		std::cout << "    --> Reading wave phases..." << std::endl;
		arma::mat phases_tmp;
		phases_tmp.load(this->filePhases_path, arma::arma_ascii);
		int num_comps_tmp = phases_tmp.n_rows;
		int num_headings_tmp = phases_tmp.n_cols;
		if ((num_comps_tmp != num_comps) || (num_headings_tmp != num_headings))
		{
			std::stringstream ss;
			ss << "Error while parsing file:" << wavePhasesFileName << "; Unexpected number of components or headings. \n";
			throw ValueError(ss.str());
		}
		else
		{
			this->phases = phases_tmp;
		}
		this->eta_FS = GetFreeSurface(this->amplitudes, this->phases, this->num_points);
		flag = CheckPhases(this->t_FS, this->eta_FS);
		if (flag == 0)
		{
			std::stringstream ss;
			ss << "It was not possible to find a realistic wave with the phases provided in file. \n";
			throw ValueError(ss.str());
		}
		else
		{
			std::cout << "    --> ... wave phases read!" << std::endl;
		}
	}
	else
	{
		std::stringstream ss;
		ss << "Error while parsing file: dataWaves.dat; Unexpected readPhases_flag. \n";
		throw ValueError(ss.str());
	}
}

int IrregularWave::CheckPhases(arma::vec t, arma::vec eta)
{
	int flag = 0; // Status flag [0: something is wrong; 1: everything is ok]

	// Compute the upcrossing periods and heights
	arma::vec T, H;
	std::tie(T, H) = upcrossing(t, eta - arma::mean(eta));

	// Extract the relevant statistics from the upcrossing periods and heights
	double Tmean = arma::as_scalar(arma::mean(T));
	double Hmax = H.max();
	arma::vec Haux = arma::sort(H, "descend");
	int nH = std::round(H.n_elem / 3);
	double Hsig = arma::as_scalar(arma::mean(Haux.rows(0, nH)));

	// Compute the relevant statistics from the wave spectrum
	double TmeanC = sqrt(zero_order_moment / second_order_moment);
	double HsigC = sqrt(zero_order_moment) * 4.00;
	double HmaxC = HsigC * 1.86;

	// Check if the wave is realistic, use Hmax only if the time series is long enough
	int n = t.n_elem;
	if ((t(n - 1) - t(0)) < 1800.0)
	{
		if ((Tmean <= TmeanC * (1 + rel_tol)) && (Tmean >= TmeanC * (1 - rel_tol)) &&
			(Hsig <= HsigC * (1 + rel_tol)) && (Hsig >= HsigC * (1 - rel_tol)))
		{
			flag = 1;
		}
		else
		{
			std::cout << "      --> Tmean: " << Tmean << " / " << TmeanC << std::endl;
			std::cout << "      --> Hsig: " << Hsig << " / " << HsigC << std::endl;
		}
	}
	else
	{
		if ((Tmean <= TmeanC * (1 + rel_tol)) && (Tmean >= TmeanC * (1 - rel_tol)) &&
			(Hsig <= HsigC * (1 + rel_tol)) && (Hsig >= HsigC * (1 - rel_tol)) &&
			(Hmax <= HmaxC * (1 + rel_tol)) && (Hmax >= HmaxC * (1 - rel_tol)))
		{
			flag = 1;
		}
		else
		{
			std::cout << "      --> Tmean: " << Tmean << " / " << TmeanC << std::endl;
			std::cout << "      --> Hsig: " << Hsig << " / " << HsigC << std::endl;
			std::cout << "      --> Hmax: " << Hmax << " / " << HmaxC << std::endl;
		}
	}

	return flag;
}

void IrregularWave::CutSpectrumZeros(void)
{
	// Crop the spectrum to avoid zeros
	std::cout << "    --> Cropping spectrum to avoid zeros..." << std::endl;

	int number_components_original = this->num_comps * this->num_headings;

	double min_period = 1.0; // Minimum period to consider a wave component [s]
	arma::uvec ind_rows = arma::find((spectral_density >= factor * arma::as_scalar(arma::mean(spectral_density))) && (periods >= min_period));
	arma::uvec ind_cols = arma::find(directional_spreading >= factor * arma::as_scalar(arma::mean(directional_spreading)));

	num_comps = ind_rows.n_elem;
	num_headings = ind_cols.n_elem;

	periods = periods.elem(ind_rows);
	freqs = freqs.elem(ind_rows);
	ang_freqs = ang_freqs.elem(ind_rows);
	headings = headings.elem(ind_cols);
	amplitudes = amplitudes.submat(ind_rows, ind_cols);
	phases = phases.submat(ind_rows, ind_cols);
	lambdas = lambdas.elem(ind_rows);
	spectral_density = spectral_density.elem(ind_rows);
	directional_spreading = directional_spreading.elem(ind_cols);
	k = k.elem(ind_rows);
	kx = kx.submat(ind_rows, ind_cols);
	ky = ky.submat(ind_rows, ind_cols);

	amplitudes_1D = amplitudes_1D.elem(ind_rows);
	phases_1D = phases_1D.elem(ind_rows);
	kx_1D = kx_1D.elem(ind_rows);
	ky_1D = ky_1D.elem(ind_rows);

	int number_components_cropped = this->num_comps * this->num_headings;
	std::cout << "    --> ... pieces spectrum cropped! From " << number_components_original << " to " << number_components_cropped << std::endl;
}

void IrregularWave::ReadWaveTimeSeries(void)
{
	std::cout << "    --> Reading and processing wave time series..." << std::endl;

	// Declare local variables
	int nn;
	arma::mat time, eta;
	double dtemp;
	char bufferLine[1000];

	// Open file
	FILE *file_pointer = fopen(file_path.c_str(), "r");
	if (file_pointer == NULL)
	{
		std::stringstream ss;
		ss << "Not possible to open the file: dataWaves.dat\n    ->Dir: " << file_path << std::endl;
		throw IOError(ss.str());
	}

	// Read the number of points from file
	fscanf(file_pointer, "%i", &nn, bufferLine);
	fscanf(file_pointer, "%[^\n]\n", bufferLine);

	// Initialize the time and eta vectors
	time = arma::zeros(nn);
	eta = arma::zeros(nn);

	// Read the time and eta vectors from file
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

	// Postprocess the time and eta vectors
	simulationTime = time.max() - time(0, 0);
	// double time_ini = arma::as_scalar(time(0, 0)); // TODO: Include wave offset wrt wind
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
	// phases = phases + ang_freqs * time_ini;  // TODO: Include wave offset wrt wind
	arma::vec psd = arma::pow(arma::abs(yf.rows(0, num_comps - 1)), 2) / df;
	psd.rows(1, num_comps - 1) = 2.0 * psd.rows(1, num_comps - 1);
	psd(0, 0) = 0.0;
	amplitudes = arma::sqrt(2.0 * psd * df);

	num_headings = 1;
	directional_spreading = arma::ones<arma::vec>(1);
	dtheta = 1.0;
	headings = arma::ones<arma::vec>(1) * heading;
	headings_1D = headings;
	amplitudes_1D = amplitudes;
	phases_1D = phases;
	spectral_density = psd;

	// TODO: Check how doing this here affects the piecewise computation
	GetWaveLengths();
	GetFreeSurface();

	std::cout << "    --> ...done!" << std::endl;
}

void IrregularWave::ReadWavePSD(void)
{
	std::cout << "    --> Reading wave PSD" << std::endl;
	std::stringstream ss;
	ss << "Method ReadWavePSD in class IrregularWave not implemented yet. \n";
	throw NotImplementedError(ss.str());
	std::cout << "    --> Wave Spectrum Read" << std::endl;
}

void IrregularWave::FindSpectrumWidth(void)
{
	std::cout << "    --> Computing wave spectrum width" << std::endl;
	// Define the search frequency vector
	arma::vec fm = arma::linspace(0.0, 1.0, 100001);
	if (specType_flag == 1)
	{
		// Search of FWHM Jonswap spectrum
		std::cout << "    --> Computing Wave JONSWAP with finer frequency step..." << std::endl;
		arma::vec S = GetJonswapSpectrum(fm, this->height, this->period, this->gamma);
		std::cout << "    --> ...done!" << std::endl;

		std::cout << "    --> Find spectrum peak..." << std::endl;
		double Smax = S.max();
		if (Smax == 0)
		{
			std::stringstream ss;
			ss << "Error while parsing file: dataWaves.dat; Zero spectrum peak. \n";
			throw ValueError(ss.str());
		}
		std::cout << "    --> ...done!" << std::endl;

		std::cout << "    --> Find first spectrum half height crossing..." << std::endl;
		arma::uvec i1 = arma::find(S > (Smax / 2), 1, "first");
		double f1 = arma::as_scalar(fm(i1));
		std::cout << "    --> ...done!" << std::endl;

		std::cout << "    --> Find last spectrum half height crossing..." << std::endl;
		arma::uvec i2 = arma::find(S > (Smax / 2), 1, "last");
		double f2 = arma::as_scalar(fm(i2));
		std::cout << "    --> ...done!" << std::endl;

		this->spectrum_width = f2 - f1;
	}
	else if (specType_flag == 2)
	{
		std::stringstream ss;
		ss << "FindSpectrumWidth is not implemented for custom wave in time domain \n";
		throw ValueError(ss.str());
	}
	else if (specType_flag == 3)
	{
		std::stringstream ss;
		ss << "FindSpectrumWidth is not implemented for custom wave in frequency domain \n";
		throw ValueError(ss.str());
	}
	else
	{
		std::stringstream ss;
		ss << "Error while parsing file: dataWaves.dat; Unexpected spectrum type. \n";
		throw ValueError(ss.str());
	}
	std::cout << "    --> ... done!" << std::endl;
}

void IrregularWave::GetPiecesNumber(void)
{
	std::cout << "    --> Computing number of pieces to divide wave time series..." << std::endl;
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
	// TODO: Review and fix the memory saving constraints
	// tmp_int = std::ceil(-(time_sim * std::sqrt(36 * std::pow(this->time_gap, 2) -
	// 										   36 * this->time_gap * time_sim +
	// 										   std::pow(time_sim, 2)) +
	// 					  (18 * this->time_gap * time_sim -
	// 					   18 * std::pow(this->time_gap, 2) -
	// 					   std::pow(time_sim, 2))) /
	// 					(18 * std::pow(this->time_gap, 2)));
	// if (n_min < tmp_int)
	// {
	// 	std::cout << "    --> Imposing memory saving constraints (n_min)..." << std::endl;
	// 	n_min = tmp_int;
	// 	std::cout << "        --> n_min: " << n_min << std::endl;
	// }
	// tmp_int = std::ceil(-(time_sim * std::sqrt(36 * std::pow(this->time_gap, 2) -
	// 										   36 * this->time_gap * time_sim +
	// 										   std::pow(time_sim, 2)) -
	// 					  (18 * this->time_gap * time_sim -
	// 					   18 * std::pow(this->time_gap, 2) -
	// 					   std::pow(time_sim, 2))) /
	// 					(18 * std::pow(this->time_gap, 2)));
	// if (n_max > tmp_int)
	// {
	// 	std::cout << "    --> Imposing memory saving constraints (n_max)..." << std::endl;
	// 	n_max = tmp_int;
	// 	std::cout << "        --> n_max: " << n_max << std::endl;
	// }

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
	this->time_piece = (this->simulationTime + (this->num_pieces - 1) * this->time_gap) / this->num_pieces;
	std::cout << "    --> ... done!" << std::endl;
}

void IrregularWave::GetPiecesWaveLengths(void)
{
	std::cout << "    --> Getting pieces wave lengths..." << std::endl;

	std::tie(this->lambdas_piece,
			 this->k_piece,
			 this->kx_piece,
			 this->ky_piece,
			 this->kx_1D_piece,
			 this->ky_1D_piece) = this->GetWaveLengths(this->periods_piece, this->headings);

	std::cout << "    --> ... pieces wave lengths computed!" << std::endl;
}

void IrregularWave::CutPiecesSpectrumZeros(void)
{
	std::cout << "    --> Cropping pieces spectrum to avoid zeros..." << std::endl;

	// Find the maximum amplitude value among the pieces
	int number_components_original = this->num_comps_piece * this->num_headings_piece;
	arma::mat amplitudes_max = arma::zeros(num_comps_piece, num_headings_piece);
	for (int ii = 0; ii < this->num_pieces; ii++)
	{
		for (int jj = 0; jj < num_comps_piece; jj++)
		{
			for (int kk = 0; kk < num_headings_piece; kk++)
			{
				amplitudes_max(jj, kk) = std::max(arma::as_scalar(amplitudes_max(jj, kk)), arma::as_scalar(this->amplitudes_piece(ii)(jj, kk)));
			}
		}
	}

	// Compute the sum of the maximum amplitudes
	double amplitudes_max_mean = arma::accu(amplitudes_max) / (num_comps_piece * num_headings_piece);

	// Find the indices of the frequencies with amplitudes greater than factor*amplitudes_max_mean
	arma::vec amplitudes_max_freqs = arma::max(amplitudes_max, 1);
	arma::uvec ind_rows = arma::find(amplitudes_max_freqs >= factor * amplitudes_max_mean);
	arma::uvec ind_cols;

	// If applicable, find the indices of the headings with amplitudes greater than factor*amplitudes_max_mean
	if (num_headings_piece > 1)
	{
		arma::vec amplitudes_max_headings = arma::max(amplitudes_max, 0);
		ind_cols = arma::find(amplitudes_max_headings >= factor * amplitudes_max_mean);
		// Make sure to include the 0 and 360 degrees headings
		ind_cols = arma::join_vert(ind_cols, arma::zeros<arma::uvec>(1));
		ind_cols = arma::join_vert(ind_cols, arma::ones<arma::uvec>(1) * (num_headings_piece - 1));
		ind_cols = arma::unique(ind_cols);
		// Keep also the indices of the zero-amplitude headings that are not surrounded by zero-amplitude headings
		arma::uvec ind_cols_aux = ind_cols;
		for (int ii = 1; ii < ind_cols_aux.n_elem - 1; ii++)
		{
			int jj = ind_cols_aux(ii);
			if ((amplitudes_max_headings(jj - 1) >= factor * amplitudes_max_mean) || (amplitudes_max_headings(jj + 1) >= factor * amplitudes_max_mean))
			{
				ind_cols = arma::join_vert(ind_cols, arma::ones<arma::uvec>(1) * jj);
			}
		}
		// Sort the indices withouth repetitions
		ind_cols = arma::unique(ind_cols);
	}
	else
	{
		ind_cols = arma::zeros<arma::uvec>(1);
	}

	// Crop the spectrum to avoid zeros
	this->num_comps_piece = ind_rows.n_elem;
	this->num_headings_piece = ind_cols.n_elem;
	int number_components_cropped = this->num_comps_piece * this->num_headings_piece;
	std::cout << "        -> cropping periods..." << std::endl;
	this->periods_piece = this->periods_piece.elem(ind_rows);
	std::cout << "        -> cropping freqs..." << std::endl;
	this->freqs_piece = this->freqs_piece.elem(ind_rows);
	std::cout << "        -> cropping ang_freqs..." << std::endl;
	this->ang_freqs_piece = this->ang_freqs_piece.elem(ind_rows);
	std::cout << "        -> cropping headings..." << std::endl;
	this->headings_piece = this->headings.elem(ind_cols);
	std::cout << "        -> cropping wave numbers..." << std::endl;
	this->k_piece = this->k_piece.elem(ind_rows);
	std::cout << "        -> cropping wave numbers 2D (X)..." << std::endl;
	this->kx_piece = this->kx_piece.elem(ind_rows, ind_cols);
	std::cout << "        -> cropping wave numbers 2D (Y)..." << std::endl;
	this->ky_piece = this->ky_piece.elem(ind_rows, ind_cols);
	std::cout << "        -> cropping wave numbers 1D (X)..." << std::endl;
	this->kx_1D_piece = this->kx_1D_piece.elem(ind_rows);
	std::cout << "        -> cropping wave numbers 1D (Y)..." << std::endl;
	this->ky_1D_piece = this->ky_1D_piece.elem(ind_rows);
	std::cout << "        -> cropping amplitudes..." << std::endl;
	for (int ii = 0; ii < this->num_pieces; ii++)
	{
		this->amplitudes_piece(ii) = this->amplitudes_piece(ii).submat(ind_rows, ind_cols);
		this->phases_piece(ii) = this->phases_piece(ii).submat(ind_rows, ind_cols);
		this->amplitudes_1D_piece(ii) = this->amplitudes_1D_piece(ii).elem(ind_rows);
		this->phases_1D_piece(ii) = this->phases_1D_piece(ii).elem(ind_rows);
	}

	std::cout << "    --> ... pieces spectrum cropped! From " << number_components_original << " to " << number_components_cropped << std::endl;
}

void IrregularWave::GetPiecesTimeIntervals(void)
{
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
}

void IrregularWave::GetPiecesSpectra(void)
{
	std::cout << "    --> Computing piecewise spectrums..." << std::endl;

	int number_components_original = this->num_comps * this->num_headings;

	this->amplitudes_piece = arma::field<arma::mat>(this->num_pieces);
	this->phases_piece = arma::field<arma::mat>(this->num_pieces);
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
	// Repeat for 1D wave
	// TODO: Review how to define the 1D wave
	this->amplitudes_1D_piece = arma::field<arma::vec>(this->num_pieces);
	this->phases_1D_piece = arma::field<arma::vec>(this->num_pieces);
	for (int ii = 0; ii < this->num_pieces; ii++)
	{
		arma::vec tmp_t = time_piece_total + arma::as_scalar(this->time_ref(ii));
		arma::vec tmp_eta;
		arma::interp1(t_FS, eta_FS, tmp_t, tmp_eta, "*linear");
		arma::cx_vec yf = arma::fft(tmp_eta) / num_points_piece;
		arma::vec psd = arma::pow(arma::abs(yf.rows(0, num_comps_piece - 1)), 2) / df_piece;
		psd.rows(1, num_comps_piece - 1) = 2.0 * psd.rows(1, num_comps_piece - 1);
		psd(0, 0) = 0.0;
		this->amplitudes_1D_piece(ii) = arma::sqrt(2.0 * psd * df_piece);
		this->phases_1D_piece(ii) = arma::atan2(arma::imag(yf.rows(0, num_comps_piece - 1)),
												arma::real(yf.rows(0, num_comps_piece - 1))) -
									ang_freqs_piece * this->time_ref(ii); // TODO: Check this
	}
	headings_1D_piece = headings_1D;

	int number_components_piecewise = this->num_comps_piece * this->num_headings_piece;

	std::cout << "    --> ...piecewise spectrums computed! From " << number_components_original << " to " << number_components_piecewise << std::endl;
}