
#include <armadillo>
#include <tuple>
#include "Wave.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"

Wave::Wave(double H, double T, double D)
{
    height = H;
    period = T;
    heading = D;
}

void Wave::CheckBreakingWave(void)
{
    std::cout << "--> Checking wave breaking limits." << std::endl;

    double lambda = solve_lambda(period);
	double hL = waterDepth/lambda;

    double hmax;

    if (hL>0.108)
    {
    	hmax = 0.142*tanh(2.0*pi*hL)*lambda;
    }
    else if (hL<=0.108)
    {
    	hmax = 0.78*waterDepth;
    }

    if (height>hmax)
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
	lambdas = arma::zeros(size(periods)); double T;
	for(int ii=0; ii<num_comps; ii++)
	{
		T = arma::as_scalar(periods(ii,0));
		lambdas(ii,0) = solve_lambda(T);
	}
	k = 2.0*pi/lambdas;
	arma::mat cos_theta = arma::cos(pi/180.0*headings.t());
	arma::mat sin_theta = arma::sin(pi/180.0*headings.t());
	kx = k*cos_theta; ky = k*sin_theta;
    kx_1D = k*cos(heading);
    ky_1D = k*sin(heading);
}


double Wave::solve_lambda(double T)
{
	double lambda = gravity*T*T/(2.0*pi);
	double f = f_lambda(lambda,T); double df = 0;
	int ii = 0; int nIterMax = 100; double atol = 1e-6;
	while ((abs(f)>atol)&&(ii<nIterMax))
	{
		df = df_lambda(lambda,T);
		lambda = lambda - f/df;
		f = f_lambda(lambda,T);
		ii++;
	}

	if (abs(f)>atol)
	{
		std::stringstream ss;
	    ss << "Convergence failed for computing wave length. \n";
	    throw ValueError(ss.str());
	}
	return lambda;
}

double Wave::f_lambda(double lambda, double T)
{
	return gravity*T*T/(2.0*pi)*tanh(2.0*pi*waterDepth/lambda)-lambda;
}

double Wave::df_lambda(double lambda, double T)
{
	return -gravity*T*T*waterDepth/(lambda*lambda*cosh(2.0*pi*waterDepth/lambda))-1.0;
}

void RegularWave::GetWaveSpectrum(void)
{
	num_comps = 1;
	num_headings = 1;
	df = 1.0/simulationTime; dw = 2.0*pi*df;
	periods = arma::ones(1,1)*period;
	ang_freqs = arma::ones(1,1)*(2.0*pi/period);
	headings = arma::ones(1,1)*heading; headings = mod(headings,2.0*pi);
	amplitudes = arma::ones(1,1)*height;
	phases = arma::zeros(1,1);
	GetWaveLengths();
}


void IrregularWave::GetWaveSpectrum(void)
{
	if (specType_flag == 1)
	{
		std::cout << "--> Computing Wave Spectrum (JONSWAP)" << std::endl;

		simulationTime = std::max(simulationTime,400.0);
		num_points = round(simulationTime/dt + 1);
		num_comps = floor(num_points/2.0)+1;
		freqs = arma::linspace(0,1.0/(2.0*dt),num_comps);
		ang_freqs = 2.0*pi*freqs;
		periods = 1.0/freqs;
		df = arma::as_scalar(freqs(1,0)); dw = 2.0*pi*df;

		GetJonswapSpectrum(); S_w(0,0) = 0.0;

		if (s>1.0)
		{
			headings = arma::regspace(heading-90,dtheta,heading+90);
			num_headings = headings.n_elem;
			GetSpreadingFunction();
		}
		else
		{
			headings = arma::ones(1,1)*heading;
			num_headings = 1;
			g_theta = arma::ones(1,1);
			dtheta = 180.0/pi;
		}
		headings = mod(headings,2.0*pi);

		amplitudes = arma::sqrt(2.0*(S_w*g_theta.t())*df*(dtheta*pi/180.0));

		// Repeat the computetion without directional spreading for the QTFs
		headings_1D = arma::ones(1,1)*heading;
	    amplitudes_1D = arma::sqrt(2.0*S_w*df);

		int flag = 1;
		int ii = 1;
		int nIterMax = 20;
		while ((flag==1)&&(ii<=nIterMax))
		{
			phases = arma::randn(size(amplitudes))*pi;
			flag = CheckPhases(); ii++;
		}
		if (flag==1)
		{
			std::stringstream ss;
		    ss << "It was not possible to find a realistic wave. \n";
		    throw ValueError(ss.str());
		}

		GetWaveLengths();
		CutSpectrumZeros();

		std::cout << "----> Wave Spectrum Computed" << std::endl;
	}
	else if (specType_flag == 2)
	{
		ReadWaveSpectrumHDF5();
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

    // Function jonswap: It returns the jonswap frequency spectra according 
    // to the IEC 61400-3 standard. It does it for certain frequencies.
    //
    //
    //  Alvaro Rodriguez Luis
    //  Feb. 2020
    //  IH Cantabria

    double c1 = 0.3125;
    double c2 = 0.287;
    
    double fp = 1.0/period;

    arma::mat sigma = arma::zeros(size(freqs));
    for(int ii=0; ii<num_comps; ii++)
	{
		if (freqs(ii,0)<fp){
			sigma(ii,0) = 0.07;
		} else {
			sigma(ii,0) = 0.09;
		}
	}
    
    double cte = c1*height*height*period*(1.0-c2*log(gamma));
    arma::mat alpha = arma::exp(-0.5*(arma::pow((freqs-fp)/(freqs%sigma),2.0)));
    arma::mat temp1 = arma::pow(freqs/fp,-5.0);
    arma::mat temp2 = arma::exp(-1.25*arma::pow(freqs/fp,-4.0));
    arma::mat temp3 = arma::exp(log(gamma)*alpha);
    S_w = cte * temp1 % temp2 % temp3;
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

    arma::mat theta = arma::regspace(-90,dtheta,90)*pi/180.0;    
    arma::mat G = arma::pow(arma::cos(theta),s);
    double inte = arma::as_scalar(arma::trapz(theta,G,0));
    g_theta = G/inte;
    inte = arma::as_scalar(arma::trapz(theta,g_theta,0));
    if (abs(inte-1)>1e-6)
    {
        g_theta = g_theta/inte;
    }
}

int IrregularWave::CheckPhases(void)
{
	int flag = 1;

	arma::mat t = arma::linspace(0.0,simulationTime,num_points);

	arma::mat Cm = arma::sum(amplitudes%arma::cos(phases),1);
	arma::mat Sm = arma::sum(amplitudes%arma::sin(phases),1);
	arma::mat Am = arma::sqrt(arma::pow(Cm,2)+arma::pow(Sm,2));
	arma::mat PHIm = arma::atan2(Sm,Cm);
	arma::cx_mat iPHIm(arma::zeros(size(PHIm)),PHIm);
	arma::cx_mat Y1 = (0.5*num_points*Am) % arma::exp(-iPHIm); Y1(0,0) = 0.0;
	arma::cx_mat Y;
	if (num_points % 2 == 0)
	{
		Y = arma::join_vert(Y1,arma::flipud(arma::conj(Y1.rows(1,num_comps-2))));
	}
	else
	{
		Y = arma::join_vert(Y1,arma::flipud(arma::conj(Y1.rows(1,num_comps-1))));
	}

	arma::cx_mat eta_cx = arma::ifft(Y);

	double imag = arma::as_scalar(arma::sum(arma::abs(arma::imag(eta_cx)),0));
	if (imag>1e-13*num_points)
	{
		std::stringstream ss;
	    ss << "Something went wrong with the ifft. \n";
	    throw ValueError(ss.str());
	}
	arma::mat eta = arma::real(eta_cx);

	arma::mat T, H;
	std::tie(T,H) = upcrossing(t,eta);

	double Tmean = arma::as_scalar(arma::mean(T));
	double Hmax = H.max();
	arma::mat Haux = arma::sort(H,"descend"); 
	int nH = std::round(H.n_elem/3);
	double Hsig = arma::as_scalar(arma::mean(Haux.rows(0,nH)));

	double TmeanC = period/(sqrt(2.0)*pow(gamma,-0.082));
	double HsigC = height;
	double HmaxC = HsigC*1.8;

	if ((Tmean<=TmeanC*(1+rel_tol))&&(Tmean>=TmeanC*(1-rel_tol))&&
		(Hsig <= HsigC*(1+rel_tol))&&(Hsig >= HsigC*(1-rel_tol))&&
		(Hmax <= HmaxC*(1+rel_tol))&&(Hmax >= HmaxC*(1-rel_tol)))
	{
		flag = 0;
		arma::cx_mat spc =  arma::fft(eta);
		phases_1D = arma::atan2(arma::imag(spc),arma::real(spc));
	}

	//std::cout << "Tmean = " << Tmean << ";  Tmean_th = " << TmeanC << std::endl;
	//std::cout << "Hsig = " << Hsig << ";  Hsig_th = " << HsigC << std::endl;
	//std::cout << "Hmax = " << Hmax << ";  Hmax_th = " << HmaxC << std::endl;

	return flag;
}

void IrregularWave::CutSpectrumZeros(void)
{
	// Quitar las componentes frequenciales y direccionales muy proximas a cero, por implementar
	arma::uvec ind_rows = arma::find(S_w    >=factor*arma::as_scalar(arma::mean(S_w    )));
	arma::uvec ind_cols = arma::find(g_theta>=factor*arma::as_scalar(arma::mean(g_theta)));
	arma::uvec ind_0 = arma::zeros<arma::uvec>(1);

	num_comps = ind_rows.n_elem;
	num_headings = ind_cols.n_elem;

	periods = periods.submat(ind_rows,ind_0);
	freqs = freqs.submat(ind_rows,ind_0);
	ang_freqs = ang_freqs.submat(ind_rows,ind_0);
	headings = headings.submat(ind_cols,ind_0);
	amplitudes = amplitudes.submat(ind_rows,ind_cols);
	phases = phases.submat(ind_rows,ind_cols);
	lambdas = lambdas.submat(ind_rows,ind_0);
	S_w = S_w.submat(ind_rows,ind_0);
	g_theta = g_theta.submat(ind_cols,ind_0);
	k = k.submat(ind_rows,ind_0);
	kx = kx.submat(ind_rows,ind_cols);
	ky = ky.submat(ind_rows,ind_cols);

	amplitudes_1D = amplitudes_1D.submat(ind_rows,ind_0);
    phases_1D = phases_1D.submat(ind_rows,ind_0);
    kx_1D = kx_1D.submat(ind_rows,ind_0);
    ky_1D = ky_1D.submat(ind_rows,ind_0);

}

void IrregularWave::ReadWaveSpectrumHDF5(void)
{
	std::cout << "--> Reading Wave Spectrum (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWaveSpectrumHDF5 in class IrregularWave not implemented yet. \n";
    throw NotImplementedError(ss.str());
    std::cout << "----> Wave Spectrum Read" << std::endl;
}