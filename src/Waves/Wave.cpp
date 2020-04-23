
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

void RegularWave::GetWaveSpectrum(void)
{
	nPeriods = 1;
	nHeadings = 1;
	df = 1/t_sim; dw = 2*pi*df;
	periods = arma::ones(1,1)*period;
	ang_freqs = arma::ones(1,1)*(pi/period);
	headings = arma::ones(1,1)*heading;
	amplitudes = arma::ones(1,1)*height;
	phases = arma::zeros(1,1);
}


void IrregularWave::GetWaveSpectrum(void)
{
	if (specType_flag == 1)
	{
		std::cout << "--> Computing Wave Spectrum (JONSWAP)" << std::endl;
		t_sim = std::max(t_sim,120.0);
		df = 1.0/t_sim; dw = 2*pi*df;
		freqs = arma::regspace(0,df,1.0/(2.0*dt));
		periods = 1.0/freqs; nPeriods = periods.n_elem;
		GetJonswapSpectrum(); S_w(0,0) = 0.0;

		if (s>1)
		{
			headings = arma::regspace(heading-90,dtheta,heading+90);
			nHeadings = headings.n_elem;
			GetSpreadingFunction();
		}
		else
		{
			headings = arma::ones(1,1)*heading;
			nHeadings = 1;
			g_theta = arma::ones(1,1);
		}

		amplitudes = sqrt(2.0*S_w*g_theta.t()*df*dtheta*pi/180.0);

		int flag = 1;
		int ii = 1;
		int nIterMax = 20;
		while ((flag==1)&&(ii<=nIterMax))
		{
			phases = amplitudes.randn()*pi;
			flag = CheckPhases(); ii++;
		}

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
    for(int ii=0; ii<nPeriods; ii++)
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

    arma::mat theta = arma::regspace(-90,dtheta,+90)*pi/180.0;    
    arma::mat G = arma::pow(arma::cos(theta),s);
    double inte = arma::as_scalar(arma::trapz(theta,G,1));
    g_theta = G/inte;
    inte = arma::as_scalar(arma::trapz(theta,g_theta,1));
    if (abs(inte-1)>1e-6)
    {
        g_theta = g_theta/inte;
    }
}

int IrregularWave::CheckPhases(void)
{
	int flag = 1;

	arma::mat t = arma::regspace(0,dt,t_sim); int M = t.n_elem;

	arma::mat Cm = arma::sum(amplitudes%arma::cos(phases),1);
	arma::mat Sm = arma::sum(amplitudes%arma::sin(phases),1);
	arma::mat Am = arma::sqrt(arma::pow(Cm,2)+arma::pow(Sm,2));
	arma::mat PHIm = arma::atan2(Sm,Cm);
	arma::cx_mat iPHIm(arma::zeros(size(PHIm)),PHIm);
	arma::cx_mat Y1 = (0.5*M*Am) % arma::exp(-iPHIm); Y1(0,0) = 0.0;
	arma::cx_mat Y;
	if (M % 2 == 0)
	{
		Y = arma::join_vert(Y1,arma::flipud(arma::conj(Y1.rows(1,nPeriods-2))));
	}
	else
	{
		Y = arma::join_vert(Y1,arma::flipud(arma::conj(Y1.rows(1,nPeriods-1))));
	}

	arma::cx_mat eta_cx = arma::ifft(Y);

	double imag = arma::as_scalar(arma::sum(arma::abs(arma::imag(eta_cx)),0));
	std::cout << "imag = " << imag << std::endl;
	if (imag>1e-13*t.n_elem)
	{
		std::stringstream ss;
	    ss << "Something went wrong with the ifft. \n";
	    throw ValueError(ss.str());
	}
	arma::mat eta = arma::real(eta_cx);

	std::cout << "eta = "<< std::endl << eta << std::endl;

	arma::mat T, H;
	std::tie(T,H) = upcrossing(t,eta);

	double Tmean = arma::as_scalar(arma::mean(T));
	double Hmax = H.max();
	arma::mat Haux = arma::sort(H,"descend"); int nH = std::round(H.n_elem/3);
	double Hsig = arma::as_scalar(arma::mean(H.rows(0,nH)));

	double TmeanC = period/(sqrt(2)*pow(gamma,-0.082));
	double HsigC = height;
	double HmaxC = HsigC*1.8;

	if ((Tmean<=TmeanC*(1+rel_tol))&&(Tmean>=TmeanC*(1-rel_tol))&&
		(Hsig <= HsigC*(1+rel_tol))&&(Hsig >= HsigC*(1-rel_tol))&&
		(Hmax <= HmaxC*(1+rel_tol))&&(Hmax >= HmaxC*(1-rel_tol)))
	{
		flag = 0;
	}

	std::cout << "Tmean = " << Tmean << ";  TmeanC = " << TmeanC << std::endl;
	std::cout << "Hsig = " << Hsig << ";  HsigC = " << HsigC << std::endl;
	std::cout << "Hmax = " << Hmax << ";  HmaxC = " << HmaxC << std::endl;
	std::cout << "flag = " << flag << ";  flag_ini = 1" << std::endl;

	std::stringstream ss;
	ss << "Done with this \n";
	throw ValueError(ss.str());

	return flag;
}

void IrregularWave::ReadWaveSpectrumHDF5(void)
{
	std::cout << "--> Reading Wave Spectrum (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWaveSpectrumHDF5 in class IrregularWave not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "----> Wave Spectrum Read" << std::endl;
}