
#include <armadillo>
#include <tuple>
#include "Wave.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"
#include <math.h>

Wave::Wave(Simulation* pSimInc,double H, double T, double D)
{
    height = H;
    period = T;
    heading = D*pi/180.0;
	pSim = pSimInc;
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
	std::cout << "--> Getting Wave Lengths" << std::endl;

	lambdas = arma::zeros(size(periods)); double T;
	k = lambdas;
	for(int ii=0; ii<num_comps; ii++)
	{
		T = arma::as_scalar(periods(ii,0));
		if (T == 0){
			lambdas(ii,0) = 0.0;
			k(ii,0) = 0.0;
		} else if (std::isinf(T)) {
			lambdas(ii,0) = 0.0;
			k(ii,0) = 0.0;
		} else {
			lambdas(ii,0) = solve_lambda(T);
			k(ii,0) = 2.0*pi/lambdas(ii,0);
		}
	}
	arma::mat cos_theta = arma::cos(headings.t());
	arma::mat sin_theta = arma::sin(headings.t());
	kx = k*cos_theta; ky = k*sin_theta;
    kx_1D = k*cos(heading);
    ky_1D = k*sin(heading);
    headings = wrapToPi(headings);
	lambda_peak = solve_lambda(period);
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
		arma::uvec ind = arma::find(T==periods);
		double H  = arma::accu(amplitudes.rows(ind));
		std::stringstream ss;
	    ss << "Convergence failed for computing wave length for period T = " << T << " s. \n";
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

void Wave::GetFreeSurface(void)
{
	std::cout << "----> Computing Free Surface" << std::endl;
	t_FS = arma::linspace(0.0,simulationTime,num_points);
	arma::mat Cm = arma::sum(amplitudes%arma::cos(phases),1);
	arma::mat Sm = arma::sum(amplitudes%arma::sin(phases),1);
	arma::mat Am = arma::sqrt(arma::pow(Cm,2)+arma::pow(Sm,2));
	arma::mat PHIm = arma::atan2(Sm,Cm);
	arma::cx_mat iPHIm(arma::zeros(size(PHIm)),PHIm);
	arma::cx_mat Y1 = (0.5*num_points*Am) % arma::exp(-iPHIm); Y1(0,0) = 0.0;
	arma::cx_mat Y;
	if (num_points % 2 == 0)
	{
		Y = arma::join_vert(arma::conj(Y1),arma::flipud(Y1.rows(1,num_comps-2)));
	}
	else
	{
		Y = arma::join_vert(arma::conj(Y1),arma::flipud(Y1.rows(1,num_comps-1)));
	}
	arma::cx_mat eta_cx = arma::ifft(Y);
	double imag = arma::as_scalar(arma::sum(arma::abs(arma::imag(eta_cx)),0));
	if (imag>1e-9*num_points)
	{
		std::stringstream ss;
	    ss << "Something went wrong with the ifft. imag = " << imag << " \n";
	    throw ValueError(ss.str());
	}
	eta_FS = arma::real(eta_cx);
	std::cout << "----> Free Surface Computed" << std::endl;
}

arma::vec Wave::GetFreeSurface(double time, arma::vec x, arma::vec y)
{

	arma::vec eta = arma::sum(amplitudes*ones(size(amplitudes.t())) * arma::cos(phases*ones(size(x.t()))
							 + kx*x.t() + ky*y.t() - time*ang_freqs*ones(size(x.t())))).t();

	return eta;
}

arma::vec Wave::GetPressure(double time, arma::vec x, arma::vec y, arma::vec z)
{
	// std::cout << "--> Getting Pressure " << std::endl;

	arma::vec dPhidt, dPhidx, dPhidy, dPhidz, pressure;
	arma::vec A = amplitudes, W = ang_freqs; double H = waterDepth;
	arma::mat onesNp = arma::ones(size(x.t()));
	arma::mat aux = phases*onesNp + kx*x.t() + ky*y.t() - time*ang_freqs*onesNp;

	dPhidt = arma::sum(-((A % arma::pow(W,2))*onesNp) % arma::cos(aux) % arma::cosh(k * (H*onesNp+z.t())) / 
	                    ((k % arma::sinh(H * k))*onesNp)).t();
	dPhidx = arma::sum( ((A % kx % W)*onesNp) % arma::cos(aux) % arma::cosh(k * (H*onesNp+z.t())) / 
						((k % arma::sinh(H * k))*onesNp)).t();
	dPhidy = arma::sum( ((A % ky % W)*onesNp) % arma::cos(aux) % arma::cosh(k * (H*onesNp+z.t())) / 
						((k % arma::sinh(H * k))*onesNp)).t();
	dPhidz = arma::sum( ((A % W)*onesNp % arma::sin(aux)) % arma::sinh(k * (H*onesNp+z.t())) /
			          	(arma::sinh(H * k)*onesNp)).t();

	pressure = -pSim->waterDensity*(pSim->gravity*z + dPhidt + 0.5*(dPhidx%dPhidx + dPhidy%dPhidy + dPhidz%dPhidz));

	return pressure;
}

void Wave::WriteOut(std::string path)
{
	char buffer1[50];
	int nn1 = sprintf(buffer1,"WaveSpectrum.txt");
	std::string file_path1 = JoinPath(path, buffer1);
	pfile_SPEC = fopen (file_path1.c_str(),"ang_freqs");
	if (pfile_SPEC == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn1 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	for(int ii=0;ii<num_comps;ii=ii+1) fprintf(pfile_SPEC, "%f    %f  \n",freqs(ii,0),S_w(ii,0));
	fclose(pfile_SPEC);

	char buffer2[50];	
	int nn2 = sprintf(buffer2,"WaveTimeSeries.txt");	
	std::string file_path2 = JoinPath(path, buffer2);	
	pfile_TIME = fopen (file_path2.c_str(),"ang_freqs");
	if (pfile_TIME == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn2 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}	

	for(int ii=0;ii<num_points;ii=ii+1) fprintf(pfile_TIME, "%f    %f  \n",t_FS(ii,0),eta_FS(ii,0));
	
	fclose(pfile_TIME);
}

void RegularWave::GetWaveSpectrum(void)
{
	num_comps = 1;
	num_headings = 1;
	df = 1.0/simulationTime; dw = 2.0*pi*df;
	periods = arma::ones(1,1)*period;
	freqs = 1.0/periods;
	ang_freqs = arma::ones(1,1)*(2.0*pi/period);
	headings = arma::ones(1,1)*heading; headings = mod(headings,2.0*pi);
	amplitudes = arma::ones(1,1)*height*0.5;
	phases = arma::zeros(1,1);
	headings_1D = headings;
	amplitudes_1D = amplitudes;
	phases_1D = phases;
	GetWaveLengths(); 
}

void IrregularWave::GetWaveSpectrum(void)
{
	dtheta = dtheta*pi/180;
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
			headings = arma::regspace(heading-0.5*pi,dtheta,heading+0.5*pi);
			num_headings = headings.n_elem;
			GetSpreadingFunction();
		}
		else
		{
			headings = arma::ones(1,1)*heading;
			num_headings = 1;
			g_theta = arma::ones(1,1);
			dtheta = 1.0;
		}
		headings = mod(headings,2.0*pi);

		amplitudes = arma::sqrt(2.0*(S_w*g_theta.t())*df*dtheta);

		// Repeat the computetion without directional spreading for the QTFs
		headings_1D = arma::ones(1,1)*heading;
	    amplitudes_1D = arma::sqrt(2.0*S_w*df);

		int flag = 1;
		int ii = 1;
		int nIterMax = 200;
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
		GetFreeSurface();
		CutSpectrumZeros();

		std::cout << "----> Wave Spectrum Computed" << std::endl;
	}
	else if (specType_flag == 2)
	{
		ReadWaveSpectrumHDF5();
	}
	else if (specType_flag == 3)
	{
		ReadWaveSpectrumASCII();
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
		Y = arma::join_vert(arma::conj(Y1),arma::flipud(Y1.rows(1,num_comps-2)));
	}
	else
	{
		Y = arma::join_vert(arma::conj(Y1),arma::flipud(Y1.rows(1,num_comps-1)));
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

void IrregularWave::ReadWaveSpectrumASCII(void)
{
	std::cout << "--> Reading Wave Spectrum (ASCII format)" << std::endl;

	int nn;
	arma::mat time, eta;
	double dtemp;

    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: dataWaves.dat\n    ->Dir: " << file_path << std::endl;
        throw IOError(ss.str());
	}
	
    char bufferLine [1000];

    fscanf(file_pointer, "%i", &nn, bufferLine);
    fscanf(file_pointer, "%[^\n]\n", bufferLine);
    time = arma::zeros(nn,1); eta = arma::zeros(nn,1);

    for (int ii=0; ii<nn; ii++)
	{
		fscanf(file_pointer, "%lf", &dtemp); time(ii,0) = dtemp;
		fscanf(file_pointer, "%lf", &dtemp); eta(ii,0) = dtemp;
		fscanf(file_pointer, "%[^\n]\n", bufferLine);
	}

    // Close file
    fclose(file_pointer);

    simulationTime = time.max()-time(0,0);
    dt = arma::as_scalar(time(1,0)-time(0,0));
    num_points = nn;
    num_comps = floor(num_points/2.0)+1;
    freqs = arma::linspace(0.0,1.0/(dt*2.0),num_comps);
    df = arma::as_scalar(freqs(1,0)-freqs(0,0));
    ang_freqs = 2.0*pi*freqs; dw = 2.0*pi*df;
	periods = 1.0/freqs; dtheta = 1.0;

    arma::cx_mat yf = arma::fft(eta)/num_points;
    phases = arma::atan2(arma::imag(yf.rows(0,num_comps-1)), arma::real(yf.rows(0,num_comps-1)));
    arma::mat psd = arma::pow(arma::abs(yf.rows(0,num_comps-1)),2)/df;
    psd.rows(1,num_comps-1) = 2.0*psd.rows(1,num_comps-1); psd(0,0) = 0.0;
    amplitudes = arma::sqrt(2.0*psd*df);

    num_headings = 1;
	g_theta = arma::ones(1,1);
	dtheta = 1.0;
	headings = arma::ones(1,1)*heading;
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

    std::cout << "----> Wave Spectrum Read" << std::endl;
}