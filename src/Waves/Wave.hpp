
#ifndef wavedef_hpp__
#define wavedef_hpp__
#include <armadillo>
#include <string>
#include <cstdio>

class Wave
{
public:
    // Declare class variables
    double pi = arma::datum::pi;
    double height;
    double period;
    double heading;
    
    double t_sim;
    double df;
    double dw;
    double dt;
    double dtheta;

    int nPeriods;
    arma::mat periods;
    arma::mat freqs;
    arma::mat ang_freqs;
    int nHeadings;
    arma::mat headings;
    arma::mat S_w;
    arma::mat g_theta;

    arma::mat amplitudes;
    arma::mat phases;

    // Variables for irregular waves
    int specType_flag;
	double gamma;
	double s;
	double rel_tol;
	std::string waveDatabaseName;

    // Declare class constructors
    Wave(double H, double T, double D);

    // Methods
    virtual void GetWaveSpectrum(void) = 0;
};

class RegularWave: public Wave
{
public:
	RegularWave(double H, double T, double D): Wave(H,T,D){};
	void GetWaveSpectrum(void);
};


class IrregularWave: public Wave
{
public:
	IrregularWave(double H, double T, double D): Wave(H,T,D){};
	void GetWaveSpectrum(void);
	void ReadWaveSpectrumHDF5(void);
	void GetJonswapSpectrum(void);
	void GetSpreadingFunction(void);
	int CheckPhases(void);
};

#endif // wavedef_hpp__