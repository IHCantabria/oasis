
#ifndef wavedef_hpp__
#define wavedef_hpp__
#include <armadillo>
#include <string>
#include <cstdio>

// Class predefinition in order to avoid class cross-linking
class Simulation;

class Wave
{
public:
    // Declare class variables
	Simulation* pSim; // Pointer to simulation instance
    double pi = arma::datum::pi;
    double height;
    double period;
    double heading;
    
    double simulationTime;
    double df;
    double dw;
    double dt;
    double dtheta;

    int num_comps;
    int num_points;
    arma::mat periods;
    arma::mat freqs;
    arma::mat ang_freqs;
    int num_headings;
    arma::mat headings;
    arma::mat S_w;
    arma::mat g_theta;

    arma::mat amplitudes;
    arma::mat phases;

    arma::mat t_FS, eta_FS;

    double gravity;
    double waterDepth;
    double lambda_peak;
    arma::mat lambdas;
    arma::mat k;
    arma::mat kx;
    arma::mat ky;

    arma::mat headings_1D;
    arma::mat amplitudes_1D;
    arma::mat phases_1D;
    arma::mat kx_1D;
    arma::mat ky_1D;

    // Variables for irregular waves
    int specType_flag;
	double gamma;
	double s;
	double rel_tol;
    double factor;
	std::string waveDatabaseName;
    std::string file_path;

    FILE* pfile_SPEC;
    FILE* pfile_TIME;

    // Declare class constructors
    Wave(Simulation* pSimInc, double H, double T, double D);

    // Methods

    virtual void GetWaveSpectrum(void) = 0;
    void CheckBreakingWave(void);
    void GetWaveLengths(void);
    double solve_lambda(double T);
    double f_lambda(double lambda, double T);
    double df_lambda(double lambda, double T);
    void GetFreeSurface(void);
    arma::vec GetFreeSurface(double time, arma::vec x, arma::vec y);
    arma::vec GetPressure(double time, arma::vec x, arma::vec y, arma::vec z);
    void WriteOut(std::string path);
};

class RegularWave: public Wave
{
public:
	RegularWave(Simulation* pSimInc, double H, double T, double D): Wave(pSimInc,H,T,D){};
	void GetWaveSpectrum(void);
};


class IrregularWave: public Wave
{
public:
	IrregularWave(Simulation* pSimInc, double H, double T, double D): Wave(pSimInc,H,T,D){};
	void GetWaveSpectrum(void);
	void ReadWaveSpectrumHDF5(void);
    void ReadWaveSpectrumASCII(void);
	void GetJonswapSpectrum(void);
	void GetSpreadingFunction(void);
    void CutSpectrumZeros(void);
	int CheckPhases(void);
};

#endif // wavedef_hpp__