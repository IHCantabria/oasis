
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
    Simulation *pSim; // Pointer to simulation instance
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
    arma::vec periods;
    arma::vec freqs;
    arma::vec ang_freqs;
    int num_headings;
    arma::vec headings;
    arma::vec spectral_density;
    arma::vec directional_spreading;

    arma::mat amplitudes;
    arma::mat phases;
    arma::mat original_phases;

    arma::vec t_FS, eta_FS;

    double gravity;
    double waterDepth;
    double lambda_peak;
    arma::vec lambdas;
    arma::vec k;
    arma::mat kx;
    arma::mat ky;

    arma::vec headings_1D;
    arma::vec amplitudes_1D;
    arma::vec phases_1D;
    arma::vec kx_1D;
    arma::vec ky_1D;

    // Variables for vectorised matrices
    arma::vec A;
	arma::vec P;
	arma::vec KX;
	arma::vec KY;
	arma::vec K;
	arma::vec W;

    // Variables for irregular waves
    int specType_flag;
    int readPhases_flag;
    double gamma;
    double s;
    double rel_tol;
    double factor;
    double zero_order_moment;
    double first_order_moment;
    double second_order_moment;
    std::string waveDatabaseName;
    std::string file_path;
    std::string wavePhasesFileName;
    std::string filePhases_path;
    FILE *pfile_SPEC;
    FILE *pfile_TIME;

    // Variables for piecewise irregular waves
    int piecewise_flag;
    int num_pieces = 1;
    double time_gap;
    double time_piece;
    double spectrum_width;
    double df_max;
    double df_piece;
    double dw_piece;
    int num_comps_piece;
    int num_points_piece;
    int num_headings_piece;
    arma::vec time_ref;
    arma::vec time_ini;
    arma::vec time_end;
    arma::vec periods_piece;
    arma::vec freqs_piece;
    arma::vec ang_freqs_piece;
    arma::vec lambdas_piece;
    arma::vec k_piece;
    arma::mat kx_piece;
    arma::mat ky_piece;
    arma::vec headings_piece;
    arma::field<arma::mat> amplitudes_piece;
    arma::field<arma::mat> phases_piece;
    arma::vec headings_1D_piece;
    arma::vec kx_1D_piece;
    arma::vec ky_1D_piece;
    arma::field<arma::vec> amplitudes_1D_piece;
    arma::field<arma::vec> phases_1D_piece;

    // Declare class constructors
    Wave(Simulation *pSimInc, double H, double T, double D);

    // Methods

    virtual void GetWaveSpectrum(void) = 0;
    void CheckBreakingWave(void);
    void GetWaveLengths(void);
    std::tuple<arma::vec, arma::vec, arma::mat, arma::mat, arma::vec, arma::vec> GetWaveLengths(arma::vec periods, arma::vec headings);
    double solve_lambda(double T);
    double f_lambda(double lambda, double T);
    double df_lambda(double lambda, double T);
    void GetFreeSurface(void);
    arma::vec GetFreeSurface(arma::mat amplitudes, arma::mat phases, int num_points);
    arma::vec GetFreeSurface(double time, arma::vec x, arma::vec y);
    arma::vec GetPressure(double time, arma::vec x, arma::vec y, arma::vec z);
    void VectoriseComponentsMatrices(void);
    void SetSinglePiece(void);
    void SetZeroHeight(void);
    void WriteOut(std::string path);
};

class RegularWave : public Wave
{
public:
    RegularWave(Simulation *pSimInc, double H, double T, double D) : Wave(pSimInc, H, T, D) {};
    void GetWaveSpectrum(void);
};

class IrregularWave : public Wave
{
public:
    IrregularWave(Simulation *pSimInc, double H, double T, double D) : Wave(pSimInc, H, T, D) {};
    void GetWaveSpectrum(void);
    void GetCheckedFreeSurface(void);
    void ReadWaveTimeSeries(void);
    void ReadWavePSD(void);
    void GetJonswapSpectrum(void);
    arma::vec GetJonswapSpectrum(arma::vec freqs, double height, double period, double gamma);
    void GetSpreadingFunction(void);
    void CutSpectrumZeros(void);
    int CheckPhases(arma::vec t, arma::vec eta);

    // Methods for piecewise irregular waves
    void FindSpectrumWidth(void);
    void GetPiecesNumber(void);
    void GetPiecesWaveLengths(void);
    void CutPiecesSpectrumZeros(void);
    void GetPiecesTimeIntervals(void);
    void GetPiecesSpectra(void);
};

#endif // wavedef_hpp__