
#ifndef WINCHIE_CONTROLLER_HORIZONTAL_FLAG
#define WINCHIE_CONTROLLER_HORIZONTAL_FLAG
#include "WinchiesController.hpp"

class WinchieControllerHorizontal : public WinchieController
{
public:
	int k = 0;		// contador
	int indBody;	// Index of the body that is controlled with the winches

	double Ac, Bc, Cc, Dc;			  // Controller state-space model parameters
	arma::mat Kc = arma::zeros(3, 1); // Lead controller gains
	arma::mat Ki = arma::zeros(3, 1); // Integral time gains
	double Af, Bf, Cf, Df;			  // Filter state-space model parameters
	double Ar, Br, Cr, Dr;			  // Reference state-space model parameters
	arma::mat ur = arma::zeros(3, 1); // Vector of reference positions x(m), y(m), heading(rad)
	arma::mat yb = arma::zeros(3, 1); // Vector of actual positions of the body
	arma::mat error;				  // Mat with timeseries with the error

	arma::mat xc = arma::zeros(3, 1); // Variables de estado del controlador
	arma::mat yc = arma::zeros(3, 1); // Output del controlador
	arma::mat xf = arma::zeros(3, 1); // Variables de estado del filtro
	arma::mat yf = arma::zeros(3, 1); // Output del filtro

	int reference_flag;				  // Flag for reference type [0: state-space; n<=1: n discrete points]
	arma::mat xr = arma::zeros(3, 1); // Variables de estado de la referencia
	arma::mat yr = arma::zeros(3, 1); // Output de la referencia
	arma::mat t_ref, x_ref, y_ref, yaw_ref;

	double Kw; // Ganancia de los winchies
	double time_ini;

	int inversor_flag;
	double T_max, T_min;

	// Straigt lines inversor atributes
	int nIterMax;
	double atol;
	arma::mat posAnchG, posFairL;

	// Coefficients inversor atributes
	arma::uvec ind_x_pos, ind_x_neg, ind_y_pos, ind_y_neg, ind_g_pos, ind_g_neg;
	arma::mat coef_x_pos, coef_x_neg, coef_y_pos, coef_y_neg, coef_g_pos, coef_g_neg;

	// Output files specific to horizontal controller
	FILE *pfile_FC;
	FILE *pfile_RP;

	WinchieControllerHorizontal(int n, Winchie **Ws, Simulation *pIncSim);
	~WinchieControllerHorizontal() override;

	void ReadPropertiesASCII(FILE *pFile) override;
	void ReadPropertiesYAML(YAML::Node node) override;
	void SetUpWinchiesController(void) override;
	void controlWinchies(double time) override;
	void inversorBlock(void);
	void OpenOutputFilesASCII(std::string path) override;
	void CloseOutputFilesASCII(void) override;
	void WriteOut(double t) override;
};

#endif
