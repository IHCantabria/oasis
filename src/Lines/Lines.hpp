
#include <armadillo>
#include <string>
#include "../BCPs/BCPs.hpp"

class Line {
private:
	double fg, fs, fd;
	arma::mat v, vt, vn;
	arma::mat FF, ff, t; 
	arma::mat drds, drdsdt, norm_drds, dedt, T;
	arma::mat e_z;
public:
	int lineType, nLine, nNodos, p, N, floor_flag, BCP_1, BCP_N;
	double L, dL, dL0, EA, beta, rho0, d, A, Cdt, Cdn, Cmn, CB, GK, GC, Gmu, Gvc, Dz, Kn;
	arma::mat ten_1 = arma::zeros(3,1), ten_N = arma::zeros(3,1);
	arma::mat pos_1 = arma::zeros(3,1), pos_N = arma::zeros(3,1);
	arma::mat pos, vel, acc, F, s, xc, zc, dxcds, dzcds, Te, roots, weights;
	arma::mat C, D, MassMatrix, MM, StiffMatrix, MSMatrix, MassMatrix_diag;
	arma::sp_mat D_sp, MassMatrix_sp, MM_sp, StiffMatrix_sp, MSMatrix_sp;
	BCP * LineBCP [2];
	double xF, zF, HF, VF, HA, VA, cosa, sina;
	void set_nLine(int n){nLine=n;}
	void leer_datosLines(void);
	void print_out(void);
	void initLine(void);
	void qs_Functions(double& ff,double& gg,double& DfDH,double& DfDV,double& DgDH,double& DgDV);
	void qs_GetTen(void);
	void qs_Solution(void);
	void write_out(double t);
	void SEM_getBaseFunctions(void);
	void SEM_coefficients(void);
	double SEM_poly(double x, int i);
	double SEM_poly_first_derivative(double x, int i);
	arma::mat SEM_get_D_local(void);
	void SEM_computeF(void);
};