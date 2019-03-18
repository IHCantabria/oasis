

#include <armadillo>
#include <string>
#include "math_lib.h"



class BCP {
public:
	int nBCP; // Indice identificador del punto de condicion de contorno
	int nLinesBCP; // Numero de lineas que confluyen en el punto
	int * BCPLineIndex; // Array con los indices identificadores de las lineas que confluyen en el punto
	int * BCPLineNode; // Array de flags que, para cada linea ii que confluye al punto, indica si la linea confluye al nodo 1 (BCPLineNode[ii]=1) o al nodo N (BCPLineNode[ii]=2)
	std::string fileName;
	arma::mat pos = arma::zeros(3,1); // Posicion del punto
	arma::mat vel = arma::zeros(3,1); // Velocidad del punto
	arma::mat acc = arma::zeros(3,1); // Aceleracion del punto
	void set_nBCP(int n){nBCP=n;}
	void leer_datosBCPs(void);
	virtual void getValues(double t) =0;
};

class AnchorBCP: public BCP {
public:	
	void getValues(double t);
};


class FairleadBCP: public BCP {
public:	
	void getValues(double t);
	arma::mat posF;
	spline x_spl, y_spl, z_spl;
};

class Line {
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

struct solver_data{
	int nLines, nSistema, nSistema2;
	Line * Lines;
};

class BDF{
private:
	double dt_max   = 1e-2;
	double dt_min   = 1e-8;
	double dt_ini   = 1e-4;
	double eta_min  = 1e-3;
	double eta_max  = 1e-2;
	double atol     = 1e-6;
	double rho      = 1e+1;
	double sigma    = 1e-2;
	int    nIterMax = 100;
public:
	int nSistema;
	double t, t_prev, tmax, dt, eta;
	arma::mat y, y_prev, y_jac, J;
	solver_data SD;
	arma::mat (*fun) (double, arma::mat, solver_data);
	arma::mat (*jac) (double, arma::mat, solver_data);
	BDF(void){nSistema=0;};
	BDF(
		double t_u,
		double tmax_u,
		arma::mat y_u, 
		arma::mat (*fun_u) (double, arma::mat, solver_data), 
		arma::mat (*jac_u) (double, arma::mat, solver_data),
		solver_data SD_u
		){
			t = t_u;
			tmax = tmax_u;
			y = y_u;
			nSistema = y_u.n_rows;
			fun = fun_u;
			jac = jac_u;
			SD = SD_u;
			dt = dt_ini;
	}
	void step(void);
	void get_next_dt(void);
};