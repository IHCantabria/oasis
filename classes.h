/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CLASE CatLine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <armadillo>

class MotherLine {
public:
	int nLine, nNodos, p, N, floor_flag;
	double L, dL, dL0, EA, beta, rho0, d, A, Cdt, Cdn, Cmn, CB, GK, GC, Gmu, Gvc, Dz;
	arma::mat tenAnch = arma::zeros(3,1), tenFair = arma::zeros(3,1);
	arma::mat posAnch = arma::zeros(3,1), posFair = arma::zeros(3,1);
	arma::mat pos, vel, acc, F, s, xc, zc, dxcds, dzcds, Te, roots, weights;
	arma::mat C, D, MassMatrix, MM, StiffMatrix, MSMatrix, MassMatrix_diag;
	double xF, zF, HF, VF, HA, VA, cosa, sina;
	void set_nLine(int n){nLine=n;}
	void leer_datosMoorings(void);
	void print_out(void);
	void initLine(void);
	void qs_Functions(double& ff,double& gg,double& DfDH,double& DfDV,double& DgDH,double& DgDV);
	void qs_GetTen(void);
	void qs_Solution(void);
	void write_out(void);
	void SEM_getBaseFunctions(void);
	void SEM_coefficients(void);
	double SEM_poly(double x, int i);
	double SEM_poly_first_derivative(double x, int i);
	arma::mat SEM_get_D_local(void);
	void SEM_computeF(void);
};

class TensorLine: public MotherLine {
public:
	void initLine(void);

};

class TowingLine: public MotherLine {
public:
	void initLine(void);

};

class MooringLine: public MotherLine {
public:
	void initLine(void);

};

