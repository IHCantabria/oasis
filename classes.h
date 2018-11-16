/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CLASE CatLine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <armadillo>

class MotherLine {
public:
	int nLine, nNodos, floor_flag;
	double L, dL, EA, beta, rho0, d, A, Cdt, Cdn, Cmn, CB, GK, GC, Gmu, Gvc, Dz;
	arma::mat tenAnch = arma::zeros(3,1), tenFair = arma::zeros(3,1);
	arma::mat posAnch = arma::zeros(3,1), posFair = arma::zeros(3,1);
	arma::mat pos, vel, acc, s, xc, zc, dxcds, dzcds, Te;
	double xF, zF, HF, VF, HA, VA, cosa, sina;
	void set_nLine(int n){nLine=n;}
	void leer_datosMoorings(void);
	void print_out(void);
	void initLine(void);
	void qs_Functions(double& ff,double& gg,double& DfDH,double& DfDV,double& DgDH,double& DgDV);
	void qs_GetTen(void);
	void qs_Solution(void);
	void write_out(void);
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

