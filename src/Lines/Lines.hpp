
#ifndef linedef_hpp__
#define linedef_hpp__
#include <armadillo>
#include <string>
#include "../BCPs/BCPs.hpp"
#include "../SeaFloor/SeaFloor.hpp"

class Line {
private:
	int id;
	double fg, fs, fd;
	arma::mat v, vt, vn;
	arma::mat FF, ff, t; 
	arma::mat drds, drdsdt, norm_drds, dedt, T;
	arma::mat e_z;
	
public:
	int numBcps=2;
	int indexBcps [2];
	int lineType, nLine, nNodos, p, N, floor_flag, BCP_1, BCP_N, flag_tension, flag_stiffness, smoothstep, frictionModel;
	double L, dL, dL0, EA, beta, rho0, d, A, Cdt, Cdn, Cmn, CB, GK, GC, Kn, dampCoef, VR, vth, ust, usn, ud, deltamax, fn;
	double paramNormal, parammuelle1, parammuelle2, paramVel;
	arma::mat strain_data, stress_data;
	arma::mat ten_1 = arma::zeros(3,1), ten_N = arma::zeros(3,1);
	arma::mat pos_1 = arma::zeros(3,1), pos_N = arma::zeros(3,1);
	arma::mat a_1 = arma::zeros(4,1);// Variables para almacenar los coeficientes polinomicos del coeficiente de friccion
	arma::mat a_2 = arma::zeros(4,1);// Variables para almacenar los coeficientes polinomicos del coeficiente de friccion
	arma::mat pos, vel, acc, F, s, xc, zc, dxcds, dzcds, Te, roots, weights, isSlip, posFriccion;
	arma::mat C, D, MassMatrix, MM, StiffMatrix, MSMatrix, MassMatrix_diag;
	arma::mat inv_MM, inv_MM_1, inv_MM_N, inv_MM_1N;
	arma::sp_mat D_sp, MassMatrix_sp, MM_sp, StiffMatrix_sp, MSMatrix_sp;
	arma::mat projectedPoints, projectionDirection, zCoordinates;
	BCP* pLineBcps [2];
	arma::uvec ind4CouplingMat;
	arma::mat F_1, F_N;
	int first_node, last_node;
	double xF, zF, HF, VF, HA, VA, cosa, sina;
	int indexSeaFloor;
	SeaFloor* pLineSeaFloor;

	double g;
	double rhoW;
	double fondo;

	FILE* pfile_xpos;
	FILE* pfile_ypos;
	FILE* pfile_zpos;
	FILE* pfile_ten;	
	FILE* pfile_ten_line;

	Line(int incId, double incG, double incRhoW, double incFondo);
	int GetId();
	void ReadPropertiesASCII(FILE* pFilePointer);
	void OpenOutputFilesASCII(std::string path);
	void CloseOutputFilesASCII (void);
	void print_out(void);
	void initLine(void);
	void qs_Functions(double& ff,double& gg,double& DfDH,double& DfDV,double& DgDH,double& DgDV);
	void qs_GetTen(void);
	void qs_Solution(void);
	void WriteOut(double t);
	void SEM_getBaseFunctions(void);
	void SEM_coefficients(void);
	double SEM_poly(double x, int i);
	double SEM_poly_first_derivative(double x, int i);
	arma::mat SEM_get_D_local(void);
	void SEM_computeF(void);
	
};


#endif