
#ifndef odedef_hpp__
#define odedef_hpp__
#include <armadillo>
#include <string>

#include "../Hydro/HydroDatabase.hpp"
#include "../Bodies/Bodies.hpp"
#include "../Lines/Lines.hpp"
#include "../Spring/Spring.hpp"
#include "../BCPs/BCPs.hpp"
#include "../BCPs/Winchies.hpp"

struct solver_data
{
	int nSistema, nSistema2;
	int nBcps;
	BCP **Bcps;
	int nLines;
	Line **Lines;
	int nSprings;
	Spring **Springs;
	int nBodies;
	Body **Bodies;
	HydroDatabase *Water;
	int nWinchies;
	Winchie **Winchies;
	int timeBufferSize = 1e2;
	arma::mat *timeBuffer;
	int *timeBufferCount;
};

class BDF
{
private:
	double h_0, h_1, h_2;
	double err, zz, sigma;
	arma::mat y_0, y_1, y_2;
	arma::mat yprime, dy;
	arma::mat I, J, LTE, EWT;
	arma::mat M, F, invM;
	int k, q;
	bool status;

public:
	double dt_max = 1e-2;
	double dt_min = 1e-9;
	double dt_ini = 1e-9;
	double atol = 1e-6;
	double rtol = 1e-3;
	int nIterMax = 10;

	int iJ = 0;

	int nSistema;
	double t, tmax, dt_out;
	arma::mat y;
	solver_data SD;
	Simulation *pSim;
	BDF(double t_u, double tmax_u, double dt_out_u, arma::mat y_u, Simulation *pIncSim);
	arma::mat fun(double, arma::mat);
	void Initialize(void);
	void step(void);
	void jac(double tt, arma::mat yy);
	arma::mat BDF2_fun(double t, arma::mat y);
};

#endif