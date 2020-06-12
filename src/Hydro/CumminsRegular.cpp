
#include <armadillo>
#include "CumminsRegular.hpp"
#include "../Simulations/Simulation.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"


arma::mat CumminsRegular::ComputeFirstWaveExcForce()
{
	// Create local variables
	double time = pSim->pTimeSolver->t;
	// Get First Order Wave exciting data from storage
	arma::mat wave_exc_mag = pWaveExcitingMag->subcube(0, numPeriodExc, numHeadingExc, 5, numPeriodExc, numHeadingExc);
	arma::mat wave_exc_pha = pWaveExcitingPha->subcube(0, numPeriodExc, numHeadingExc, 5, numPeriodExc, numHeadingExc);
	// Calculate wave force
	double time_slope = 1/(*pFrequencies)(0, numPeriodExc);
	arma::mat wave_force = arma::zeros(6, 1);
	double angular_freq = 2*M_PI*(*pFrequencies)(0, numPeriodExc);
	for (int i=0; i<6; i++)
	{
		wave_force(i, 0) = wave_exc_mag(i, 0)*waveAmplitude*cos(angular_freq*time + wave_exc_pha(i, 0));

		if (time < time_slope)
		{
			wave_force(i, 0) = wave_force(i, 0)*time/time_slope;
		}
	}
	//std::cout << "Wave Amplitude: " << waveAmplitude << std::endl;

	return wave_force;
}


CumminsRegular::CumminsRegular(int incId, int incIdBody, Body** incBodies, Simulation* pIncSim, double incWaveHeight, double incWavePeriod, double incWaveHeading): HydroDatabase(incId, incIdBody, incBodies, pIncSim)
{
    waveHeight = incWaveHeight;
    wavePeriod = incWavePeriod;
    waveHeading = incWavePeriod;
}


arma::mat CumminsRegular::CalculateHydrodynamicForces(double time)
{
    // Get hydrostatic forces
    // arma::mat hydrostatic_force = this->CalculateHydrostaticForces();

    // Calculate Damping effects
    arma::mat radiation_force = this->ComputeRadiationForces();

    // Calculate Wave exciting forces

}