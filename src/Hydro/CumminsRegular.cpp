
#include <armadillo>
#include "CumminsRegular.hpp"


CumminsRegular::CumminsRegular(int incId, Body** incBodies, Simulation* pIncSim, double incWaveHeight, double incWavePeriod, double incWaveHeading): HydroDatabase(incId, incBodies, pIncSim)
{
    waveHeight = incWaveHeight;
    wavePeriod = incWavePeriod;
    waveHeading = incWavePeriod;
}


double CumminsRegular::CalculateHydrodynamicForces(double time)
{
    // Get hydrostatic forces
    arma::mat hydrostatic_force = this->ComputeHydrostaticForces();

    // Calculate Damping effects
    arma::mat radiation_force = this->ComputeRadiationForces();

    // Calculate Wave exciting forces
    this->


}