#include "WinchiesControllerConstantTension.hpp"
#include "../Simulations/Simulation.hpp"
#include "../Exceptions/Exception.hpp"

WinchieControllerConstantTension::WinchieControllerConstantTension(int n, Winchie **Ws, Simulation *pIncSim)
	: WinchieController(n, Ws, pIncSim)
{
	targetTension = 0.0;
}

WinchieControllerConstantTension::~WinchieControllerConstantTension() {}

void WinchieControllerConstantTension::ReadPropertiesASCII(FILE *pFile)
{
	// NOTE: The 3 header lines and controller type have already been read by the factory.
	// We read the target tension directly.
	char buffer_line[1000];
	fscanf(pFile, "%lf %[^\n]\n", &targetTension, buffer_line);
}

void WinchieControllerConstantTension::SetUpWinchiesController(void)
{
	T = arma::ones(nWinchies, 1) * targetTension;
	std::cout << "Winches controller: Constant tension = " << targetTension << " N" << std::endl;
	applyTensions();
}

void WinchieControllerConstantTension::controlWinchies(double time)
{
	// Tension is constant — just reapply each step
	applyTensions();
}
