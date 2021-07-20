
#ifndef WINDTURBINE_FLAG
#define WINDTURBINE_FLAG

#include <armadillo>
#include <string>
#include <OpenFAST.H>
#include "../Bodies/Bodies.hpp"

class Simulation;
class Body;

class WindTurbine {
public:

    int nTurbine; // Turbine number
    Simulation* pSim; // Pointer to simulation instance
    Body* pBody; // Pointer to corresponding body

    fast::OpenFAST FAST;
    fast::fastInputs fi;

    arma::mat forceOnBase = arma::zeros(6,1);

    WindTurbine(int n,Simulation* pSim){nTurbine=n;pSim = pSim;} // Inicializa un objeto de clase turbina dandole el indice
	void ReadPropertiesASCII(FILE*& pFile); // Lee inputs de las turbinas
    void SetUp(void); // Configura el objeto turbina
    void Finalize(void); // Configura el objeto turbina

    void ComputeForceOnBase(void);
    void UpdateNodesVelocities(void);

};

#endif