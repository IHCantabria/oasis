
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

    int numWindTurbines; // Turbine number
    Simulation* pSim; // Pointer to simulation instance
    Body** pBodies; // Pointers to corresponding bodies

    fast::OpenFAST FAST;
    fast::fastInputs fi;

    arma::mat forceOnBase = arma::zeros(6,1);

    WindTurbine(int n,Simulation* pSimInp){numWindTurbines=n; pSim = pSimInp;} // Inicializa un objeto de clase turbina dandole el indice
	void ReadPropertiesASCII(FILE*& pFile); // Lee inputs de las turbinas
    void Initialize(void); // Configura el objeto turbina
    void Step(void);
    void Finalize(void); // Cierra el caso

    void GetBaseForces(void);
    void SetBaseMovements(void);

};

#endif