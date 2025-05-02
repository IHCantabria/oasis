
#ifndef OWC_FLAG
#define OWC_FLAG

#include <armadillo>
#include <string>

class Simulation;
class Body;

class OWC
{
public:
    int idOWC;          // OWC number
    Simulation *pSim;   // Pointer to simulation instance
    Body *pBodyOWC;     // Pointer to corresponding body
    Body *pBodyFloater; // Pointer to linked body

    int idBodyOWC; // OWC body ID
    int idBodyFloater; // Floater body ID

    double waterplane_area;      // Waterplane area of the OWC
    double reference_air_volume; // Reference air volume of the OWC
    double hole_area;            // Hole area of the OWC
    double hole_discharge_coef;  // Hole discharge coefficient of the OWC
    int turbine_type;            // Type of turbine [-1: none, 0: hole, n>0: id of turbine in dataOWCTurbine.dat]
    arma::mat pos_local = arma::zeros(3,1); // Local position of the OWC in the floater body

    double displacement; // Current displacement of the OWC
    double velocity;     // Current velocity of the OWC
    double floater_displacement; // Current displacement of the floater
    double floater_velocity;     // Current velocity of the floater
    double air_pressure;    // Current pressure of the OWC
    double rel_pressure; // Current pressure difference of the OWC
    double rel_pressure_dot; // Current pressure difference derivative of the OWC
    double air_volume; // Current air volume of the OWC
    double air_volume_dot; // Current air volume derivative of the OWC
    double air_density; // Current air density of the OWC
    double air_mass; // Current air mass of the OWC
    double air_mass_dot; // Current air mass derivative of the OWC

    FILE *pfile_output; // Pointer to output file

    OWC(int n, Simulation *pSimInp)
    {
        idOWC = n + 1;
        pSim = pSimInp;
    }
    void ReadPropertiesASCII(FILE *pFile);
    void Initialize(FILE *pFile);
    void Finalize(void);
    void WriteOut(double t);

    void ComputeForces(double time);
    void ComputePressure(double time);
    double ComputeAirDensity(double pressure);
    double ComputeHoleMassFlowRate(double pressure_diff, double stagnation_density, double hole_area, double hole_discharge_coef);
    double ComputeRelPressureRate(double air_pressure, double mass_flow_rate, double air_mass, double air_volume_dot, double air_volume);
};

#endif