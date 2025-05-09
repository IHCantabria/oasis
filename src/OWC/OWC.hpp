
#ifndef OWC_FLAG
#define OWC_FLAG

#include <armadillo>
#include <string>

class Simulation;
class Body;

class OWCTurbineType
{
public:
    // OWC turbine type number
    int idOWCTurbineType;

    // Pointer to the simulation instance
    Simulation *pSim;

    // Properties of the turbine type
    int valve_type;                       // Type of valve [0: None, 1: Bypass, 2: Throttle, 3: Bypass + Throttle]
    double inertia;                       // Inertia of the turbine
    double damping;                       // Damping of the turbine
    double rotor_diameter;                // Rotor diameter of the turbine
    double bypass_valve_area;             // Bypass valve area of the turbine
    double bypass_valve_discharge_coef;   // Bypass valve discharge coefficient of the turbine
    double throttle_gap_air_volume;       // Throttle gap air volume of the turbine
    double throttle_valve_area;           // Throttle valve area of the turbine
    double throttle_valve_discharge_coef; // Throttle valve discharge coefficient of the turbine
    int controller_type;                  // Type of controller [0: None, 1: Algebraic, 2: Custom]
    double generator_rated_power;         // Rated power of the generator
    double generator_max_torque;          // Maximum torque of the generator
    // Buckling curves
    arma::vec buck_curve_adim_pressure;       // Adimensional pressure curve of the turbine
    arma::vec buck_curve_adim_power;          // Adimensional power curve of the turbine
    arma::vec buck_curve_adim_mass_flow_rate; // Adimensional mass flow rate curve of the turbine
    arma::vec buck_curve_efficiency;          // Efficiency curve of the turbine
    // Custom controller parameters
    // TODO: Implement custom controller parameters

    // Preprocessed data
    double adim_pressure_bep;            // Adimensional pressure at best efficiency point
    double adim_power_bep;               // Adimensional power at best efficiency point
    double adim_critical_pressure;       // Adimensional critical pressure of the turbine
    double adim_critical_mass_flow_rate; // Adimensional critical mass flow rate of the turbine
    // Algebraic controller parameters
    double param_a; // Parameter a of the controller
    double param_b; // Parameter b of the controller

    // Class methods
    OWCTurbineType(int n, Simulation *pSimInp);

    void ReadPropertiesASCII(FILE *pFile);
    void Initialize(FILE *pFile);
    void Finalize(void);
};

class OWCTurbine
{
public:
    OWCTurbineType *pOWCTurbineType; // Pointer to the corresponding turbine type
    Simulation *pSim;                // Pointer to the simulation instance

    double mass_flow_rate;          // Mass flow rate of the turbine
    double angular_velocity;        // Angular velocity of the turbine
    double angular_acceleration;    // Angular acceleration of the turbine
    double air_torque;              // Air torque of the turbine
    double generator_torque;        // Generator torque of the turbine
    double power;                   // Power of the turbine
    double critical_pressure;       // Critical pressure of the turbine
    double critical_mass_flow_rate; // Critical mass flow rate of the turbine
    double gap_rel_pressure;        // Relative air pressure in the gap of the turbine
    double gap_rel_pressure_dot;    // Relative air pressure derivative in the gap of the turbine
    double gap_mass_flow_rate;      // Mass flow rate in the gap of the turbine
    double gap_air_pressure;        // Air pressure in the gap of the turbine
    double gap_air_density;         // Air density in the gap of the turbine
    double gap_air_mass;            // Air mass in the gap of the turbine

    // Class methods
    OWCTurbine(Simulation *pSimInp, int id_OWCTurbineType);

    void Initialize(double initial_angular_velocity);

    bool IsAdimPressureInRange(double adim_pressure);
    void ComputeMassFlowRate(double pressure_diff, double stagnation_density);
    void ComputeAirTorque(double pressure_diff, double stagnation_density);
    void ComputeCriticalPressure(double stagnation_density);
    void ComputeCriticalMassFlowRate(double stagnation_density);
    void ComputeAngularAcceleration(double pressure_diff, double stagnation_density);
    void ComputeGenTorque(void);
};

class OWC
{
public:
    int idOWC;               // OWC number
    Simulation *pSim;        // Pointer to simulation instance
    Body *pBodyOWC;          // Pointer to corresponding body
    Body *pBodyFloater;      // Pointer to linked body
    OWCTurbine *pOWCTurbine; // Pointer to the corresponding turbine

    int idBodyOWC;     // OWC body ID
    int idBodyFloater; // Floater body ID

    double waterplane_area;                  // Waterplane area of the OWC
    double reference_air_volume;             // Reference air volume of the OWC
    double hole_area;                        // Hole area of the OWC
    double hole_discharge_coef;              // Hole discharge coefficient of the OWC
    int turbine_type;                        // Type of turbine [-1: none, 0: hole, n>0: id of turbine in dataOWCTurbine.dat]
    double turb_initial_angular_velocity;    // Initial angular velocity of the turbine
    arma::mat pos_local = arma::zeros(3, 1); // Local position of the OWC in the floater body

    double displacement;         // Current displacement of the OWC
    double velocity;             // Current velocity of the OWC
    double floater_displacement; // Current displacement of the floater
    double floater_velocity;     // Current velocity of the floater
    double air_pressure;         // Current pressure of the OWC
    double rel_pressure;         // Current pressure difference of the OWC
    double rel_pressure_dot;     // Current pressure difference derivative of the OWC
    double air_volume;           // Current air volume of the OWC
    double air_volume_dot;       // Current air volume derivative of the OWC
    double air_density;          // Current air density of the OWC
    double air_mass;             // Current air mass of the OWC
    double air_mass_dot;         // Current air mass derivative of the OWC

    double gap_pressure_diff;    // Throttle valve gap pressure difference of the OWC
    double gap_rel_pressure;     // Throttle valve gap relative pressure of the OWC
    double gap_rel_pressure_dot; // Throttle valve gap relative pressure derivative of the OWC

    FILE *pfile_output; // Pointer to output file

    OWC(int n, Simulation *pSimInp);

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