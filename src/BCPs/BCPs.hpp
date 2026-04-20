
#ifndef bcpdef_hpp__
#define bcpdef_hpp__
#include <armadillo>
#include <string>
#include <cstdio>
#include <yaml-cpp/yaml.h>

// Class predefinition in order to avoid class cross-linking
class Body;
class Line;

// Define class BCP and subclasses
class BCP
{
private:
    int id;
    int typeBcp = 0;

public:
    // Conexion attributes
    int countBody = 0;
    int countLine = 0;
    Line** pLines;        // Array of pointers to the lines connected at this point
    Body** pBodies;       // Array of pointers to the bodies connected at this point
    int numBodiesBcp = 0; // Number of bodies connected to the BCP
    int numLinesBcp = 0;  // Number of lines connected at this point
    int winchId = 0;      // Id of the winch connected to the instance of the BCP

    // Attributes common to all BCPs
    double tBCP = 0.0;
    int* pBcpLineIndex; // Array of line identifier indices connected at this point
    int* pBcpLineNode;  // Array of flags: for each line ii connected to this point, indicates whether the line connects
                       // at nodo 1 (pBcpLineNode[ii]=0) o al nodo N (pBcpLineNode[ii]=1)
    arma::mat pos = arma::zeros(3, 1); // Point position
    arma::mat vel = arma::zeros(3, 1); // Point velocity
    arma::mat acc = arma::zeros(3, 1); // Point acceleration

    // Fairlead-specific attributes
    std::string actuatorFileName;
    arma::mat tF, posF, velF, accF;
    arma::mat pos0; // Initial position of the point
    int ni = 0;
    double dt;

    // Joint-specific attributes
    int iLJ = 0;
    arma::mat posLines;
    arma::mat velLines;
    arma::mat JointForce = arma::zeros(1, 3);
    double mass_Joint, vol_Joint, rad_Joint, sec_Joint;
    double g, rhoW, seabedDepth;
    int flag_assigned = 0;
    int flag_counted = 0;
    arma::uword couplingMatIndex;

    // Body BCP-specific attributes
    arma::mat posWrtCdgLocal = arma::zeros(3, 1);  // Brazo COG-bcp en local
    arma::mat posWrtCdgGlobal = arma::zeros(3, 1); // Brazo COG-bcp en global
    arma::mat rotMat = arma::eye(3, 3);            // Rotation matrix of the body
    arma::mat rotMat_dot = arma::zeros(3, 3);      // Time derivative of the body rotation matrix
    arma::mat posG_BCP = arma::zeros(3, 1);        // BCP position in global frame (3 DOF)
    arma::mat velG_BCP = arma::zeros(3, 1);        // BCP velocity in global frame (3 DOF)
    arma::mat accG_BCP = arma::zeros(3, 1);        // BCP acceleration in global frame (3 DOF)
    arma::mat forceBcp = arma::zeros(6, 1);        // Fuerzas y momentos que actuan sobre el BCP en 6 dofs

    arma::mat temp = arma::zeros(6, 1); // Para pintado de variables y debugeo

    // Methods
    BCP(int incId);
    int GetId(void);
    virtual int GetType(void);
    virtual void GetValues(double t) = 0;
    virtual void Initialize(void);
    virtual void Print(void);
    void ReadPropertiesASCII(FILE*& pFilePointer);
    void ReadPropertiesYAML(YAML::Node node);
    virtual void UpdateBoundary();
};

class AnchorBCP : public BCP
{
private:
    int typeBcp = 1;

public:
    AnchorBCP(int incId)
        : BCP(incId) {};
    int GetType(void);
    void GetValues(double t);
};

class BodyBCP : public BCP
{
private:
    int typeBcp = 4;

public:
    BodyBCP(int incId)
        : BCP(incId) {};
    int GetType(void);
    void GetValues(double t);
};

class FairleadBCP : public BCP
{
private:
    int typeBcp = 2;

public:
    FairleadBCP(int incId)
        : BCP(incId) {};
    int GetType(void);
    void GetValues(double t);
    void Initialize(std::string folder_path);
    void ReadPropertiesASCII(FILE*& pFilePointer, std::string inputFilePath);
    void ReadPropertiesYAML(YAML::Node node, std::string inputFilePath);
    void Print(void);
};

class JointBCP : public BCP
{
private:
    int typeBcp = 3;

public:
    JointBCP(int incId)
        : BCP(incId) {};
    int GetType(void);
    void GetValues(double t);
    void Initialize(double incG, double incRhoW, double incFondo);
};

class ElasticAnchorBCP : public BCP
{
private:
    int typeBcp = 5;

public:
    // Exponential restoring force parameters: F = -c*(1 - exp(-k*x)) toward reference position
    arma::vec pos_ref = arma::zeros(3, 1); // Reference/equilibrium position [m]
    double c_param = 0.0;                  // Ultimate holding capacity [N]
    double k_param = 0.0;                  // Rate parameter [1/m]
    double anchor_mass = 0.0;              // Anchor mass [kg]
    double anchor_vol = 0.0;               // Anchor volume [m³] (for drag calculation)
    double rad_anchor = 0.0;               // Anchor radius [m] (computed from volume)
    double sec_anchor = 0.0;               // Cross-sectional area [m²] (for drag)
    double g = 0.0;                        // Gravity [m/s²]
    double rhoW = 0.0;                     // Water density [kg/m³]

    // Constructor
    ElasticAnchorBCP(int incId)
        : BCP(incId) {};

    // Methods
    int GetType(void);
    void GetValues(double t);
    void Initialize(double incG, double incRhoW);
    arma::vec ComputeRestoringForce(void); // Compute F = -c*(1 - exp(-k*x)) restoring force
    void Print(void);
};

#endif // bcp_hpp__