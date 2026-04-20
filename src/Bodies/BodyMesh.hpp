
#ifndef bodymeshdef_hpp__
#define bodymeshdef_hpp__

#include <armadillo>
#include <string>
#include <cstdio>

// Class predefinition in order to avoid class cross-linking
class Body;
class Simulation;

// Define class BodyMesh and subclasses
class BodyMesh
{
private:
    int id;
    int typeMesh = 0;

public:
    // Main attributes
    int numElemNodes;                  // Number of nodes per element
    int numVertices;                   // Initial number of vertices of the body
    int maxEdgesLength;                // Maximum edges length
    int iniNumNodes, numNodes;         // Initial and final number of nodes of the body
    int iniNumElems, numElems;         // Initial and final number of elements of the body
    arma::mat iniNodes, nodes;         // Initial and final coordinates of nodes in the local body frame
    arma::umat iniElems, elems;        // Initial and final correspondence between nodes and elements
    arma::mat iniNormals, normals;     // Initial and final normal vectors to each element
    arma::vec iniJacobians, jacobians; // Initial and final jacobians of the transformations associated to each element
    arma::mat weightsJacNormal;        // Normal vectors times touching element's jacobian and corresponding weights

    Body* pBody;              // Pointer to the body linked to this mesh
    Simulation* pSim;         // Pointer to simulation instance
    std::string meshFileName; // Filename containing mesh data

    // Methods
    BodyMesh(int incId, std::string incMeshFileName, Body* incpBody);
    int GetId(void);
    virtual int GetType(void);
    virtual void Preprocess(void) = 0;
    void ReadPropertiesASCII(void);
    void TransformMesh(void);
    void CutMesh(double time);
    void IntegrateMesh(void);
};

class BodyTri2DMesh : public BodyMesh
{
private:
    int typeMesh = 1;

public:
    BodyTri2DMesh(int incId, std::string incMeshFileName, Body* incpBody)
        : BodyMesh(incId, incMeshFileName, incpBody) {};
    int GetType(void);
    void Preprocess(void);
};

#endif // bodymesh_hpp__