// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef seafloordef_hpp__
#define seafloordef_hpp__

#include <armadillo>
#include <string>
#include <yaml-cpp/yaml.h>

class SeaFloor
{
private:
    int id;               // Body index, used as unique identifier and in the constructor
    int seaFloorType = 0; // Placeholder value; will be set by the concrete subclass

    // Atributes
public:
    std::string meshFileName; // Mesh filename for variable bathymetry
    // For variable bathymetry
    int numCloudPoints;
    int numTriangles;
    int flagBarycenter; // proximamente a borrar
    double seabedDepth;

    // For the inclined plane
    arma::mat p1 = arma::zeros(1, 3);
    arma::mat p2 = arma::zeros(1, 3);
    arma::mat p3 = arma::zeros(1, 3);

    // Methods
    SeaFloor(int incId);
    int GetId(void);
    virtual int GetType(void) = 0;
    virtual arma::field<arma::mat> projectPoints(arma::mat nodes) = 0;
    // Methods that must be defined for bathymetry
    void ReadPropertiesASCII(FILE*& pFilePointer);
    void ReadPropertiesYAML(YAML::Node node);
};

class Bathymetry : public SeaFloor
{
private:
    int seaFloorType = 3;
    // Private method only used internally within the class
    arma::mat vertexCoordinates(int triangle);

public:
    // atributos
    int numCloudPoints;
    int numTriangles;
    int flagBarycenter;
    arma::mat pointMatrix;
    arma::umat triangleMatrix;
    arma::mat vertexNormals;
    arma::mat barycenter;
    arma::field<arma::mat> projectionMatrix;
    arma::field<arma::mat> changeFrameMatrix;
    arma::field<arma::mat> normalsTriangle;

    // Metodos
    Bathymetry(int incId)
        : SeaFloor(incId) {};
    void ReadPropertiesASCII(FILE*& pFilePointer, std::string inputFilePath);
    void ReadPropertiesYAML(YAML::Node node, std::string inputFilePath);
    void getVertexNormals(void);
    void getProjectionMatrix(void);
    arma::field<arma::mat> projectPoints(arma::mat nodes);
    arma::uvec closerTriangles(arma::mat point);
    int GetType(void);
};

class Inclined : public SeaFloor
{
private:
    int seaFloorType = 2;

public:
    // atributes
    arma::mat planeNormal = arma::zeros(1, 3);
    double a, b, c, d; // Coefficients defining the plane equation
    // methods
    Inclined(int incId)
        : SeaFloor(incId) {};
    void getPlaneEquation(void);
    arma::field<arma::mat> projectPoints(arma::mat nodes);
    int GetType(void);
};

class Flat : public SeaFloor
{
private:
    int seaFloorType = 1;

public:
    Flat(int incId)
        : SeaFloor(incId) {};
    arma::field<arma::mat> projectPoints(arma::mat nodes);
    int GetType(void);
};

#endif