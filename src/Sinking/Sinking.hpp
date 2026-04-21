// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef sinkingdef_hpp__
#define sinkingdef_hpp__

#include <armadillo>
#include <string>
#include <cstdio>
#include <yaml-cpp/yaml.h>
#include "../Hydro/HydroForce.hpp"

// Attribute class objects forward declaration
class Body;
class Simulation;

class Sinking
{
public:
    double rhoW;

    Simulation* pSim;   // Pointer to simulation instance
    Body* pSinkingBody; // Array to store the sinking body object

    int indBody; // index of the body in the general database
    int numGroups;
    arma::field<arma::mat> groupsPoints;        // [numGroups, 1] [numPoints,2]
    arma::mat groupsAreas;                      // [numGroups, 1]
    arma::mat groupsCenters;                    // [numGroups, 3]
    arma::mat groupsIx;                         // [numGroups, 1]
    arma::mat groupsIy;                         // [numGroups, 1]
    arma::field<arma::mat> groupsFillingTimes;  // [numGroups, 1] [numTimes 1]
    arma::field<arma::mat> groupsFillingStates; // [numGroups, 1] [numTimes 1]

    int numHDBs;
    arma::mat InterpMasses; // Vector of filling masses for which the hydroforce objects are provided
    HydroDatabase** pHydro; // Vector of hydroforce objects for the different filling masses
    int indHydro1, indHydro2;
    double hydroInterpCoef;

    double totalFillingMass;
    double bodyReferenceMass;
    arma::mat bodyReferenceMassMat = arma::zeros(6, 6);
    arma::mat groupsMasses; // [numGroups, 1]
    arma::mat groupsCOG = arma::zeros(3, 1);
    arma::mat groupsInertia = arma::zeros(6, 6);

    FILE* pfile_FillingCOG;
    FILE* pfile_FillingCOG_csv = nullptr;

    // Declare constructors
    Sinking(int n, Simulation* pSim_inp);

    // Methods definition
    void ReadPropertiesASCII(FILE* filePointer, std::string inputFolderPath); // Read sinking properties
    void ReadPropertiesYAML(YAML::Node node, std::string inputFolderPath);
    void UpdateGroupsFillingState(double t);
    void UpdateInterpHydro(void);
    void UpdateBodyProperties(void);
    void UpdateSinkingHydrodynamics(double t);
    void UpdateSinkingHydrostatics(double t);

    void OpenOutputFilesASCII(std::string path);
    void OpenOutputFilesCSV(std::string path);
    void CloseOutputFilesASCII(void);
    void CloseOutputFilesCSV(void);
    void OpenOutputFiles(std::string path);
    void CloseOutputFiles(void);
    void WriteOut(double t); // Escribir datos a fichero
};

#endif // sinkingdef_hpp__