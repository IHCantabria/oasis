
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>

#include "SeaFloor.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"

SeaFloor::SeaFloor(int incId)
{
    id = incId;
}

int SeaFloor::GetId(void)
{
    return this->id;
}

void SeaFloor::ReadPropertiesASCII(FILE *&pFilePointer)
{
    // Declare local variables
    std::string header_check;
    char buffer_line[1000];
    char cmeshFileName[1000];

    // Ignore header lines
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), pFilePointer);
        header_check = buffer_line;
        if (header_check.substr(0, 3).compare("///"))
        {
            std::stringstream ss;
            ss << "Error while parsing file: dataSeaFloor.dat - HINT FLOOR ID: " << this->GetId() << " - Please check";
            throw IOError(ss.str());
        }
    }
    if (this->GetType() == 3)
    {
        fscanf(pFilePointer, "%s %[^\n]\n", &cmeshFileName, buffer_line);
        meshFileName = cmeshFileName;
    }
    if (this->GetType() == 2)
    {
        // Save the positions of the three points defining the plane
        fscanf(pFilePointer, "%lf %lf %lf %[^\n]\n", &p1(0, 0), &p1(0, 1), &p1(0, 2), buffer_line);
        fscanf(pFilePointer, "%lf %lf %lf %[^\n]\n", &p2(0, 0), &p2(0, 1), &p2(0, 2), buffer_line);
        fscanf(pFilePointer, "%lf %lf %lf %[^\n]\n", &p3(0, 0), &p3(0, 1), &p3(0, 2), buffer_line);
    }
    if (this->GetType() == 1)
    {
        fscanf(pFilePointer, "%lf %[^\n]\n", &fondo, buffer_line);
    }
}

// Flaf floor methods

int Flat::GetType(void)
{
    return this->seaFloorType;
}

arma::field<arma::mat> Flat::projectPoints(arma::mat nodos)
{
    arma::field<arma::mat> output(1, 3);
    // INICIALIZACION DE NODOS
    int numNodos = nodos.n_rows;
    arma::mat projected_points = nodos; // aqui se devolveran los puntos ya proyectados
    // INICIALIZACION DE Z
    arma::mat zCoordinates = arma::zeros(numNodos, 1);
    // INICIALIZACION DE LAS NORMALES  RETORNAR
    arma::mat normales = arma::zeros(numNodos, 3);
    normales.col(2) = arma::ones(numNodos, 1);
    projected_points.col(2) = fondo * arma::ones(numNodos, 1);
    zCoordinates = nodos.col(2) - fondo * arma::ones(numNodos, 1);
    output(0, 0) = projected_points;
    output(0, 2) = normales;
    output(0, 1) = zCoordinates;
    return output;
}

// INCLINED DEFINITION
// ¿No seria mejor que los puntos estuvieran en el constructor, que este calculara la normal?
void Inclined::getPlaneEquation(void)
{
    arma::mat Lado1 = p2 - p1;
    arma::mat Lado2 = p3 - p1;
    arma::mat prodCross = arma::cross(Lado1, Lado2);
    normalPlano = prodCross / arma::norm(prodCross, 2);
    double tol = 1e-10;
    if (arma::norm(prodCross, 2) < tol)
    {
        std::stringstream ss;
        ss << "The input points are aligned\n";
        throw ValueError(ss.str());
    }
    if (normalPlano(2) < -1e-10)
    {
        normalPlano = -normalPlano;
    }
    this->a = arma::as_scalar(normalPlano(0));
    this->b = arma::as_scalar(normalPlano(1));
    this->c = arma::as_scalar(normalPlano(2));
    this->d = a * arma::as_scalar(p1(0)) + b * arma::as_scalar(p1(1)) + c * arma::as_scalar(p1(2));
}

int Inclined::GetType(void)
{
    return this->seaFloorType;
}

arma::field<arma::mat> Inclined::projectPoints(arma::mat nodos)
{
    arma::field<arma::mat> aRetornar(1, 3);
    // INICIALIZACION DE NODOS
    int numNodos = nodos.n_rows;
    arma::mat projected_points = arma::zeros(numNodos, 3); // aqui se devolveran los puntos ya proyectados
    // INICIALIZACION DE Z
    arma::mat zCoordinates = arma::zeros(numNodos, 1);
    // INICIALIZACION DE LAS NORMALES  RETORNAR
    arma::mat normales = arma::zeros(numNodos, 3);
    double t; // parametricas de la recta que une el punto con su proyeccion
    for (int i = 0; i < numNodos; i++)
    {
        t = (d - a * arma::as_scalar(nodos(i, 0)) - b * arma::as_scalar(nodos(i, 1)) - c * arma::as_scalar(nodos(i, 2))) / (a * a + b * b + c * c);
        projected_points(i, 0) = arma::as_scalar(nodos(i, 0)) + this->a * t;
        projected_points(i, 1) = arma::as_scalar(nodos(i, 1)) + this->b * t;
        projected_points(i, 2) = arma::as_scalar(nodos(i, 2)) + this->c * t;
        arma::mat vector = nodos.row(i) - projected_points.row(i); // vector del puntoPro. al nodo
        if (vector(2) > -1e-10)
        {
            // porque la normal apunta hacia arriba y se proyecta en esa direcciomm
            zCoordinates(i) = arma::norm(vector, 2); // esto si esta por encima
        }
        else
        {
            zCoordinates(i) = -arma::norm(vector, 2);
        }
    }
    normales.col(0) = normalPlano(0) * arma::ones(numNodos, 1);
    normales.col(1) = normalPlano(1) * arma::ones(numNodos, 1);
    normales.col(2) = normalPlano(2) * arma::ones(numNodos, 1);
    aRetornar(0, 0) = projected_points;
    aRetornar(0, 2) = normales;
    aRetornar(0, 1) = zCoordinates;
    return aRetornar;
}

// BATHYMETRY DEFINITION
void Bathymetry::ReadPropertiesASCII(FILE *&pFilePointer, std::string inputFilePath)
{

    // Read properties from file
    SeaFloor::ReadPropertiesASCII(pFilePointer);
    char buffer_line[1000];
    std::string header_check;
    // Check if the mesh file exists
    std::string meshFilePath = JoinPath(inputFilePath, meshFileName);
    std::cout << "meshFileName = " << std::endl
              << meshFileName << std::endl;
    FILE *file_pointer = fopen(meshFilePath.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file in which the bathymetry data is\n    ->Dir: " << inputFilePath << std::endl;
        throw IOError(ss.str());
    }
    fscanf(file_pointer, "%d %[^\n]\n", &numPuntosNube, buffer_line);
    std::cout << "Number of points (mesh) = " << std::endl
              << numPuntosNube << std::endl;
    pointMatrix = arma::zeros(numPuntosNube, 3);

    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), file_pointer);
        header_check = buffer_line;
        if (header_check.substr(0, 3).compare("///"))
        {
            std::stringstream ss;
            ss << "Error while parsing file HINT SEAFLOOR ID: " << this->GetId() << " - Please check that each type of seafloor has its correct number of inputs.";
            throw IOError(ss.str());
        }
    }
    // no lee bin esto y no se por que
    // voy a cambiar pFilePointer por file_pointer OK, ESTA BIEN
    // bufferline da un salto, pero se están acumulando dos saltos. VOY A QUITARLO
    for (int i = 0; i < numPuntosNube; i++)
    {
        fscanf(file_pointer, "%lf %lf %lf %\n", &pointMatrix(i, 0), &pointMatrix(i, 1), &pointMatrix(i, 2), buffer_line);
    }
    // LECTURA DE TRIANGULOS
    // se salta 3 lineas de nuevo
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), file_pointer);
        header_check = buffer_line;
        if (header_check.substr(0, 3).compare("///"))
        {
            std::stringstream ss;
            ss << "Error. SeaFloor ID:" << this->GetId() << " - Please check that each type of SeaFloor has its correct number of inputs.";
            throw IOError(ss.str());
        }
    }
    fscanf(file_pointer, "%d %[^\n]\n", &numTriangulos, buffer_line);
    std::cout << "Number of triangles = " << std::endl
              << numTriangulos << std::endl;
    // se salta 3 lineas de nuevo
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), file_pointer);
        header_check = buffer_line;
        if (header_check.substr(0, 3).compare("///"))
        {
            std::stringstream ss;
            ss << "Error. SeaFlor ID: " << this->GetId() << " - Please check that each type of SeaFloor has its correct number of inputs.";
            throw IOError(ss.str());
        }
    }
    triangleMatrix.zeros(numTriangulos, 3);
    for (int i = 0; i < numTriangulos; i++)
    {
        fscanf(file_pointer, "%d %d %d %\n", &triangleMatrix(i, 0), &triangleMatrix(i, 1), &triangleMatrix(i, 2), buffer_line);
    }
    // se salta 3 lineas de nuevo
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), file_pointer);
        header_check = buffer_line;
        if (header_check.substr(0, 3).compare("///"))
        {
            std::stringstream ss;
            ss << "Error. SeaFloor ID:" << this->GetId() << " - Please check that each type of SeaFloor has its correct number of inputs.";
            throw IOError(ss.str());
        }
    }
    fscanf(file_pointer, "%d %[^\n]\n", &flagBarycenter, buffer_line);
    fclose(file_pointer);
}
int Bathymetry::GetType(void)
{
    return this->seaFloorType;
}

void Bathymetry::getVertexNormals(void)
{
    // se usa las matrices leidas pointMatrix y triangleMatrix
    // se itera con cada triangulo y se sacan las coordenadas de sus vertices
    // OJO si triangleMatrix(i)= 1 2 3, las filas donde hay que buscar las coordenadas son 1-1 2-1 3-1 por venir de Matlab
    arma::mat normaTriangulo = arma::ones(numTriangulos, 3);
    barycenter = arma::ones(numTriangulos, 3);
    for (int i = 0; i < numTriangulos; i++)
    {
        arma::mat V0 = arma::ones(1, 3);
        arma::mat V1 = arma::ones(1, 3);
        arma::mat V2 = arma::ones(1, 3);
        V0 = vertexCoordinates(i).row(0);
        V1 = vertexCoordinates(i).row(1);
        V2 = vertexCoordinates(i).row(2);
        // coordenadas del baricentro
        barycenter.row(i) = (V0 + V1 + V2) / 3.0; // media de las coordenadas de los vertices
        // lados
        arma::mat Lado1 = V1 - V0;
        arma::mat Lado2 = V2 - V0;
        // norma del triangulo
        arma::mat prodCross = arma::cross(Lado1, Lado2);
        normaTriangulo.row(i) = prodCross / arma::norm(prodCross, 2);
        // baricentro del triangulo
    }
    vertexNormals = arma::zeros(numPuntosNube, 3);
    for (int i = 0; i < numPuntosNube; i++)
    {
        arma::mat sumaNormal = arma::zeros(1, 3);
        int numTrianguloVertice = 0;
        for (int j = 0; j < numTriangulos; j++)
        {
            if ((triangleMatrix(j, 0) == (i + 1)) || (triangleMatrix(j, 1) == (i + 1)) || (triangleMatrix(j, 2) == (i + 1)))
            {
                sumaNormal = sumaNormal + normaTriangulo.row(j);
                numTrianguloVertice++;
            }
        }
        vertexNormals.row(i) = sumaNormal / numTrianguloVertice;
        vertexNormals.row(i) = vertexNormals.row(i) / arma::norm(vertexNormals.row(i), 2);
    }
}
void Bathymetry::getProjectionMatrix(void)
{
    // inicializacion
    projectionMatrix = arma::field<arma::mat>(numTriangulos, 7);
    changeFrameMatrix = arma::field<arma::mat>(numTriangulos, 4);
    normalsTriangle = arma::field<arma::mat>(numTriangulos, 3);
    for (int k = 0; k < numTriangulos; k++)
    {
        arma::mat V0 = arma::ones(1, 3);
        arma::mat V1 = arma::ones(1, 3);
        arma::mat V2 = arma::ones(1, 3);
        V0 = vertexCoordinates(k).row(0);
        V1 = vertexCoordinates(k).row(1);
        V2 = vertexCoordinates(k).row(2);
        arma::mat bigMatrix = triangleChangeFrame(V0, V1, V2);
        arma::mat M = bigMatrix.submat(0, 0, 3, 3);    // de 3d a 2d
        arma::mat invM = bigMatrix.submat(4, 0, 7, 3); // de 2d a 3d
        changeFrameMatrix(k, 0) = M;
        changeFrameMatrix(k, 1) = invM;
        arma::mat V0n = bigMatrix.submat(8, 1, 8, 3);   // nuevas coordenadas del punto V0
        arma::mat V1n = bigMatrix.submat(9, 1, 9, 3);   // nuevas coordenadas del punto V0
        arma::mat V2n = bigMatrix.submat(10, 1, 10, 3); // nuevas coordenadas del punto V0
        changeFrameMatrix(k, 2) = V1n;
        changeFrameMatrix(k, 3) = V2n;
        // ahora hay que ver las normales de los vertices en este nuevo sistema
        arma::mat n0 = arma::zeros(4, 1);
        arma::mat n1 = arma::zeros(4, 1);
        arma::mat n2 = arma::zeros(4, 1);
        n0(0) = 0;
        n1(0) = 0;
        n2(0) = 0;
        for (int m = 1; m < 4; m++)
        {
            n0(m) = vertexNormals(triangleMatrix(k, 0) - 1, m - 1);
            n1(m) = vertexNormals(triangleMatrix(k, 1) - 1, m - 1);
            n2(m) = vertexNormals(triangleMatrix(k, 2) - 1, m - 1);
        }
        // ya tengo los vectores 4x4
        arma::mat n0n = M * n0;
        arma::mat n1n = M * n1;
        arma::mat n2n = M * n2;
        // los guardo
        normalsTriangle(k, 0) = n0n;
        normalsTriangle(k, 1) = n1n;
        normalsTriangle(k, 2) = n2n;
        // betas para cada triangulo
        double b00 = n0n(1) / n0n(3);
        double b01 = n0n(2) / n0n(3);
        double b10 = n1n(1) / n1n(3);
        double b11 = n1n(2) / n1n(3);
        double b20 = n2n(1) / n2n(3);
        double b21 = n2n(2) / n2n(3);
        // matrices
        arma::mat U = arma::ones(2, 2);
        U(0, 0) = b10 - b00;
        U(0, 1) = V1n(0);
        U(1, 0) = b11 - b01;
        U(1, 1) = 0;
        arma::mat V = arma::ones(2, 2);
        V(0, 0) = b20 - b00;
        V(0, 1) = V2n(0);
        V(1, 0) = b21 - b01;
        V(1, 1) = V2n(1);
        arma::mat W = arma::ones(2, 4);
        W(0, 0) = 1;
        W(0, 1) = 0;
        W(0, 2) = -b00;
        W(0, 3) = 0;
        W(1, 0) = 0;
        W(1, 1) = 1;
        W(1, 2) = -b01;
        W(1, 3) = 0;
        arma::mat UU = arma::strans(U) * U;
        arma::mat UV = arma::strans(U) * V;
        arma::mat VV = arma::strans(V) * V;
        arma::mat UW = arma::strans(U) * W;
        arma::mat VW = arma::strans(V) * W;
        projectionMatrix(k, 0) = UU;
        projectionMatrix(k, 1) = UV;
        projectionMatrix(k, 2) = VV;
        projectionMatrix(k, 3) = UW;
        projectionMatrix(k, 4) = VW;
    }
}
arma::uvec Bathymetry::closerTriangles(arma::mat point)
{
    arma::mat distanceMatrix = arma::ones(numTriangulos, 1);
    // std::cout << "hola pepsi" << std::endl;
    for (int k = 0; k < numTriangulos; k++)
    {
        // encuentra las distsncias del baricentro al punto
        distanceMatrix(k) = arma::norm(barycenter.row(k) - point, 2); // distancia del baricentro al punto en XY porque hemos cambiado de sist.ref. Punto tiene tamaño 4x1
    }
    arma::uvec indicesOrdenados = arma::sort_index(distanceMatrix);
    // estaria bien que devolviera no solo la distancia sino el punto en el sist. de ref 2D

    return indicesOrdenados;
}

arma::field<arma::mat> Bathymetry::projectPoints(arma::mat nodos)
{
    // se inicializa lo que se devuelve
    // el primer elemeto seran los puntos proyectados
    // el segundo elemento seran las z en 2d, que es un vector
    // el tercer elemento seran las normales en cada punto
    arma::field<arma::mat> aRetornar(1, 3);
    // INICIALIZACION DE NODOS
    int numNodos = nodos.n_rows;
    arma::mat estaProyectado = arma::zeros(numNodos, 1);   // valdra 0 si no ha sido aun proyectado
    arma::mat projected_points = arma::zeros(numNodos, 3); // aqui se devolveran los puntos ya proyectados
    // INICIALIZACION DE Z
    arma::mat zCoordinates = arma::zeros(numNodos, 1);
    // INICIALIZACION DE LAS NORMALES  RETORNAR
    arma::mat normales = arma::zeros(numNodos, 3);
    double tol = 1e-10;
    if (flagBarycenter == 1)
    {
        int k;
        // std::cout << "hola cocacola" << std::endl;
        for (int i = 0; i < numNodos; i++)
        {
            arma::uvec triangulosCerca = closerTriangles(nodos.row(i));
            for (int j = 0; j < numTriangulos; j++)
            {
                if (estaProyectado(i) == 0)
                {
                    k = triangulosCerca(j);
                    arma::mat punto = arma::zeros(4, 1);
                    punto(0) = 1;
                    punto(1) = nodos(i, 0);
                    punto(2) = nodos(i, 1);
                    punto(3) = nodos(i, 2);
                    arma::mat puntoGirado = changeFrameMatrix(k, 0) * punto;
                    // la z ya esta en el buen sistema de referencia
                    double z = puntoGirado(3);
                    zCoordinates(i) = z; // almacena la coordenada z
                    arma::mat puntobien = arma::zeros(4, 1);
                    puntobien(0) = puntoGirado(1);
                    puntobien(1) = puntoGirado(2);
                    puntobien(2) = z;
                    puntobien(3) = 1;
                    arma::mat zvec = arma::ones(2, 1);
                    zvec(0) = z;
                    zvec(1) = 1;
                    arma::mat uu = arma::strans(zvec) * projectionMatrix(k, 0) * zvec;
                    arma::mat vv = arma::strans(zvec) * projectionMatrix(k, 2) * zvec;
                    arma::mat uv = arma::strans(zvec) * projectionMatrix(k, 1) * zvec;
                    arma::mat uw = arma::strans(zvec) * projectionMatrix(k, 3) * puntobien;
                    arma::mat vw = arma::strans(zvec) * projectionMatrix(k, 4) * puntobien;
                    // coordenadas baricentricas
                    double s = arma::as_scalar(uv * vw - vv * uw) / arma::as_scalar(uv * uv - uu * vv);
                    double t = arma::as_scalar(uv * uw - uu * vw) / arma::as_scalar(uv * uv - uu * vv);
                    // dan valores razonables
                    if ((t >= -tol) && (s >= -tol) && (s + t <= (1.0 + tol)))
                    {
                        estaProyectado(i) = 1;
                        arma::mat puntoProyectar = arma::zeros(1, 3);
                        puntoProyectar = changeFrameMatrix(k, 2) * s + changeFrameMatrix(k, 3) * t;
                        arma::mat puntoProyectarBien = arma::zeros(4, 1);
                        puntoProyectarBien(0) = 1;
                        puntoProyectarBien(1) = puntoProyectar(0);
                        puntoProyectarBien(2) = puntoProyectar(1);
                        puntoProyectarBien(3) = puntoProyectar(2);
                        arma::mat puntoSol = changeFrameMatrix(k, 1) * puntoProyectarBien;
                        projected_points(i, 0) = puntoSol(1);
                        projected_points(i, 1) = puntoSol(2);
                        projected_points(i, 2) = puntoSol(3);
                        // calculo de las normales
                        arma::mat normalNoBien = arma::zeros(1, 3);
                        // normals triangle te las devuelve bien
                        normalNoBien = normalsTriangle(k, 0) + (normalsTriangle(k, 1) - normalsTriangle(k, 0)) * s + (normalsTriangle(k, 2) - normalsTriangle(k, 0)) * t;
                        arma::mat normal3D = changeFrameMatrix(k, 1) * normalNoBien;
                        normales(i, 0) = normal3D(1);
                        normales(i, 1) = normal3D(2);
                        normales(i, 2) = normal3D(3);
                    }
                }
            }
        }
    } // el de flagBarycenter

    if (flagBarycenter == 0)
    {
        for (int k = 0; k < numTriangulos; k++)
        {
            for (int i = 0; i < numNodos; i++)
            {
                if (estaProyectado(i) == 0)
                {
                    arma::mat punto = arma::zeros(4, 1);
                    punto(0) = 1;
                    punto(1) = nodos(i, 0);
                    punto(2) = nodos(i, 1);
                    punto(3) = nodos(i, 2);
                    arma::mat puntoGirado = changeFrameMatrix(k, 0) * punto;
                    // la z ya esta en el buen sistema de referencia
                    double z = puntoGirado(3);
                    zCoordinates(i) = z; // almacena la coordenada z
                    arma::mat puntobien = arma::zeros(4, 1);
                    puntobien(0) = puntoGirado(1);
                    puntobien(1) = puntoGirado(2);
                    puntobien(2) = z;
                    puntobien(3) = 1;
                    arma::mat zvec = arma::ones(2, 1);
                    zvec(0) = z;
                    zvec(1) = 1;
                    arma::mat uu = arma::strans(zvec) * projectionMatrix(k, 0) * zvec;
                    arma::mat vv = arma::strans(zvec) * projectionMatrix(k, 2) * zvec;
                    arma::mat uv = arma::strans(zvec) * projectionMatrix(k, 1) * zvec;
                    arma::mat uw = arma::strans(zvec) * projectionMatrix(k, 3) * puntobien;
                    arma::mat vw = arma::strans(zvec) * projectionMatrix(k, 4) * puntobien;
                    // coordenadas baricentricas
                    double s = arma::as_scalar(uv * vw - vv * uw) / arma::as_scalar(uv * uv - uu * vv);
                    double t = arma::as_scalar(uv * uw - uu * vw) / arma::as_scalar(uv * uv - uu * vv);
                    // dan valores razonables
                    if ((t >= -tol) && (s >= -tol) && (s + t <= (1.0 + tol)))
                    {
                        estaProyectado(i) = 1;
                        arma::mat puntoProyectar = arma::zeros(1, 3);
                        puntoProyectar = changeFrameMatrix(k, 2) * s + changeFrameMatrix(k, 3) * t;
                        arma::mat puntoProyectarBien = arma::zeros(4, 1);
                        puntoProyectarBien(0) = 1;
                        puntoProyectarBien(1) = puntoProyectar(0);
                        puntoProyectarBien(2) = puntoProyectar(1);
                        puntoProyectarBien(3) = puntoProyectar(2);
                        arma::mat puntoSol = changeFrameMatrix(k, 1) * puntoProyectarBien;
                        projected_points(i, 0) = puntoSol(1);
                        projected_points(i, 1) = puntoSol(2);
                        projected_points(i, 2) = puntoSol(3);
                        // calculo de las normales
                        arma::mat normalNoBien = arma::zeros(1, 3);
                        // normals triangle te las devuelve bien
                        normalNoBien = normalsTriangle(k, 0) + (normalsTriangle(k, 1) - normalsTriangle(k, 0)) * s + (normalsTriangle(k, 2) - normalsTriangle(k, 0)) * t;
                        arma::mat normal3D = changeFrameMatrix(k, 1) * normalNoBien;
                        normales(i, 0) = normal3D(1);
                        normales(i, 1) = normal3D(2);
                        normales(i, 2) = normal3D(3);
                    }
                }
            }
        }
    }

    aRetornar(0, 0) = projected_points;
    aRetornar(0, 1) = zCoordinates;
    aRetornar(0, 2) = normales;
    int hayError = 0; // para lanzar un aviso por lo del baricentro mas cercano
    for (int i = 0; i < numNodos; i++)
    {
        if (projected_points(i, 0) == 0 && projected_points(i, 1) == 0 && projected_points(i, 2) == 0)
        {
            std::cout << "     WARNING: Probably, node" << i + 1 << " could not be projected. Check the floor size. " << std::endl;
            hayError = 1;
        }
    }
    if (hayError = 1 && flagBarycenter == 1)
    {
        std::stringstream ss;
        ss << "If it is not a floor size problem, try not to use the closest barycenter instead. \n";
        throw ValueError(ss.str());
    }
    return aRetornar;
}

arma::mat Bathymetry::vertexCoordinates(int triangle)
{
    // triangle es la fila de la matriz de triangulos de la que queremos saber los vertices
    // se leen los indices que son las filas de la matriz de puntos
    int indexV0 = triangleMatrix(triangle, 0) - 1;
    int indexV1 = triangleMatrix(triangle, 1) - 1;
    int indexV2 = triangleMatrix(triangle, 2) - 1;
    // coordenadas de los vertices del triangulo
    arma::mat V0 = pointMatrix.row(indexV0);
    arma::mat V1 = pointMatrix.row(indexV1);
    arma::mat V2 = pointMatrix.row(indexV2);

    // devuelve los vertices en un array 3x3
    arma::mat unionVertices = arma::join_cols(V0, V1);
    arma::mat coordenadas = arma::join_cols(unionVertices, V2);
    return coordenadas;
}