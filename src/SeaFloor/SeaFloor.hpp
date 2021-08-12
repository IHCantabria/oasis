#ifndef seafloordef_hpp__
#define seafloordef_hpp__

#include <armadillo>
#include <string>

class SeaFloor {
private:
	int id; //esto sirve para luego, y para el constructor
	int seaFloorType=0; //valor no coherente que luego cambiara

 //Atributes
public:
 	std::string meshFileName; //para leer el mallado en batimetria variable
	//para batimetria variable
	int numPuntosNube;
	int numTriangulos;
	double fondo;

 
 	//para el plano inclinado
 	arma::mat point1 = arma::zeros(3,1);
	arma::mat point2 = arma::zeros(3,1);
 	arma::mat point3= arma::zeros(3,1);


 	//Methods
 	SeaFloor(int incId);
 	int GetId(void);
 	virtual int GetType(void) = 0;
	virtual arma::field<arma::mat> projectPoints(arma::mat nodos) = 0;
 	//los de batimetria que tienen que estar definidos
 	void ReadPropertiesASCII(FILE* &pFilePointer);


 

};


class Bathymetry: public SeaFloor
{
private:
	int seaFloorType=3;
	//se define un metodo privado solo utilizado dentro de la propia clase
	arma::mat vertexCoordinates(int triangle);

public:
	//atributos
	int numPuntosNube;
	int numTriangulos;
	arma::mat pointMatrix;
	arma::umat triangleMatrix;
	arma::mat vertexNormals;
	arma::field<arma::mat> projectionMatrix;
	arma::field<arma::mat> changeFrameMatrix;
	arma::field<arma::mat> normalsTriangle;
	//Metodos
	Bathymetry(int incId): SeaFloor(incId) {};
	void ReadPropertiesASCII(FILE* &pFilePointer, std::string inputFilePath);
	void getVertexNormals(void);
	void getProjectionMatrix(void);
	arma::field<arma::mat> projectPoints(arma::mat nodos);
	int GetType(void);
};

class Inclined: public SeaFloor
{
private:
	int seaFloorType=2;


public:
	Inclined(int incId): SeaFloor(incId) {};
	arma::mat getNormal(arma::mat p1, arma::mat p2, arma::mat p3);
	arma::field<arma::mat> projectPoints(arma::mat nodos);
	int GetType(void);

};

class Flat: public SeaFloor
{
private:
	int seaFloorType=1;
public:
	Flat(int incId): SeaFloor(incId) {};
	arma::field<arma::mat> projectPoints(arma::mat nodos);
	int GetType(void);
	
	
};


#endif