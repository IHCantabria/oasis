
#ifndef bcp_hpp__
#define bcp_hpp__
#include <armadillo>
#include <string>
#include <cstdio>

class BCP
{
private:
	int id;
	int typeBcp=0;
public:
	//Atributos comunes a todos los BCPs
	double tBCP;
	int nLinesBCP; // Numero de lineas que confluyen en el punto
	int * BCPLineIndex; // Array con los indices identificadores de las lineas que confluyen en el punto
	int * BCPLineNode; // Array de flags que, para cada linea ii que confluye al punto, indica si la linea confluye al nodo 1 (BCPLineNode[ii]=1) o al nodo N (BCPLineNode[ii]=2)
	arma::mat pos; // Posicion del punto
	arma::mat vel; // Velocidad del punto
	arma::mat acc; // Aceleracion del punto

	//Necesario para Fairlead
	std::string actuatorFileName;
	arma::mat tF, posF, velF, accF;
	arma::mat pos0; // Posicion inicial del punto
	int ni = 0;
	double dt;

	//Necesario para Joint
	int iLJ = 0;
	arma::mat posLines;
	arma::mat velLines;
	arma::mat accLines;

	//Necesario para Body
	//arma::mat posL = arma::zeros(3,1); // Posicion del punto en el cuerpo en local en 3 dofs
	arma::mat posLinLocal = arma::zeros(3,1); // Posicion del punto en el cuerpo en local en 3 dofs
	//arma::mat posG = arma::zeros(3,1); // Posicion del punto en el cuerpo en global en 3 dofs
	arma::mat posLinGlobal = arma::zeros(3,1); // Posicion del punto en el cuerpo en global en 3 dofs
	arma::mat rotMat = arma::zeros(3,3); // Matriz de rotacion del cuerpo
	arma::mat posGlobal = arma::zeros(6,1); // Posicion del punto en global en 6 dofs
	arma::mat velGlobal = arma::zeros(6,1); // Velocidad del punto en global en 6 dofs
	arma::mat accGlobal = arma::zeros(6,1); // Aceleración del punto en global en 6 dofs
	arma::mat forceBcp = arma::zeros(6,1); // Fuerzas y momentos que actuan sobre el BCP en 6 dofs

	// Methods	
	BCP(int incId);
	virtual int GetType(void);
	virtual void GetValues(double t) =0;
	void ReadPropertiesASCII(FILE* pFilePointer);
};

class AnchorBCP: public BCP
{
private:
	int typeBcp=1;
public:
	AnchorBCP(int incId);
	virtual int GetType(void);
	void GetValues(double t);
};

class BodyBCP: public BCP
{
private:
	int typeBcp=4;
public:
	BodyBCP(int incId);
	virtual int GetType(void);
	void GetValues(double t);
};

class FairleadBCP: public BCP
{
private:
	int typeBcp=2;
public:
	FairleadBCP(int incId);
	virtual int GetType(void);
	void GetValues(double t);
};

class JointBCP: public BCP
{
private:
	int typeBcp=3;
public:
	JointBCP(int incId);
	virtual int GetType(void);
	void GetValues(double t);
};

#endif // bcp_hpp__