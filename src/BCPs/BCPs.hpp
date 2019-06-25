
#ifndef bcpdef_hpp__
#define bcpdef_hpp__
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
	int numLinesBcp; // Numero de lineas que confluyen en el punto
	int* pBcpLineIndex; // Array con los indices identificadores de las lineas que confluyen en el punto
	int* pBcpLineNode; // Array de flags que, para cada linea ii que confluye al punto, indica si la linea confluye al nodo 1 (BCPLineNode[ii]=1) o al nodo N (BCPLineNode[ii]=2)
	arma::mat pos = arma::zeros(3, 1); // Posicion del punto
	arma::mat vel = arma::zeros(3, 1); // Velocidad del punto
	arma::mat acc = arma::zeros(3, 1); // Aceleracion del punto

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
	arma::mat posWrtCdgLocal = arma::zeros(3,1); // Brazo cdg-bcp en local
	arma::mat posWrtCdgGlobal = arma::zeros(3,1); // Brazo cdg-bcp en global
	arma::mat rotMat = arma::zeros(3,3); // Matriz de rotacion del cuerpo
	arma::mat posG_BCP = arma::zeros(6,1); // Posicion del punto en global en 6 dofs
	arma::mat velG_BCP = arma::zeros(6,1); // Velocidad del punto en global en 6 dofs
	arma::mat accG_BCP = arma::zeros(6,1); // Aceleración del punto en global en 6 dofs
	arma::mat forceBcp = arma::zeros(6,1); // Fuerzas y momentos que actuan sobre el BCP en 6 dofs

	// Methods	
	BCP(int incId);
	virtual int GetType(void);
	virtual void GetValues(double t) = 0;
	virtual void Initialize();
	void ReadPropertiesASCII(FILE*& pFilePointer);
	virtual void UpdateBoundary();
};

class AnchorBCP: public BCP
{
private:
	int typeBcp=1;
public:
	AnchorBCP(int incId): BCP(incId) {};
	void GetValues(double t);
};

class BodyBCP: public BCP
{
private:
	int typeBcp=4;
public:
	BodyBCP(int incId): BCP(incId) {};
	void GetValues(double t);
};

class FairleadBCP: public BCP
{
private:
	int typeBcp=2;
public:
	FairleadBCP(int incId): BCP(incId) {};
	void GetValues(double t);
	void Initialize(std::string folder_path);
};

class JointBCP: public BCP
{
private:
	int typeBcp=3;
public:
	JointBCP(int incId): BCP(incId) {};
	void GetValues(double t);
};

#endif // bcp_hpp__