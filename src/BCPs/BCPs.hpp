
#ifndef BCP_FLAG
#define BCP_FLAG
#include <armadillo>
#include <string>

class BCP {
public:
	//Atributos comunes a todos los BCPs
	double tBCP;
	int typeBCP;
	int nBCP; // Indice identificador del punto de condicion de contorno
	int nLinesBCP; // Numero de lineas que confluyen en el punto
	int * BCPLineIndex; // Array con los indices identificadores de las lineas que confluyen en el punto
	int * BCPLineNode; // Array de flags que, para cada linea ii que confluye al punto, indica si la linea confluye al nodo 1 (BCPLineNode[ii]=1) o al nodo N (BCPLineNode[ii]=2)
	arma::mat pos; // Posicion del punto
	arma::mat vel; // Velocidad del punto
	arma::mat acc; // Aceleracion del punto

	//Necesario para Fairlead
	std::string fileName;
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
	arma::mat posL = arma::zeros(3,1); // Posicion del punto en el cuerpo en local en 3 dofs
	arma::mat posG = arma::zeros(3,1); // Posicion del punto en el cuerpo en global en 3 dofs
	arma::mat RotMat = arma::zeros(3,3); // Matriz de rotacion del cuerpo
	arma::mat posG_BCP = arma::zeros(6,1); // Posicion del punto en global en 6 dofs
	arma::mat velG_BCP = arma::zeros(6,1); // Velocidad del punto en global en 6 dofs
	arma::mat accG_BCP = arma::zeros(6,1); // Aceleración del punto en global en 6 dofs
	arma::mat ForceBCP = arma::zeros(6,1); // Fuerzas y momentos que actuan sobre el BCP en 6 dofs

	//Metodos
	void set_nBCP(int n){
		nBCP=n;
		pos = arma::zeros(3,1); // Posicion del punto
		vel = arma::zeros(3,1); // Velocidad del punto
		acc = arma::zeros(3,1); // Aceleracion del punto
	}
	void leer_datosBCPs(void);
	virtual void getValues(double t) =0;
};

class AnchorBCP: public BCP {
public:
	void getValues(double t);
};


class FairleadBCP: public BCP {
public:
	void getValues(double t);
};

class JointBCP: public BCP {
public:
	void getValues(double t);
};

class BodyBCP: public BCP {
public:
	void getValues(double t);
};



#endif