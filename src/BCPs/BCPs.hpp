

#include <armadillo>
#include <string>
#include "../SEM_math/SEM_math.hpp"

class BCP {
public:
	int nBCP; // Indice identificador del punto de condicion de contorno
	int nLinesBCP; // Numero de lineas que confluyen en el punto
	int * BCPLineIndex; // Array con los indices identificadores de las lineas que confluyen en el punto
	int * BCPLineNode; // Array de flags que, para cada linea ii que confluye al punto, indica si la linea confluye al nodo 1 (BCPLineNode[ii]=1) o al nodo N (BCPLineNode[ii]=2)
	std::string fileName;
	arma::mat pos = arma::zeros(3,1); // Posicion del punto
	arma::mat vel = arma::zeros(3,1); // Velocidad del punto
	arma::mat acc = arma::zeros(3,1); // Aceleracion del punto
	void set_nBCP(int n){nBCP=n;}
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
	arma::mat posF;
	spline x_spl, y_spl, z_spl;
};