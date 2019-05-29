#include <armadillo>
#include <string>
#include "../BCPs/BCPs.hpp"

class Body {
public:

	int nBody; // Indice identificador del cuerpo
	int nDOFs; // Numero de grados de libertad que se consideran en el cuerpo
	int* DOFs; // Array de grados de libertad que se consideran en el cuerpo

	int nBCPs; // Numero de BCPs en el cuerpo
	int* index_BCPs;
	BCP* BodyBCPs; // Array de pointers a los puntos de condicion de contorno

	arma::mat pos = arma::zeros(6,1); // Posicion del cuerpo
	arma::mat vel = arma::zeros(6,1); // Velocidad del cuerpo
	arma::mat acc = arma::zeros(6,1); // Aceleracion del punto

	arma::mat inertia = arma::zeros(6,6); // matriz de inercia del cuerpo

	arma::mat RotMat = arma::zeros(3,3); // Matriz de rotación

	void set_nBody(int n){nBody=n;} // Inicializa un objeto de clase cuerpo dandole el indice
	void leer_datosBody(void); // Leer datos de los cuerpos
	void getRotMat(void); // Obten la matriz de rotación y pasasela a los BCPs, junto con la posicion, velocidad y aceleración
};