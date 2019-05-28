#include <armadillo>
#include <string>

class Body {
public:

	int nBody; // Indice identificador del cuerpo
	int nDOFs; // Numero de grados de libertad que se consideran en el cuerpo
	int* DOFs; // Array de grados de libertad que se consideran en el cuerpo

	int nBCPs; // Numero de BCPs en el cuerpo
	BCP* BodyBCPs; // Array de pointers a los puntos de condicion de contorno

	arma::mat pos(6,1); // Posicion del cuerpo
	arma::mat vel(6,1); // Velocidad del cuerpo
	arma::mat acc(6,1); // Aceleracion del punto

	double mass; // Masa del cuerpo
	double vol; // volumen del cuerpo
	arma::mat inertia(6,6); // matriz de inercia del cuerpo

	arma::mat RotMat(3,3); // Matriz de rotación

	void leer_datosBody(void); // Leer datos de los cuerpos
	void getRotMat(void); // Obten la matriz de rotación y pasasela a los BCPs, junto con la posicion, velocidad y aceleración
};