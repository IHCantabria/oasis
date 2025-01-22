
#include <armadillo>
#include <string>
#include "Winchies.hpp"

// Leer inputs
void Winchie::ReadPropertiesASCII(FILE *pFile)
{
	// Declare variables
	char buffer_line[1000];

	// Ignoro las tres primeras lineas, donde pone "New winchie"
	for (int ii = 0; ii < 3; ii++)
	{
		fgets(buffer_line, sizeof(buffer_line), pFile);
	}

	// Leo todo
	fscanf(pFile, "%d %[^\n]\n", &nLine, buffer_line);
	fscanf(pFile, "%d %[^\n]\n", &LineBCP, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &inertia, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &radius, buffer_line);
	fscanf(pFile, "%lf %[^\n]\n", &drag, buffer_line);
}

// Obten las aceleraciones del winchie y el cambio de longitud de las lineas
void Winchie::computeWinchie(void)
{

	double F;

	if (LineBCP == 1)
	{
		F = arma::norm(LineW->ten_1);
	}
	else if (LineBCP == 2)
	{
		F = arma::norm(LineW->ten_N);
	}
	else
	{
		std::cout << std::endl
				  << "ERROR: Options are 1 or 2." << std::endl;
		throw std::exception();
	}

	alpha = (radius * F - tau) / inertia - drag * omega; // Damping harcodeado

	LineW->dL = LineW->dL0 * (LineW->L + radius * theta) / LineW->L;

	// std::cout << "    theta =  " << theta << std::endl << std::endl;
	// std::cout << "    tau =  " << tau << std::endl << std::endl;
}
