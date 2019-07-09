
#include <iostream>
#include <sstream>
#include <cstdio>
#include <armadillo>
#include "MathTools.hpp"



arma::mat arange(double a, double b, double stepSize)
{
	// Calculate number of points
	int numPoints = (b-a)/stepSize + 1;
	
	// Allocate solution vector
	arma::mat vec = arma::zeros(1, numPoints);
	
	// Fill solution vector
	for (int i=0; i<numPoints; i++)
	{
		vec[i] = a + i*stepSize;
	}
	
	return vec;
}


arma::mat interp1(arma::mat x, arma::mat y, arma::mat xi)
{
	// Declare and allocate solution vector
	arma::mat y1 = arma::zeros(1, xi.n_cols);
	
	// Loop over coordinates to interpolate
	int m = 0;
	int index = 0;
	double dx = 0;
	for (int i=0; i<xi.n_cols; i++)
	{
		index = -1;
		// Look interval to perform interpolation
		for (int j=m; j<x.n_cols-1; j++)
		{
			if ((xi[i] >= x[j]) && (xi[i] <= x[j+1]))
			{
				index = j;
				break;
			}
		}
		
		if (index == -1)
		{
			std::stringstream ss;
			ss << "Not possible to find an interval for the abcissa value: " << xi[i];
			perror(ss.str().c_str());
		}
		dx = x[index+1] - x[index];
		y1[i] = y[index+1]*(xi[i]-x[index])/dx + y[index]*(x[index+1]-xi[i])/dx;
		m = index;
	}
	
	return y1;
}


arma::mat linspace(double a, double b, int numPoints)
{
	// Allocate solution vector
	arma::mat vec = arma::zeros(1, numPoints);

	// Calculate step size
	double stepSize = (b-a)/(numPoints-1);
	
	// Fill vector
	for (int i=0; i<numPoints; i++)
	{
		vec[i] = a + i*stepSize;
	}
	
	return vec;
}


double trapz(arma::mat y, double h)
{
	double int_value = 0.0;
	for (int i=0; i<y.n_cols-1; i++)
	{
		int_value += (y[i+1]+y[i])/2.0;
	}
	int_value *= h;
	
	return int_value;
}


double trapzi(arma::mat t, arma::mat y)
{
	double int_value = 0.0;
	double h;
	for (int i=0; i<y.n_cols-1; i++)
	{
		h = t[i+1] - t[i];
		int_value += h*(y[i+1]+y[i])/2.0;
	}
	
	return int_value;
}