
#include <iostream>
#include <sstream>
#include <cstdio>
#include <armadillo>
#include <tuple>
#include "MathTools.hpp"
#include "Exceptions/Exception.hpp"



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
			x.print();
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


std::tuple<arma::mat,arma::mat> upcrossing(arma::mat t, arma::mat u)
{
    
    // Función que, para una serie temporal de altura de ola (t,u), calcula los 
    // tiempos entre upcrossing y upcrossing (periodos T) y las alturas de ola 
    // (H) en tales periodos.
    //
    //
    //  Alvaro Rodriguez Luis
    //  Feb. 2020
    //  IH Cantabria

    // Inicio los outputs
    arma::mat T; arma::mat H;
    // Tolerancia para considerar un valor como cero.
    double tol = 1e-14;
    // Reordeno como vectores columna los inputs
    int nt = t.n_elem; int nu = u.n_elem;
    if(nt!=nu)
   	{
   		std::stringstream ss;
	    ss << "upcrossing: t and u must be the same length. \n";
	    throw ValueError(ss.str());
   	}

    u = arma::reshape(u,nu,1); t = arma::reshape(t,nu,1);

    // Busco los upcrossings y los tiempos en los que se producen
    arma::mat pos = arma::zeros(size(u)); 
    pos.elem(arma::find(u > 0.0)) += 1;
    arma::uvec ind = arma::find(arma::diff(pos)>0);
    arma::mat times = t.elem(ind) - u.elem(ind)%(t.elem(ind+1)-t.elem(ind))/(u.elem(ind+1)-u.elem(ind));
    int N = times.n_elem;

    if (N>1)
    {
    	T = arma::diff(times);
    	H = arma::zeros(size(T));
    	arma::mat range;
    	for(int ii=0; ii<N-1; ii++)
		{
			range = u.rows(ind(ii),ind(ii+1));
			H(ii) = range.max()-range.min();
		}
    }
    else
    {
    	std::stringstream ss;
	    ss << "upcrossing: could not find enough upcrossings. \n";
	    throw ValueError(ss.str());
    }

    return std::make_tuple(T,H);
}