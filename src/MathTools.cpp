
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


/*
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
*/


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


arma::mat mod(arma::mat a, double x)
{
	return a - arma::floor(a/x)*x;
}

arma::mat interp1(arma::mat x, arma::mat y, arma::mat xi)
{
	// Check that x and xi are vectors, reshape to columns
	if((x.n_cols>1)&&(x.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: x must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	if((xi.n_cols>1)&&(xi.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: xi must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	x = arma::reshape(x,x.n_elem,1);
	xi = arma::reshape(xi,xi.n_elem,1);

	// Check that the number of columns in y is the same as in x
	if(y.n_rows!=x.n_rows)
	{
		std::stringstream ss;
	    ss << "Error in interp1: number of rows in y must be equal to x length. \n";
	    throw ValueError(ss.str());
	}

	// Initiallize the output matrix
	arma::mat yi = arma::zeros(xi.n_rows,y.n_cols);
	
	arma::mat temp_input, temp_output;

	for(int ii=0; ii<y.n_cols; ii++)
	{
		temp_input = y(arma::span::all,arma::span(ii));
		arma::interp1(x,temp_input,xi,temp_output,"linear",0);
		yi(arma::span::all,arma::span(ii)) = temp_output;
	}

	return yi;
}

arma::cube interp1(arma::mat x, arma::cube y, arma::mat xi)
{
	// Check that x and xi are vectors, reshape to columns
	if((x.n_cols>1)&&(x.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: x must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	if((xi.n_cols>1)&&(xi.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: xi must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	x = arma::reshape(x,x.n_elem,1);
	xi = arma::reshape(xi,xi.n_elem,1);

	// Check that the number of columns in y is the same as in x
	if(y.n_rows!=x.n_rows)
	{
		std::stringstream ss;
	    ss << "Error in interp1: number of rows in y must be equal to x length. \n";
	    throw ValueError(ss.str());
	}

	// Initiallize the output matrix
	arma::cube yi = arma::zeros(xi.n_rows,y.n_cols,y.n_slices);

	arma::mat temp_input, temp_output;

	for(int ii=0; ii<y.n_cols; ii++)
	{
		for(int jj=0; jj<y.n_slices; jj++)
		{
			temp_input = y(arma::span::all,arma::span(ii),arma::span(jj));
			arma::interp1(x,temp_input,xi,temp_output,"linear",0);
			yi(arma::span::all,arma::span(ii),arma::span(jj)) = temp_output;
		}
	}

	return yi;
}


arma::cube interp2(arma::mat x, arma::mat y, arma::cube z, arma::mat xi, arma::mat yi)
{
	// Check that x, y, xi and yi are vectors, reshape to columns
	if((x.n_cols>1)&&(x.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: x must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	if((xi.n_cols>1)&&(xi.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: xi must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	if((y.n_cols>1)&&(y.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: y must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	if((yi.n_cols>1)&&(yi.n_rows>1))
	{
		std::stringstream ss;
	    ss << "Error in interp1: yi must be a column or row vector. \n";
	    throw ValueError(ss.str());
	}
	x = arma::reshape(x,x.n_elem,1);
	xi = arma::reshape(xi,xi.n_elem,1);
	y = arma::reshape(x,x.n_elem,1);
	yi = arma::reshape(xi,xi.n_elem,1);


	// Check that the number of rows in z is the same as in x,
	// and that the number of columns in z is the same as in y.
	if((z.n_rows!=x.n_rows)||(z.n_cols!=y.n_rows))
	{
		std::stringstream ss;
	    ss << "Error in interp2: number of rows in z must be equal to x length, and the number of columns in z must be equal to y length. \n";
	    throw ValueError(ss.str());
	}

	// Initiallize the output matrix
	arma::cube zi = arma::zeros(xi.n_rows,yi.n_rows,z.n_slices);

	arma::mat temp_input, temp_output;

	for(int ii=0; ii<z.n_slices; ii++)
	{
		temp_input = z(arma::span::all,arma::span::all,arma::span(ii));
		arma::interp2(x,y,temp_input,xi,yi,temp_output,"linear",0);
		zi(arma::span::all,arma::span::all,arma::span(ii)) = temp_output;
	}

	return zi;
}

arma::cube permute(arma::cube x, int ind)
{
	arma::cube y;
	int n1 = x.n_rows, n2 = x.n_cols, n3 = x.n_slices;
	switch (ind)
	{
		case 123:
		{
			y = x;
		}
        break;
		case 132:
		{
			y = arma::zeros(n1,n3,n2);
			for (int i1 = 0; i1 < n1; i1++)
                for (int i2 = 0; i2 < n2; i2++)
                    for (int i3 = 0; i3 < n3; i3++)
                        y(i1,i3,i2) = x(i1,i2,i3);
		}
        break;
		case 231:
		{
			y = arma::zeros(n2,n3,n1);
			for (int i1 = 0; i1 < n1; i1++)
                for (int i2 = 0; i2 < n2; i2++)
                    for (int i3 = 0; i3 < n3; i3++)
                        y(i2,i3,i1) = x(i1,i2,i3);
		}
        break;
		case 213:
		{
			y = arma::zeros(n2,n1,n3);
			for (int i1 = 0; i1 < n1; i1++)
                for (int i2 = 0; i2 < n2; i2++)
                    for (int i3 = 0; i3 < n3; i3++)
                        y(i2,i1,i3) = x(i1,i2,i3);
		}
        break;
		case 312:
		{
			y = arma::zeros(n3,n1,n2);
			for (int i1 = 0; i1 < n1; i1++)
                for (int i2 = 0; i2 < n2; i2++)
                    for (int i3 = 0; i3 < n3; i3++)
                        y(i3,i1,i2) = x(i1,i2,i3);
		}
        break;
        break;
		case 321:
		{
			y = arma::zeros(n3,n2,n1);
			for (int i1 = 0; i1 < n1; i1++)
                for (int i2 = 0; i2 < n2; i2++)
                    for (int i3 = 0; i3 < n3; i3++)
                        y(i3,i2,i1) = x(i1,i2,i3);
		}
        break;
        default:
        {
        	std::stringstream ss;
		    ss << "Error in permute: ind must be a integer permutation of 123. \n";
		    throw ValueError(ss.str());
        }
        break;
	}
	return y;
}