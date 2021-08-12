
#include <iostream>
#include <sstream>
#include <cstdio>
#include <armadillo>
#include <tuple>
#include "MathTools.hpp"
#include "Exceptions/Exception.hpp"
#include <math.h>       /* fmod */

arma::mat wrapTo180(arma::mat x)
{
	double a;
	arma::mat y = arma::zeros(x.n_rows,x.n_cols);
	for(int ii=0; ii<x.n_rows; ii++)
	{
		for(int jj=0; jj<x.n_cols; jj++)
		{
			a = arma::as_scalar(x(ii,jj));
			a = fmod(a, 360.0);
			if(a>180.0)
			{
				a = a - 360.0;
			}
			y(ii,jj) = a;
		}
	}
	return y;
}

arma::mat wrapToPi(arma::mat x)
{
	double a;
	arma::mat y = arma::zeros(x.n_rows,x.n_cols);
	for(int ii=0; ii<x.n_rows; ii++)
	{
		for(int jj=0; jj<x.n_cols; jj++)
		{
			a = arma::as_scalar(x(ii,jj));
			a = fmod(a, 2*arma::datum::pi);
			if(a>arma::datum::pi)
			{
				a = a - 2*arma::datum::pi;
			}
			y(ii,jj) = a;
		}
	}
	return y;
}


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
	int num_points = y.n_rows;
	if (y.n_cols > num_points)
	{
		num_points = y.n_cols;
	}
	for (int i=0; i<num_points-1; i++)
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
	y = arma::reshape(y,y.n_elem,1);
	yi = arma::reshape(yi,yi.n_elem,1);

	
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
		arma::interp2(x,y,temp_input.t(),xi,yi,temp_output,"linear",0);
		zi(arma::span::all,arma::span::all,arma::span(ii)) = temp_output.t();
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

arma::uvec comp_ind(int n,arma::uvec ind)
{
	arma::mat a = arma::ones(n,1);
	a.rows(ind) = 0*a.rows(ind);
	arma::uvec ind_out = arma::find(a);
	return ind_out;
}

arma::umat comb_n_k(int n,int k)
{
	int nCombRep, nr, nn, i1, i2;
	nCombRep = pow(n,k);
	arma::umat ind(nCombRep,k); ind.fill(0);
	for(int ik=1; ik<=k; ik=ik+1){
		nr = pow(n,k-ik);
		nn = pow(n,ik-1);
		arma::umat temp(nn,1); temp.fill(1);
		for(int j=1; j<=nr; j=j+1){
			for(int in=1; in<=n; in=in+1){
				i1 = (j-1)*n*nn + (in-1)*nn;
				i2 = (j-1)*n*nn + in*nn - 1;
				ind(arma::span(i1,i2),k-ik) = in*temp - 1;
			}
		}
		delete &temp;
	}
	arma::uvec ind_order;
	for(int ik=1; ik<=k-1; ik=ik+1){
		ind_order = arma::find(ind.col(k-ik)>ind.col(k-ik-1));
		ind = ind.rows(ind_order);
		//arma::umat ind_temp = ind.rows(ind_order);
		//delete &ind;
		//arma::umat ind = ind_temp;
		//delete &ind_temp;
	}
	return ind;
}

double step(double x, double x0, double h0, double x1, double h1)
{
	if (x<x0){
		return h0;
	} else if (x>=x0 && x<=x1){
		return h0 + (h1-h0)*pow((x-x0)/(x1-x0),2)*(3-2*(x-x0)/(x1-x0));
	} else {
		return h1;
	}
}

arma::mat triangleChangeFrame(arma::mat V0, arma::mat V1, arma::mat V2){
	//INPUT vectores 1x3 
	//V0, las coordenadas del vertice en 3D
	//V1, las coordenadas V1 en 3D
	//V20, las coordenadas en 3D

	//OUTPUT es una matriz 11x4 con:
	//de la fila 0 a la 3, matriz M, para cambiar de 3D a 2D en las triangulaciones
	//de la fila 4 a la 7, matriz invM, para cambiar de 2D a 3D en las triangulaciones
	//la fila 8, las coordenadas del vertice V0 en 2D
	//la fila 9, las coordenadas del vertice V1 en 2D
	//la fila 10, las coordenadas del vertice V2 en 2D

	//este metodo rota un triangulo de vertices V0 V1 V2 a otro sobre el plano XY, con V0 en el origen y V1 sobre ejeX

	//traslación
	arma::mat T = arma::eye(4,4);
	T(1,0)= -V0(0);
	T(2,0)= -V0(1);
	T(3,0)= -V0(2);
	arma::mat V00= arma::ones(4,1);
	arma::mat V10= arma::ones(4,1);
	arma::mat V20= arma::ones(4,1);
	V00(1)=V0(0); V00(2)=V0(1); V00(3)=V0(2);
	V10(1)=V1(0); V10(2)=V1(1); V10(3)=V1(2);
	V20(1)=V2(0); V20(2)=V2(1); V20(3)=V2(2);
	arma::mat V01= arma::ones(4,1);
	arma::mat V11= arma::ones(4,1);
	arma::mat V21= arma::ones(4,1);
	V01= T*V00;
	V11=T*V10;
	V21=T*V20;

	//rotacion1_1, ejeY
	double theta= atan2(V11(3), V11(1));
	arma::mat R11 = arma::eye(4,4);
	R11(1,1)=cos(theta);
	R11(1,3) = sin(theta);
	R11(3,1) = -sin(theta);
	R11(3,3)=cos(theta);
	arma::mat V02= arma::ones(4,1);
	arma::mat V12= arma::ones(4,1);
	arma::mat V22= arma::ones(4,1);
	V02= R11*V01;
	V12= R11*V11;
	V22= R11*V21;

	//rotacion1_2, ejeZ
	double theta2= atan2(V12(2),V12(1));
	arma::mat R12= arma::eye(4,4);
	R12(1,1)=cos(theta2);
	R12(1,2)=+sin(theta2);
	R12(2,1)=-sin(theta2);
	R12(2,2)=cos(theta2);
	arma::mat V03= arma::ones(4,1);
	arma::mat V13= arma::ones(4,1);
	arma::mat V23= arma::ones(4,1);

	V03= R12*V02;
	V13= R12*V12;
	V23=R12*V22;

	//rotacion2 eje X
	double theta3= atan2(V23(3),V23(2));
	arma::mat R2=arma::eye(4,4);
	R2(2,2)=cos(theta3);
	R2(2,3)=sin(theta3);
	R2(3,2)=-sin(theta3);
	R2(3,3)=cos(theta3);
	arma::mat V04= arma::ones(4,1);
	arma::mat V14= arma::ones(4,1);
	arma::mat V24= arma::ones(4,1);
	V04= R2*V03;
	V14= R2*V13;
	V24=R2*V23;	

	arma::mat M = arma::ones(4,4);
	M= R2*R12*R11*T;

	//paso para obtener la inversa
	arma::mat invM = arma::ones(4,4);

	arma::mat Tinv= arma::eye(4,4);
	Tinv(1,0)=V0(0);
	Tinv(2,0)= V0(1);
	Tinv(3,0)= V0(2);
	invM=Tinv* arma::strans(R11) * arma::strans(R12) * arma::strans(R2);

	arma::mat union1 = arma::join_cols(M, invM);
	arma::mat union2 =arma:: join_cols(arma::strans(V14), arma::strans(V24));
	arma::mat union3 =arma:: join_cols(union1, arma::strans(V04));
	return arma:: join_cols(union3, union2);
}	