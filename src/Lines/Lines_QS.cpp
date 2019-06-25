
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "Lines.hpp"

extern double PI;
extern double g;
extern double rhoW;
extern double fondo;

void Line::qs_Functions(double& ff,double& gg,double& DfDH,double& DfDV,double& DgDH,double& DgDV){

	/*
    This function evaluates functions f and g definned by:
    f(H_{F},V_{F}) = x_{F}(H_{F},V_{F}) - x_{end} 
    g(H_{F},V_{F}) = z_{F}(H_{F},V_{F}) - z_{end}
    It also evaluates their derivatives with respect to H_{F} and V_{F}.
    Taken from qs_GetTen which replicates Jonkman's PhD [Jonkman-2007]
    equations (2-35a) - (2.40)
	*/

	double om, temp, temp2, tempSq, tempLg, tempOm, tempOm2, tempOmSq, tempOmLg, LmVFpOm;

	om = (rho0-rhoW*A)*g;
	temp = VF/HF;
	tempOm = ( VF - om*L) / HF;
	temp2 = temp*temp;
	tempOm2 = tempOm*tempOm;
	tempSq = sqrt(1.0 + temp2);
	tempOmSq = sqrt(1.0 + tempOm2);
	tempLg = log(temp+tempSq);
	tempOmLg = log(tempOm+tempOmSq);

	if((floor_flag < 0.0) || (om < 0.0) || (VF - om*L > 0.0)){

		ff = HF*L/EA + (HF/om) * (tempLg-tempOmLg)- xF;
		gg = (1.0/EA) * (VF * L - 0.5*om*L*L) + (HF/om) * (tempSq - tempOmSq)-zF;

		DfDH = L/EA + (tempLg - tempOmLg) / om  - ((temp + temp2/tempSq) / (temp + tempSq)  -  (tempOm + tempOm2/tempOmSq) / (tempOm + tempOmSq) ) / om;
		DfDV = ((1.0 + temp/tempSq) / (temp + tempSq)  - (1.0 + tempOm/tempOmSq) / (tempOm + tempOmSq)) / om;
		DgDH = (tempSq - tempOmSq) / om - ((temp2/tempSq) - (tempOm2 / tempOmSq)) / om;
		DgDV = L/EA + ( (temp / tempSq) - (tempOm / tempOmSq)) / om;

	}else if(-CB*(VF - om*L) < HF){

		LmVFpOm = L - (VF/om);

		ff = HF*L/EA + ((HF/om) * tempLg) + LmVFpOm +  CB * om *0.5 /(EA) * ( -LmVFpOm*LmVFpOm)-xF;
		gg = (1.0/EA) * (VF * L - 0.5*om*L*L) + (HF/om) * (tempSq - 1.0)-zF;

		DfDH = L/EA + (tempLg) / om  - ((temp + temp2/tempSq) / (temp + tempSq) ) / om;
		DfDV = ((1.0 + temp/tempSq) / (temp + tempSq) ) / om  + (CB/EA) * LmVFpOm - 1.0/om;
		DgDH = (tempSq - 1.0 - (temp2/tempSq) ) / om;
		DgDV = L/EA + (temp / tempSq) / om;

	}else{

		LmVFpOm = L - (VF/om);

		ff = HF*L/(EA) + (HF/om) * tempLg + LmVFpOm  + CB * om *0.5 /(EA) * ( -LmVFpOm*LmVFpOm  + (LmVFpOm - HF/(CB*om)) * (LmVFpOm - HF/(CB*om)))-xF;
		gg = (1.0/EA) * (VF * L - 0.5*om*L*L) + (HF/om) * (tempSq - 1.0)-zF;


		DfDH = L/(EA) + (tempLg) / om  - ((temp + temp2/tempSq) / (temp + tempSq) ) / om  - (LmVFpOm - (HF/(CB*om))) / (EA);
		DfDV = ((1.0 + temp/tempSq) / (temp + tempSq) ) / om  + HF/(om*EA)  - 1.0/om;
		DgDH = (tempSq - 1.0 - (temp2/tempSq) ) / om; 
		DgDV = L/EA + (temp / tempSq) / om;

	};
}

void Line::qs_GetTen(void
){

	/*
    This function provides the fairlead's tension of the catenary.
    The method of Newton-Raphson is employed here. Functions f and g, and their
    derivatives are evaluated by function qs_Functions.
    Taken from GetFairlairTensions.m which replicates Jonkman's PhD [Jonkman-2007]
    equations (2-35a) - (2.40)
	*/

	double lambda, dH, deter, dV;
	double tol = 1.0e-5;
	int nIter=0;
	int nMaxIter=1000;
	double om = (rho0-rhoW*A)*g;
	double ff, gg, DfDH, DfDV, DgDH, DgDV;

	//Initial condition
	lambda=sqrt(3.0*((L*L-zF*zF)/(xF*xF)-1.0));


	HF=std::abs(om*xF*0.5/lambda);

	VF=0.5*om*(zF/tanh(lambda)+L);

	//Using Newton-Raphso
	dH=2*tol;

	while((nIter<nMaxIter)&&(std::abs(dH)>tol)){

		nIter=nIter+1;
		this->qs_Functions(ff,gg,DfDH,DfDV,DgDH,DgDV);

		//Compute the determinant of the Jacobian matrix
		deter = DfDH * DgDV - DfDV * DgDH;
		if (std::abs(deter) < 1e-12){
			std::cout << "WARNING: In QS method, singular Jacobian"<<std::endl<<std::endl;
			break;
		};

		//Apply that the increment in the iterant is \De x_{n} = - inv(Jac) * f(x_{n})
		dH =   -( DgDV * ff - DfDV * gg ) / deter;
		dV =   ( DgDH * ff - DfDH * gg ) / deter;

		//Update the iterant 
		HF = HF + dH;
		VF = VF + dV;
		//To avoid problems, we impose Tol as the lower limit
		HF = std::max(HF, tol);
		VF = std::max(VF, tol);
	};

	//std::cout<<"nIter QS:  "<<nIter<<std::endl<<std::endl;
	//std::cout<<std::endl<<"Tension obtenida:  "<<HF<<"   "<<VF<<std::endl<<std::endl;
}

void Line::qs_Solution(void){

	int ii, jj, nIter, nNodosFondo;
	double om = (rho0-rhoW*A)*g;
	arma::mat temp1 = arma::zeros(N,1);
	double temp2;

	if ( (floor_flag < 0.0) || (om < 0.0) || (VF > om*L) ) {
		HA=HF;
		VA=VF;
		for(ii=0;ii<N;ii=ii+1){
			temp1(ii,0) = ( VF - om*L + om*s(ii,0)) / HF;
			temp2 = ( VF - om*L) / HF;
			xc(ii,0) = HF*s(ii,0)/EA + (HF/om) * (log(temp1(ii,0) + sqrt(1.0+temp1(ii,0)*temp1(ii,0)))-log(temp2 + sqrt(1.0+temp2*temp2)));
			zc(ii,0) = s(ii,0)/(EA) * (VF - om*L + 0.5*om*s(ii,0)) +  (HF/om) * (sqrt(1.0+temp1(ii,0)*temp1(ii,0)) - sqrt(1.0+temp2*temp2));
			dxcds(ii,0) = HF/EA + (HF/om)*( (om/HF + om*temp1(ii,0)/(HF*sqrt(1.0+temp1(ii,0)*temp1(ii,0))) ) / (temp1(ii,0)+sqrt(1.0+temp1(ii,0)*temp1(ii,0))));
			dzcds(ii,0) = (VF - om*L + om*s(ii,0))/EA + (HF/om)*(om*temp1(ii,0)/(HF*sqrt(1.0+temp1(ii,0)*temp1(ii,0))));
			Te(ii,0) = sqrt(HF*HF + (VF - om*L + om*s(ii,0))*(VF - om*L + om*s(ii,0)));
		}
	} else if (-CB * (VF - om*L) < HF) {
		HA = HF + CB * (VF - om*L);
		VA = 0.0;
		for(ii=0;ii<N;ii=ii+1){
			temp1(ii,0) = ( VF - om*L + om*s(ii,0)) / HF;
			temp2 = ( VF - om*L) / HF;
			if( s(ii,0) <= L - VF/om ) {
				xc(ii,0) = s(ii,0) + (s(ii,0) / EA) * (HF + CB * (VF-om*L +0.5*om*s(ii,0)*CB));
				zc(ii,0) = 0.0;
				dxcds(ii,0) = 1.0 +(HF + CB*( VF - om*L + om*s(ii,0)*CB) )/EA;
				dzcds(ii,0) = 0.0;
				Te(ii,0) = HF + CB * (VF-om*L  + om*s(ii,0));
			} else {
				xc(ii,0) = HF*s(ii,0)/EA + (HF/om) * (log(temp1(ii,0) + sqrt(1.0+temp1(ii,0)*temp1(ii,0))))  + L - VF/om - 0.5*CB*( VF - om*L)*( VF - om*L)/(om*EA);
				zc(ii,0) = (HF/om) * (-1.0+sqrt(1.0+temp1(ii,0)*temp1(ii,0))) + s(ii,0)/(EA) * (VF - om*L + 0.5*om*s(ii,0)) + 0.5 * ( VF - om*L)*( VF - om*L)/(om*EA);
				dxcds(ii,0) = HF/EA + (HF/om)*( (om/HF + om*temp1(ii,0)/(HF*sqrt(1.0+temp1(ii,0)*temp1(ii,0))) ) / (temp1(ii,0)+sqrt(1.0+temp1(ii,0)*temp1(ii,0))) );
				dzcds(ii,0) = (HF/om)*(om*temp1(ii,0)/(HF*sqrt(1.0+temp1(ii,0)*temp1(ii,0)))) + (VF - om*L + om*s(ii,0))/EA;
				Te(ii,0) = sqrt(HF*HF + (VF - om*L + om*s(ii,0))*(VF - om*L + om*s(ii,0)));
			}
		}
	} else {
		HA = 0.0;
		VA = 0.0;
		for(ii=0;ii<N;ii=ii+1){
			temp1(ii,0) = ( VF - om*L + om*s(ii,0)) / HF;
			temp2 = ( VF - om*L) / HF;
			if (s(ii,0) <= L - VF/om - HF/(om*CB)) {
				xc(ii,0)=s(ii,0);
				zc(ii,0) = 0.0;
				dxcds(ii,0) = 1.0;
				dzcds(ii,0) = 0.0;
				Te(ii,0) = 0.0;
			} else if (s(ii,0) <= L - VF/om) {
				xc(ii,0) = s(ii,0) - (L - (VF/om) -0.5*HF/(om*CB) ) * (HF/EA) + (s(ii,0) / EA ) * (HF + CB * (VF-om*L) +0.5*om*s(ii,0)*CB) + 0.5 * CB*(VF-om*L)*(VF-om*L)/(om*EA);
				zc(ii,0) = 0.0;				
				dxcds(ii,0) = 1.0 + (HF + CB*( VF - om*L ) + om*s(ii,0)*CB )/EA;
				dzcds(ii,0) = 0.0;
				Te(ii,0) = HF + CB * (VF-om*L  + om*s(ii,0));
			} else {
				xc(ii,0) = HF*s(ii,0)/(EA) + (HF/om) * (log(temp1(ii,0) + sqrt(1.0+temp1(ii,0)*temp1(ii,0)))) + L - VF/om - ( L - VF/om - 0.5 * HF / (om*CB))* HF/(EA);
				zc(ii,0) = (HF/om) * (-1.0+sqrt(1.0+temp1(ii,0)*temp1(ii,0))) + s(ii,0)/(EA) * (VF - om*L + 0.5*om*s(ii,0)) + 0.5 * ( VF - om*L)*( VF - om*L)/(om*EA);
				dxcds(ii,0) = HF/EA + (HF/om)*( (om/HF + om*temp1(ii,0)/(HF*sqrt(1.0+temp1(ii,0)*temp1(ii,0))) ) / (temp1(ii,0)+sqrt(1.0+temp1(ii,0)*temp1(ii,0))) );
				dzcds(ii,0) = (HF/om)*(om*temp1(ii,0)/(HF*sqrt(1.0+temp1(ii,0)*temp1(ii,0)))) + (VF - om*L + om*s(ii,0))/EA;
				Te(ii,0) = sqrt(HF*HF + (VF - om*L + om*s(ii,0))*(VF - om*L + om*s(ii,0)));
			}
		}
	}
}
