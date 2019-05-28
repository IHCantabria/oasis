#include <armadillo>
#include <string>
#include "Spring.hpp"

// Lee inputs de los muelles
void Spring::leer_datosSprings(void){
	
	int ii, jj, kk, ll, temp_N; 
	std::string Dummy;
	const int nInored=147; // numero de lineas que se leen para cada nueva linea
	arma::mat temp_vec(3,1);

	SpringVectors.set_size(3,1);
	data_StressStrain.set_size(6,3);

	//Abro el fichero
	std::ifstream datosSprings ("input/datossprings.dat");

	//Ignoro la primera linea del fichero, que contiene el numero de muelles a estudiar
	datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea

	//Ignoro las lineas que ya se han leido
	for(ii=1;ii<nSpring;ii=ii+1){
		for(jj=1;jj<=nInored;jj=jj+1){
			datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	//Ignoro las tres primeras lineas, donde pone "New spring"
	for(ii=1;ii<=3;ii=ii+1){
		datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	//Leo todo
	datosSprings >> stressModelFlag; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> dampingFlag; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_1; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_2; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> L; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for(ii=0;ii<3;ii=ii+1){
		temp_vec = 0.0;
		datosSprings >> temp_vec(0,0); datosSprings >> temp_vec(1,0); datosSprings >> temp_vec(2,0); datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		SpringVectors(ii,0) = temp_vec;
	}

	for(ii=0;ii<2;ii=ii+1){
		datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		for(jj=0;jj<6;jj=jj+1){
			for(kk=0;kk<6;kk=kk+1){
				datosSprings >> SpringMatrix_K(jj,kk,ii);
			}
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	for(ii=0;ii<2;ii=ii+1){
		datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		for(jj=0;jj<6;jj=jj+1){
			for(kk=0;kk<6;kk=kk+1){
				datosSprings >> SpringMatrix_D(jj,kk,ii);
			}
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	for(ii=0;ii<2;ii=ii+1){
		for(jj=0;jj<6;jj=jj+1){
			datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			datosSprings >> temp_N; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			n_StressStrain[jj][ii] = temp_N;
			arma::mat temp_vec2 = arma::zeros(temp_N,1);
			arma::mat temp_mat = arma::zeros(temp_N,6);
			temp_vec2 = 0.0;
			for(ll=0;ll<temp_N;ll=ll+1){
				datosSprings >> temp_vec2(ll,0);
			}
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			data_StressStrain(jj,0) = temp_vec2;
			for(kk=0;kk<6;kk=kk+1){
				temp_mat = 0.0;
				for(ll=0;ll<temp_N;ll=ll+1){
					datosSprings >> temp_mat(ll,kk);
				}
				datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			}
			data_StressStrain(jj,ii+1) = temp_mat;
		}
	}

	//Cierro el fichero
	datosSprings.close();

}

// Calcula las fuerzas que aplica el muelle en los BCPs y las guarda en estos
void Spring::computeSpringForces(void){

	// Calculo el vector en global que une los BCPs
	arma::mat SpringVectorG = SpringBCP[1]->posG_BCP - SpringBCP[0]->posG_BCP;
	arma::mat SpringVectorG_dot = SpringBCP[1]->velG_BCP - SpringBCP[0]->velG_BCP;

	// Calculo los vectores unitarios del muelle en global para cada cuerpo
	arma::field<arma::mat> SpringVectorsG;
	SpringVectorsG.set_size(3,2);
	for(int ii=0;ii<2;ii=ii+1){
		for(int jj=0;jj<3;jj=jj+1){
			SpringVectorsG(jj,ii) = SpringBCP[ii]->RotMat*SpringVectors(jj,0);
		}
	}

	arma::mat SpringStrains = arma::zeros(6,2); // Deformacion en el muelle
	arma::mat SpringStrains_dot = arma::zeros(6,2); // Derivada temporal de la deformacion en el muelle
	arma::mat temp_pos = SpringVectorG_dot.rows(0,2);
	for(int ii=0;ii<2;ii=ii+1){
		for(int jj=0;jj<3;jj=jj+1){
			SpringStrains(jj,ii) = arma::dot(SpringVectorG.rows(0,2),SpringVectorsG(jj,ii)) - (2*(ii == 0)-1)*L*(jj == 0);
			SpringStrains(jj+3,ii) = SpringVectorG(jj+3,0);
			if(dampingFlag == 1){
				SpringStrains_dot(jj,ii) = arma::dot(temp_pos,SpringVectorsG(jj,ii)) - (2*(ii == 0)-1)*L*(jj == 0);
				SpringStrains_dot(jj+3,ii) = SpringVectorG_dot(jj+3,0);
			}
		}
	}

	arma::mat tempF_L = arma::zeros(6,1);
	arma::mat tempF_G = arma::zeros(6,1);
	double tempS;
	int tempN;
	int tempI;
	arma::mat temp_strainData;
	arma::mat temp_stressData;
	double temp_dS;
	arma::mat temp_slope;

	for(int ii=0;ii<2;ii=ii+1){

		tempF_L = 0.0;

		if(stressModelFlag==1){

			temp_strainData = (2*(ii == 0)-1) * SpringStrains.col(ii);
			tempF_L = - SpringMatrix_K.slice(ii)*temp_strainData;

		}else if(stressModelFlag==2){

			for(int jj=0;jj<6;jj=jj+1){

				tempS = SpringStrains(jj,ii);
				tempN = n_StressStrain[jj][ii] - 1;
				temp_strainData = data_StressStrain(jj,0);
				temp_stressData = data_StressStrain(jj,ii+1);

				if((tempS<=arma::max(arma::max(temp_strainData)))&&(tempS>=arma::min(arma::min(temp_strainData)))){

					tempI = arma::as_scalar(arma::find(temp_strainData>=tempS,1,"first")) - 1;
					temp_dS = tempS - temp_strainData(tempI-1,0);
					temp_slope = (temp_stressData.col(tempI)-temp_stressData.col(tempI-1))/(temp_strainData(tempI,0)-temp_strainData(tempI-1,0));
					tempF_L = tempF_L + temp_stressData.col(tempI-1) + temp_dS*temp_slope;


				} else if (tempS>arma::max(arma::max(temp_strainData))){ 

					temp_dS = tempS-temp_strainData(tempN,0);
					temp_slope = (temp_stressData.col(tempN)-temp_stressData.col(tempN-1))/(temp_strainData(tempN,0)-temp_strainData(tempN-1,0));
					tempF_L = tempF_L + temp_stressData.col(tempN) + temp_dS*temp_slope;

				}else if (tempS<arma::min(arma::min(temp_strainData))){

					temp_dS = tempS-temp_strainData(0,0);
					temp_slope = (temp_stressData.col(1)-temp_stressData.col(0))/(temp_strainData(1,0)-temp_strainData(0,0));
					tempF_L = tempF_L + temp_stressData.col(0) + temp_dS*temp_slope;
				}
			}
		}

		if(dampingFlag == 1){
			temp_strainData = (2*(ii == 0)-1) * SpringStrains_dot.col(ii);
			tempF_L = tempF_L - arma::abs(tempF_L) % (SpringMatrix_D.slice(ii)*temp_strainData);
		}

		tempF_G = tempF_L(0,0)*SpringVectorsG(0,ii) + tempF_L(1,0)*SpringVectorsG(1,ii) + tempF_L(2,0)*SpringVectorsG(2,ii);
		tempF_G.rows(3,5) = tempF_L.rows(3,5);
		SpringBCP[ii]->ForceBCP = SpringBCP[ii]->ForceBCP + tempF_G;
	}

}