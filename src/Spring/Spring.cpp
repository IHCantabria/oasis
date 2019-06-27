
#include <armadillo>
#include <string>
#include "Spring.hpp"

// Lee inputs de los muelles
void Spring::ReadPropertiesASCII(std::string file_path){
	
	int ii, jj, kk, ll, temp_N; 
	std::string Dummy;
	const int nInored=146; // numero de lineas que se leen para cada nueva linea
	arma::mat temp_vec = arma::zeros(3,1);

	SpringVectors.set_size(3,1);
	data_StressStrain.set_size(6,3);

	//Abro el fichero
	std::ifstream datosSprings (file_path);

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
	BCP_1 -= 1;
	BCP_2 -= 1;

	for(ii=0;ii<3;ii=ii+1){
		temp_vec = arma::zeros(3,1);
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
			for(ll=0;ll<temp_N;ll=ll+1){
				datosSprings >> temp_vec2(ll,0);
			}
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			data_StressStrain(jj,0) = temp_vec2;
			for(kk=0;kk<6;kk=kk+1){
				temp_mat = arma::zeros(temp_N,6);
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
void Spring::computeSpringForces(void)
{

	// Calculo los vectores unitarios del muelle en global para cada cuerpo
	arma::field<arma::mat> SpringVectorsG;
	SpringVectorsG.set_size(3,2);
	for(int ii=0;ii<2;ii=ii+1){ // Bucle sobre los dos BCPs
		for(int jj=0;jj<3;jj=jj+1){ // Bucle sobre los tres vectores del muelle
			SpringVectorsG(jj,ii) = SpringBCP[ii]->rotMat*SpringVectors(jj,0);
		}
	}

	// Calculo el vector en global que une los BCPs
	arma::mat SpringVectorG = SpringBCP[1]->posG_BCP - SpringBCP[0]->posG_BCP;
	arma::mat SpringVectorG_dot = SpringBCP[1]->velG_BCP - SpringBCP[0]->velG_BCP;

	arma::mat SpringStrains = arma::zeros(6,2); // Deformacion en el muelle
	arma::mat SpringStrains_dot = arma::zeros(6,2); // Derivada temporal de la deformacion en el muelle

	arma::mat temp_pos = SpringVectorG.rows(0,2); // Posicion relativa de BCPs en 3dof, en global
	arma::mat temp_pos_dot = SpringVectorG_dot.rows(0,2);//Derivada temporal de la posicion relativa de BCPs en 3dof

	for(int ii=0;ii<2;ii=ii+1){ // Bucle sobre los dos BCPs
		for(int jj=0;jj<3;jj=jj+1){ // Bucle sobre los tres vectores del muelle

			SpringStrains(jj,ii) = arma::dot(temp_pos,SpringVectorsG(jj,ii)); // Proyecto la posicion relativa global sobre los vectores del muelle globales
			SpringStrains(jj+3,ii) = SpringVectorG(jj+3,0);

			// Lo mismo para la derivada temporal
			if(dampingFlag == 1)
			{
				SpringStrains_dot(jj,ii) = arma::dot(temp_pos_dot,SpringVectorsG(jj,ii));
				SpringStrains_dot(jj+3,ii) = SpringVectorG_dot(jj+3,0);
			}
		}
	}

	arma::mat tempF_L = arma::zeros(6,1); // Fuerza registrada en el muelle en coordenadas del cuerpo
	arma::mat tempF_G = arma::zeros(6,1); // Fuerza registrada en el muelle en coordenadas globales

	// Variables temporales necesarias mas adelante
	double tempS, temp_dS;
	int tempN, tempI;
	arma::mat temp_strainData, temp_stressData, temp_slope;

	for(int ii=0;ii<2;ii=ii+1)
	{ // Bucle sobre los dos BCPs

		tempF_L = arma::zeros(6,1); // Fuerza registrada en el muelle en coordenadas del cuerpo, la reseteo a cero

		if (stressModelFlag==1) { // Modelo de muelle lineal y simetrico

			// Strain calculado previamente, cambiandole el signo dependiendo del BCP en el que se este
			// debido a como se calcula SpringVectorG
			temp_strainData =  (2 * (ii == 0) - 1) * SpringStrains.col(ii);
			// Fuerza registrada en el muelle en coordenadas del cuerpo, multiplicando matriz y strain
			tempF_L = SpringMatrix_K.slice(ii)*temp_strainData;

			// Lo mismo para efectos viscosos
			if(dampingFlag == 1){
				temp_strainData = (2*(ii == 0)-1) * SpringStrains_dot.col(ii);
				tempF_L = tempF_L + SpringMatrix_K.slice(ii) * (SpringMatrix_D.slice(ii)*temp_strainData);
			}

		} else if (stressModelFlag==2) { // Modelo de muelle no lineal, input es curva de deformacion-fuerza.

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

			if(dampingFlag == 1){
				temp_strainData = (2*(ii == 0)-1) * SpringStrains_dot.col(ii);
				tempF_L = tempF_L - arma::abs(tempF_L) % (SpringMatrix_D.slice(ii)*temp_strainData);
			}
		}

		// Paso de coordenadas locales a coordenadas globales la fuerza del muelle
		tempF_G.rows(0,2) = tempF_L(0,0)*SpringVectorsG(0,ii) + tempF_L(1,0)*SpringVectorsG(1,ii) + tempF_L(2,0)*SpringVectorsG(2,ii);
		tempF_G.rows(3,5) = SpringBCP[ii]->rotMat*tempF_L.rows(3,5);

		SpringBCP[ii]->forceBcp = SpringBCP[ii]->forceBcp + tempF_G; // Acumulo la fuerza obtenida en el BCP

	}

}