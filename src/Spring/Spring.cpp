
#include <armadillo>
#include <string>
#include "Spring.hpp"

// Lee inputs de los muelles
void Spring::ReadPropertiesASCII(std::string file_path){
	
	int ii, jj, kk, ll, temp_N; 
	std::string Dummy;
	const int nInored=89; // numero de lineas que se leen para cada nueva linea
	arma::mat temp_vec = arma::zeros(3,1);

	SpringVectors.set_size(3,2);
	data_StressStrain.set_size(6,2);

	//Abro el fichero
	std::ifstream datosSprings (file_path);

	//Ignoro la primera linea del fichero, que contiene el numero de muelles a estudiar
	datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');  // El ignore sirve para ignorar el texto de la linea

	//Ignoro las lineas que ya se han leido
	for(ii=0;ii<nSpring;ii=ii+1){
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
	datosSprings >> frictionFlag; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_1; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_2; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	BCP_1 -= 1;
	BCP_2 -= 1;

	//std::cout << "Spring " << nSpring+1 << " : BCP_1 = " << BCP_1 << "; and BCP_2 = " << BCP_2 << std::endl;

	for(jj=0;jj<2;jj=jj+1){
		for(ii=0;ii<3;ii=ii+1){
			temp_vec = arma::zeros(3,1);
			datosSprings >> temp_vec(0,0); datosSprings >> temp_vec(1,0); datosSprings >> temp_vec(2,0); datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			SpringVectors(ii,jj) = temp_vec;
		}
	}

	datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for(jj=0;jj<6;jj=jj+1){
		for(kk=0;kk<6;kk=kk+1){
			datosSprings >> SpringMatrix_K(jj,kk);
		}
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}
	datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for(jj=0;jj<6;jj=jj+1){
		for(kk=0;kk<6;kk=kk+1){
			datosSprings >> SpringMatrix_M(jj,kk);
		}
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for(jj=0;jj<6;jj=jj+1){
		datosSprings >> SpringMatrix_D(jj,0); datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}


	for(jj=0;jj<6;jj=jj+1){			
		datosSprings >> Dummy; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		datosSprings >> temp_N; datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		n_StressStrain[jj] = temp_N;
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
		data_StressStrain(jj,1) = temp_mat;
	}

	//Cierro el fichero
	datosSprings.close();

}

// Calcula las fuerzas que aplica el muelle en los BCPs y las guarda en estos
void Spring::computeSpringForces(void)
{

	double arg1, arg2;

	// Calculo los vectores unitarios del muelle en global para cada cuerpo
	arma::field<arma::mat> SpringVectorsG;
	arma::field<arma::mat> SpringVectorsG_dot;
	SpringVectorsG.set_size(3,2);
	SpringVectorsG_dot.set_size(3,2);
	for(int ii=0;ii<2;ii=ii+1){ // Bucle sobre los dos BCPs
		for(int jj=0;jj<3;jj=jj+1){ // Bucle sobre los tres vectores del muelle
			SpringVectorsG(jj,ii) = SpringBCP[ii]->rotMat*SpringVectors(jj,ii);
			SpringVectorsG_dot(jj,ii) = arma::cross(SpringBCP[ii]->velG_BCP.rows(3,5),SpringVectors(jj,ii));
		}
	}

	// Calculo el vector en global que une los BCPs
	arma::mat L12 = SpringBCP[1]->posG_BCP.rows(0,2) - SpringBCP[0]->posG_BCP.rows(0,2);
	arma::mat L12_dot = SpringBCP[1]->velG_BCP.rows(0,2) - SpringBCP[0]->velG_BCP.rows(0,2);
	arma::mat L21 = -L12;
	arma::mat L21_dot = -L12_dot;

	arma::mat x1 = SpringVectorsG(0,0); arma::mat y1 = SpringVectorsG(1,0); arma::mat z1 = SpringVectorsG(2,0);
	arma::mat x2 = SpringVectorsG(0,1); arma::mat y2 = SpringVectorsG(1,1); arma::mat z2 = SpringVectorsG(2,1);
	arma::mat x1_dot = SpringVectorsG_dot(0,0); arma::mat y1_dot = SpringVectorsG_dot(1,0); arma::mat z1_dot = SpringVectorsG_dot(2,0);
	arma::mat x2_dot = SpringVectorsG_dot(0,1); arma::mat y2_dot = SpringVectorsG_dot(1,1); arma::mat z2_dot = SpringVectorsG_dot(2,1);


	arma::mat SpringStrains = arma::zeros(6,2); // Deformacion en el muelle
	arma::mat SpringStrains_dot = arma::zeros(6,2); // Derivada temporal de la deformacion en el muelle

	SpringStrains(0,0) = arma::dot(L12,x1);
	SpringStrains(1,0) = arma::dot(L12,y1);
	SpringStrains(2,0) = arma::dot(L12,z1);
	SpringStrains(3,0) = atan2(arma::dot(y2,z1), arma::dot(z2,z1));
	SpringStrains(4,0) = asin(arma::dot(x2,z1));
	SpringStrains(5,0) = atan2(-arma::dot(x2,y1), -arma::dot(x2,x1));

	SpringStrains(0,1) = arma::dot(L21,x2);
	SpringStrains(1,1) = arma::dot(L21,y2);
	SpringStrains(2,1) = arma::dot(L21,z2);	
	SpringStrains(3,1) = -atan2(arma::dot(y1,z2), arma::dot(z1,z2));
	SpringStrains(4,1) = -asin(arma::dot(x1,z2));
	SpringStrains(5,1) = -atan2(-arma::dot(x1,y2), -arma::dot(x1,x2));

	// Lo mismo para la derivada temporal
	if(dampingFlag == 1 || frictionFlag == 1)
	{
		SpringStrains_dot(0,0) = arma::dot(L12_dot,x1) + arma::dot(L12,x1_dot);
		SpringStrains_dot(1,0) = arma::dot(L12_dot,y1) + arma::dot(L12,y1_dot);
		SpringStrains_dot(2,0) = arma::dot(L12_dot,z1) + arma::dot(L12,z1_dot);
		SpringStrains_dot(3,0) = ((arma::dot(y2_dot,z1)+arma::dot(y2,z1_dot))*arma::dot(z2,z1) + arma::dot(y2,z1)*(arma::dot(z2_dot,z1)+arma::dot(z2,z1_dot)))/(pow(arma::dot(z2,z1),2)+pow(arma::dot(y2,z1),2));
		SpringStrains_dot(4,0) = (arma::dot(x2_dot,z1)+arma::dot(x2,z1_dot))/sqrt(1.0-pow(arma::dot(x2,z1),2));
		SpringStrains_dot(5,0) = ((arma::dot(x2_dot,y1)+arma::dot(x2,y1_dot))*arma::dot(x2,x1) + arma::dot(x2,y1)*(arma::dot(x2_dot,x1)+arma::dot(x2,x1_dot)))/(pow(arma::dot(x2,y1),2)+pow(arma::dot(x2,x1),2));

		SpringStrains_dot(0,1) = arma::dot(L21_dot,x2) + arma::dot(L21,x2_dot);
		SpringStrains_dot(1,1) = arma::dot(L21_dot,y2) + arma::dot(L21,y2_dot);
		SpringStrains_dot(2,1) = arma::dot(L21_dot,z2) + arma::dot(L21,z2_dot);	
		SpringStrains_dot(3,1) = -((arma::dot(y1_dot,z2)+arma::dot(y1,z2_dot))*arma::dot(z1,z2) + arma::dot(y1,z2)*(arma::dot(z1_dot,z2)+arma::dot(z1,z2_dot)))/(pow(arma::dot(z1,z2),2)+pow(arma::dot(y1,z2),2));
		SpringStrains_dot(4,1) = -(arma::dot(x1_dot,z2)+arma::dot(x1,z2_dot))/sqrt(1.0-pow(arma::dot(x1,z2),2));
		SpringStrains_dot(5,1) = -((arma::dot(x1_dot,y2)+arma::dot(x1,y2_dot))*arma::dot(x1,x2) + arma::dot(x1,y2)*(arma::dot(x1_dot,x2)+arma::dot(x1,x2_dot)))/(pow(arma::dot(x1,y2),2)+pow(arma::dot(x1,x2),2));
	}

	arma::mat tempF_L = arma::zeros(6,1); // Fuerza registrada en el muelle en coordenadas del cuerpo
	arma::mat tempF_G = arma::zeros(6,1); // Fuerza registrada en el muelle en coordenadas globales

	// Variables temporales necesarias mas adelante
	double tempS, temp_dS;
	int tempN, tempI;
	arma::mat temp_strainData, temp_stressData, temp_slope;
	arma::mat tempF_K = arma::zeros(6,1);

	for(int ii=0;ii<2;ii=ii+1)
	{ // Bucle sobre los dos BCPs

		tempF_L = arma::zeros(6,1); // Fuerza registrada en el muelle en coordenadas del cuerpo, la reseteo a cero

		if (stressModelFlag==1) { // Modelo de muelle lineal y simetrico

			// Strain calculado previamente, cambiandole el signo dependiendo del BCP en el que se este
			// debido a como se calcula SpringVectorG
			temp_strainData =  SpringStrains.col(ii);
			// Fuerza registrada en el muelle en coordenadas del cuerpo, multiplicando matriz y strain
			tempF_K = SpringMatrix_K*temp_strainData;
			tempF_L = tempF_K;			

		} else if (stressModelFlag==2) { // Modelo de muelle no lineal, input es curva de deformacion-fuerza.

			
			for(int jj=0;jj<6;jj=jj+1){

				tempS = SpringStrains(jj,ii);
				tempN = n_StressStrain[jj] - 1;
				temp_strainData = data_StressStrain(jj,0);
				temp_stressData = data_StressStrain(jj,1);

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

		// Fuerza de damping
		if(dampingFlag == 1){
			temp_strainData = SpringStrains_dot.col(ii);
			tempF_L = tempF_L - (SpringMatrix_D % arma::abs(tempF_K)) % temp_strainData;
		}

		// Fuerza de friccion
		if(frictionFlag == 1){
			temp_strainData = SpringStrains_dot.col(ii);
			tempF_L = tempF_L - (SpringMatrix_M * arma::abs(tempF_K)) % temp_strainData;
		}

		// Paso de coordenadas locales a coordenadas globales la fuerza del muelle
		tempF_G = arma::zeros(6,1);
		tempF_G.rows(0,2) = tempF_L(0,0)*SpringVectorsG(0,ii) + tempF_L(1,0)*SpringVectorsG(1,ii) + tempF_L(2,0)*SpringVectorsG(2,ii);
		tempF_G.rows(3,5) = tempF_L(3,0)*SpringVectorsG(0,ii) + tempF_L(4,0)*SpringVectorsG(1,ii) + tempF_L(5,0)*SpringVectorsG(2,ii);

		SpringBCP[ii]->forceBcp = SpringBCP[ii]->forceBcp + tempF_G; // Acumulo la fuerza obtenida en el BCP

		//std::cout << "Spring " << nSpring+1 << ", force at BCP " << ii+1 << " in global frame is:   "  << std::endl << SpringBCP[ii]->forceBcp << std::endl;
		//std::cout << "Spring " << nSpring+1 << ", position of BCP " << ii+1 << " in global frame is:   "  << std::endl << SpringBCP[ii]->posG_BCP << std::endl;
		//std::cout << "Spring " << nSpring+1 << ", strain at BCP " << ii+1 << " is:   "  << std::endl << SpringStrains.col(ii) << std::endl;
		//std::cout << "Spring " << nSpring+1 << " BCP " << ii+1 << " is BCP "  <<  SpringBCP[ii]->GetId() << std::endl;

	}

	//std::cout << "Spring " << nSpring+1 << ", L12 = " << std::endl << L12 << std::endl;
	//std::cout << "Spring " << nSpring+1 << ", L21 = " << std::endl << L21 << std::endl;

	
	//std::cout << "Spring " << nSpring+1 << ", position increment of BCPs 1 -> 2 in global frame is:   "  << std::endl << SpringBCP[1]->posG_BCP-SpringBCP[0]->posG_BCP << std::endl;
	//std::cout << "Spring " << nSpring+1 << ", strain at BCP 1 is:   "  << std::endl << SpringStrains.col(0) << std::endl;
	//std::cout << "Spring " << nSpring+1 << ", force at BCP 1 in global frame is:   "  << std::endl << SpringBCP[0]->forceBcp << std::endl << std::endl;


	//std::cout << "Spring " << nSpring+1 << ", position increment of BCPs 2 -> 1 in global frame is:   "  << std::endl << SpringBCP[0]->posG_BCP-SpringBCP[1]->posG_BCP << std::endl;
	//std::cout << "Spring " << nSpring+1 << ", strain at BCP 2 is:   "  << std::endl << SpringStrains.col(1) << std::endl;
	//std::cout << "Spring " << nSpring+1 << ", force at BCP 2 in global frame is:   "  << std::endl << SpringBCP[1]->forceBcp << std::endl << std::endl;

}