
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdio>
#include <cmath>
#include <armadillo>
#include "../Simulations/Simulation.hpp"
#include "Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../os_tools.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"
#include "../Exceptions/Exception.hpp"


Body::Body(int n, Simulation* pIncSim)
{
	id = n;
	pSim = pIncSim;
}


// Calcula el efecto de las fuerzas sobre los BCPs sobre su CDG
void Body::ComputeBcpForces(void)
{

	// Variables locales necesarias usadas a continuacion
	arma::mat posG_temp;
	arma::mat ForceBCP_temp;
	arma::mat F_M;
	arma::mat M_F;
	double M_norm, det;
	double rx, ry, rz, Mx, My, Mz;

	for(int ii=0; ii<numBcps; ii++){ // Bucle sobre los BCPs

		posG_temp = pBodyBcps[ii]->posWrtCdgGlobal; // Guardo en variable temporal la posicion global del BCP
		ForceBCP_temp = pBodyBcps[ii]->forceBcp; // Guardo en variable temporal las fuerzas y momentos sobre el BCP

		F_M = arma::zeros(3,1); // Inicio a cero al fuerza en el cdg causada por el momento sobre el BCP 
		M_norm = arma::norm(ForceBCP_temp.rows(3,5)); // Calculo la norma del momento sobre el BCP

		// Si hay momento en el BCP, calculo la fuerza que produce sobre le cdg
		if (M_norm > 1e-8) {

			// Datos necesarios usados a continuación
			rx = posG_temp(0,0); ry = posG_temp(1,0); rz = posG_temp(2,0);
			Mx = ForceBCP_temp(3,0); My = ForceBCP_temp(4,0); Mz = ForceBCP_temp(5,0);

			det = My*ry*ry + Mx*rx*ry + Mz*ry*rz;

			if (det > 1e-8) {
				// Inversa de la matriz del sistema para resover la ecuacion M = rxF, donde F es desconocido
				//[           Mz*rx, rx*ry, My*ry + Mz*rz]
				//[           Mz*ry,  ry^2,        -Mx*ry]
				//[ - Mx*rx - My*ry, ry*rz,        -Mx*rz] / (My*ry^2 + Mx*rx*ry + Mz*ry*rz)
				arma::mat matrix = arma::zeros(3,3);
				matrix(0,0) =           Mz*rx; matrix(0,1) = rx*ry; matrix(0,2) = My*ry + Mz*rz;
				matrix(1,0) =           Mz*ry; matrix(1,1) = ry*ry; matrix(1,2) =        -Mx*ry;
				matrix(2,0) = - Mx*rx - My*ry; matrix(2,1) = rx*rz; matrix(2,2) =        -Mx*rz;
				matrix = matrix/det;

				// Vector del sistema para resover la ecuacion M = rxF, donde F es desconocido (M - Mt)
				arma::mat vector = ForceBCP_temp.rows(3,5) - arma::dot(posG_temp,ForceBCP_temp.rows(3,5))/arma::dot(posG_temp,posG_temp);
				vector(1,0) = 0.0;

				// Obtengo la fuerza en el cdg causada por el momento en el bcp, ya en global, 
				// por que tanto el brazo posG, como la fuerza ForceBCP, estan en global.
				F_M = matrix*vector;
			}
		}

		M_F = arma::cross(posG_temp,ForceBCP_temp.rows(0,2)); // Momento sobre el cdg causado por la fuerza en el bcp, en global

		bcpForces.rows(0,2) = bcpForces.rows(0,2) + ForceBCP_temp.rows(0,2) + F_M; // Acumulo la fuerza total sobre el cdg en global.
		bcpForces.rows(3,5) = bcpForces.rows(3,5) + rotMat.t() * (ForceBCP_temp.rows(3,5) + M_F); // Acumulo el momento total sobre el cdg en local.
	}

	if (bcpForces.has_nan()){
		std::cout << std::endl << "ERROR: NaN Detected on body with id = " << id << std::endl;
		throw std::exception();
	}
}


int Body::GetId(void)
{
	return id;
}


void Body::LoadHydrodynamicDatabase(Body** hydroDatabaseBodies)
{
	std::cout << "--> Reading Hydrodynamics Properties (HDF5 format)" << std::endl;

	// File path
	std::string file_path = JoinPath(this->pSim->inputFolderPath, this->hydroDatabaseName);

	// Load hydrodynamic database
	this->pHydro = new HydroDatabase(this->hydroDatabaseIndex, this->id, hydroDatabaseBodies, this->pSim);
	std::cout << "    --> Loading hydrodynamic database" << std::endl;
	this->pHydro->LoadHydrodynamicData(file_path);

	std::cout << "    --> Loading and checking the Hydrodynamic C.O.G position" << std::endl;
	// Load and check the Hydrodynamic C.O.G position
	if (this->takeCOGHydroDatabase == 1)
	{
		// Load inital position
		this->pos_eq.rows(0, 2) = this->pHydro->GetCog().t();

		// Add initial position to the global position
		this->pos = this->pos + this->pos_eq;
	}
	else if ((this->takeCOGHydroDatabase == 0) && (this->pHydro->GetNumBodies() > 1))
	{
		// Declare local variables
		bool xcond, ycond, zcond;

		// Check if the input C.O.G is in accordance with the hydrodynamic database
		double cog_tol = 1e-6;
		arma::mat cog = this->pHydro->GetCog();
		xcond = fabs(this->pos_eq(0, 0) - cog(0, 0)) > cog_tol;
		ycond = fabs(this->pos_eq(1, 0) - cog(0, 1)) > cog_tol;
		zcond = fabs(this->pos_eq(2, 0) - cog(0, 2)) > cog_tol;
		
		if (xcond || ycond || zcond)
		{
			std::stringstream ss;
			ss << "Specified Center of Gravity for body: " << this->GetId();
			ss << " mismatch with the C.O.G value in the Hydrodynamic database.\n";
			throw ValueError(ss.str());
		}
	}
	else if (this->takeCOGHydroDatabase == 0)
	{
		this->pos = this->pos + this->pos_eq;
	}
	
	std::cout << "----> Hydrodynamic Properties Read" << std::endl;
}


// Leer datos de los cuerpos
void Body::ReadPropertiesASCII(FILE* pFile)
{
	// Declare variables
	char buffer_line [1000];
	fpos_t carriage_init;
	char cHydroDatabaseName [1000];
	double dtemp;
	int itemp;

	// Read flag to take COG from Hydrodynamic database
	if (fscanf(pFile, "%d %[^\n]\n", &takeCOGHydroDatabase, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "Body: " << this->GetId() <<" - Not possible to read flag to take COG from Hydrodynamic database." << ".\n";
		throw ValueError(ss.str());
	}

	// Read dofs considered
	fgetpos(pFile, &carriage_init);
	while (fscanf(pFile, "%d", &itemp) == 1)
	{
		this->numDofs++;
	}
	fsetpos(pFile, &carriage_init);
	
	this->pDofs = new int [this->numDofs];
	for (int ii=0; ii<this->numDofs; ii++)
	{
		if (fscanf(pFile, "%d", &itemp) != 1) 
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the deegres of freedom for body: " << this->GetId() <<"\n";
			throw ValueError(ss.str());
		}
		this->pDofs[ii] = itemp - 1;
		isDofActive(itemp - 1,0) = 1.0;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read the boundary condition points in the body
	fgetpos(pFile, &carriage_init);
	while (fscanf(pFile, "%d", &itemp) == 1)
	{
		this->numBcps++;
	}
	fsetpos(pFile, &carriage_init);

	this->pIndexBcps = new int[this->numBcps];
	fgetpos(pFile, &carriage_init);
	fgets(buffer_line, sizeof(buffer_line), pFile);
	fsetpos(pFile, &carriage_init);
	std::cout << buffer_line << std::endl;
	for(int ii=0; ii<this->numBcps; ii++)
	{
		if (fscanf(pFile, "%d", &itemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the boundary condition points of the body: " << this->GetId() << " " << "\n";
			throw ValueError(ss.str());
		}
		this->pIndexBcps[ii] = itemp - 1;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read Initial position
	for (int ii=0; ii<6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the initial position of the body: " << this->GetId() << "\n";
			throw ValueError(ss.str());
		}
		
		if (this->takeCOGHydroDatabase == 0)
		{
			this->pos_eq(ii, 0) = dtemp;
			this->pos(ii, 0) = dtemp;
		}
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read Initial displacement from reference position
	for (int ii=0; ii<6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the initial displacement of the body: " << this->GetId() << "\n";
			throw ValueError(ss.str());
		}
		pos(ii, 0) +=  dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read hydrodynamic database filename
	if (fscanf(pFile, "%s %[^\n]\n", cHydroDatabaseName, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the hydrodynamic database name of the body: " << this->GetId() << "\n";
		throw ValueError(ss.str());
	}
	this->hydroDatabaseName = cHydroDatabaseName;

	// Read body index in the associated database
	if (fscanf(pFile, "%d %[^\n]\n", &(this->hydroDatabaseIndex), buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "An error ocurred when trying to read the index of the body in the hydrodynamic database: " << this->GetId() << "\n";
		throw ValueError(ss.str());
	}
	this->hydroDatabaseIndex--;


	
	// Read flag for blocking the body
	if (fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "Body: " << this->GetId() <<" - Not possible to read flag for blocking body." << ".\n";
		throw ValueError(ss.str());
	}
	if (itemp==0){
		flag_blocked = 1;
	} else if (itemp==1){
		flag_blocked = 0;
	} else {
		std::stringstream ss;
		ss << "Body: " << this->GetId() <<" - Flag for blocking body not available, must be 0 or 1." << ".\n";
		throw ValueError(ss.str());
	}

	// Read flag for first order excitation force
	if (fscanf(pFile, "%d %[^\n]\n", &firstOrderExcitationFlag, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "Body: " << this->GetId() <<" - Not possible to read flag for the first order excitation." << ".\n";
		throw ValueError(ss.str());
	}

	// Read flag for second order excitation force
	if (fscanf(pFile, "%d %[^\n]\n", &secondOrderExcitationFlag, buffer_line) != 2)
	{
		std::stringstream ss;
		ss << "Body: " << this->GetId() <<" - Not possible to read flag for the second order excitation." << ".\n";
		throw ValueError(ss.str());
	}


	// Read viscous added mass coefficients
	for (int ii=0; ii<6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the viscous added mass of body: " << this->GetId() << "\n";
			throw ValueError(ss.str());
		}
		A_visc(ii, 0) +=  dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read viscous damping coefficients
	for (int ii=0; ii<6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the viscous damping of body: " << this->GetId() << "\n";
			throw ValueError(ss.str());
		}
		B_visc(ii, 0) +=  dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);

	// Read viscous damping coefficients
	for (int ii=0; ii<6; ii++)
	{
		if (fscanf(pFile, "%lf", &dtemp) != 1)
		{
			std::stringstream ss;
			ss << "An error ocurred when trying to read the viscous damping 2 of body: " << this->GetId() << "\n";
			throw ValueError(ss.str());
		}
		B_visc2(ii, 0) +=  dtemp;
	}
	fscanf(pFile, "%[^\n]\n", buffer_line);





	// Generate array of pointers in order to storage the BCPs pointers
	this->pBodyBcps = new BCP* [this->numBcps];
}


void Body::StoreVelocities(bool restoreMatrix)
{
	if (!restoreMatrix)
	{
		velBuffer.submat(0, pSim->timeBufferCount, 5, pSim->timeBufferCount) = vel;
	}
	else
	{
		int num_points_irf = this->pHydro->GetNumPointsIrf();
		arma::mat velBufferNew = arma::zeros(6, pSim->timeBufferSize);
		velBufferNew.cols(0, num_points_irf-1) = velBuffer.cols(pSim->timeBufferSize-num_points_irf, pSim->timeBufferSize-1);
		velBuffer = velBufferNew;
		velBufferCount = num_points_irf-1;
	}
}


// Actualiza valores del BCP
void Body::UpdateBcps(void)
{
	// Datos necesarios para las matrices de rotacion
	double roll = arma::as_scalar(pos(3,0));
	double pitch = arma::as_scalar(pos(4,0));
	double yaw = arma::as_scalar(pos(5,0));
	double roll_dot = arma::as_scalar(vel(3,0));
	double pitch_dot = arma::as_scalar(vel(4,0));
	double yaw_dot = arma::as_scalar(vel(5,0));
	double roll_dot2 = arma::as_scalar(acc(3,0));
	double pitch_dot2 = arma::as_scalar(acc(4,0));
	double yaw_dot2 = arma::as_scalar(acc(5,0));
	double cr = cos(roll); double sr = sin(roll);
	double cp = cos(pitch); double sp = sin(pitch);
	double cy = cos(yaw); double sy = sin(yaw);

	// Matriz de rotacion calculada como Mz*My*Mx
	arma::mat R1 = arma::eye(3,3); 
	R1(1,1) = cr; R1(1,2) = -sr;
	R1(2,1) = sr; R1(2,2) = cr;
	arma::mat R2 = arma::eye(3,3);
	R2(2,2) = cp; R2(2,0) = -sp;
	R2(0,2) = sp; R2(0,0) = cp; 
	arma::mat R3 = arma::eye(3,3);
	R3(0,0) = cy; R3(0,1) = -sy;
	R3(1,0) = sy; R3(1,1) = cy;

	rotMat = R3*R2*R1;

	arma::mat R1_dot = arma::zeros(3,3); 
	R1_dot(1,1) = -sr; R1_dot(1,2) = -cr;
	R1_dot(2,1) = cr; R1_dot(2,2) = -sr;
	R1_dot = roll_dot*R1_dot;
	arma::mat R2_dot = arma::zeros(3,3);
	R2_dot(2,2) = -sp; R2_dot(2,0) = -cp;
	R2_dot(0,2) = cp; R2_dot(0,0) = -sp;
	R2_dot = pitch_dot*R2_dot;
	arma::mat R3_dot = arma::zeros(3,3);
	R3_dot(0,0) = -sy; R3_dot(0,1) = -cy;
	R3_dot(1,0) = cy; R3_dot(1,1) = -sy;
	R3_dot = yaw_dot*R3_dot;

	rotMat_dot = R3_dot*R2*R1 + R3*R2_dot*R1 + R3*R2*R1_dot;

	arma::mat R1_dot2 = -R1; R1_dot2(0,0) = 0;
	R1_dot2 = roll_dot2*R1_dot + roll_dot*roll_dot*R1_dot2;
	arma::mat R2_dot2 = -R2; R2_dot2(1,1) = 0;
	R2_dot2 = pitch_dot2*R2_dot + pitch_dot*pitch_dot*R2_dot2;
	arma::mat R3_dot2 = -R3; R3_dot2(2,2) = 0;
	R3_dot2 = yaw_dot2*R3_dot + yaw_dot*yaw_dot*R3_dot2;

	rotMat_dot2 = R3_dot2*R2*R1 + R3_dot*R2_dot*R1 + R3_dot*R2*R1_dot + 
	              R3_dot*R2_dot*R1 + R3*R2_dot2*R1 + R3*R2_dot*R1_dot + 
	              R3_dot*R2*R1_dot + R3*R2_dot*R1_dot + R3*R2*R1_dot2;

	// Variable temporal
	arma::mat posG_temp, posL_temp;
	for(int ii=0; ii<numBcps; ii++)
	{
		// Bucle sobre todos los BCPs
		posL_temp = pBodyBcps[ii]->posWrtCdgLocal;
		posG_temp = rotMat*posL_temp; // Brazo cdg-bcp en global
		pBodyBcps[ii]->posWrtCdgGlobal = posG_temp; // Brazo cdg-bcp en global
		pBodyBcps[ii]->rotMat = rotMat; // Matriz de rotacion
		pBodyBcps[ii]->rotMat_dot = rotMat_dot; // Matriz de rotacion
		// Posicion del BCP en global, lo mismo para vel y acc.
		pBodyBcps[ii]->posG_BCP.rows(0,2) = pos.rows(0,2) + posG_temp; // Posicion del BCP en global
		pBodyBcps[ii]->velG_BCP.rows(0,2) = vel.rows(0,2) + rotMat_dot*posL_temp;
		pBodyBcps[ii]->accG_BCP.rows(0,2) = acc.rows(0,2) + rotMat_dot2*posL_temp;
	}
}


// Resetea fuerzas BCPs
void Body::ResetBcps(void)
{
	for(int ii=0; ii<numBcps; ii++)
	{
		// Reseteo a cero la fuerza sobre el BCP
		pBodyBcps[ii]->forceBcp = arma::zeros(6,1);
		pBodyBcps[ii]->temp = arma::zeros(6,1);
	}
	// Reseteo a cero la fuerza total de todos los BCPs
	bcpForces = arma::zeros(6,1);
}


void Body::OpenOutputFilesASCII (std::string path)
{
	char buffer1[50], buffer2[50], buffer3[50], buffer4[50], buffer5[50], buffer6[50], buffer7[50], buffer8[50], buffer9[50], buffer10[50];

	int nn1 = sprintf(buffer1,"DOF_1_Body_%d.txt", GetId());
	int nn2 = sprintf(buffer2,"DOF_2_Body_%d.txt", GetId());
	int nn3 = sprintf(buffer3,"DOF_3_Body_%d.txt", GetId());
	int nn4 = sprintf(buffer4,"DOF_4_Body_%d.txt", GetId());
	int nn5 = sprintf(buffer5,"DOF_5_Body_%d.txt", GetId());
	int nn6 = sprintf(buffer6,"DOF_6_Body_%d.txt", GetId());
	int nn7 = sprintf(buffer7,"HydroStiffnessForce_Body_%d.txt", GetId());
	int nn8 = sprintf(buffer8,"WaveRadiationForce_Body_%d.txt", GetId());
	int nn9 = sprintf(buffer9,"BCPForce_Body_%d.txt", GetId());
	int nn10 = sprintf(buffer10,"WaveExcitationForce_Body_%d.txt", GetId());

	std::string file_path1 = JoinPath(path, buffer1);
	std::string file_path2 = JoinPath(path, buffer2);
	std::string file_path3 = JoinPath(path, buffer3);
	std::string file_path4 = JoinPath(path, buffer4);
	std::string file_path5 = JoinPath(path, buffer5);
	std::string file_path6 = JoinPath(path, buffer6);
	std::string file_path7 = JoinPath(path, buffer7);
	std::string file_path8 = JoinPath(path, buffer8);
	std::string file_path9 = JoinPath(path, buffer9);
	std::string file_path10 = JoinPath(path, buffer10);

	pfile_DOF_1 = fopen (file_path1.c_str(),"w");
	if (pfile_DOF_1 == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn1 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_DOF_2 = fopen (file_path2.c_str(),"w");
	if (pfile_DOF_2 == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn2 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_DOF_3 = fopen (file_path3.c_str(),"w");
	if (pfile_DOF_3 == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn3 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_DOF_4 = fopen (file_path4.c_str(),"w");
	if (pfile_DOF_4 == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn4 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_DOF_5 = fopen (file_path5.c_str(),"w");
	if (pfile_DOF_5 == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn5 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_DOF_6 = fopen (file_path6.c_str(),"w");
	if (pfile_DOF_6 == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn6 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_HSF = fopen (file_path7.c_str(),"w");
	if (pfile_HSF == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn7 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_WRF = fopen (file_path8.c_str(),"w");
	if (pfile_WRF == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn8 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_BCPF = fopen (file_path9.c_str(),"w");
	if (pfile_BCPF == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn9 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
	pfile_WEF = fopen (file_path10.c_str(),"w");
	if (pfile_WEF == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: "<< nn10 <<"\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
	}
}


void Body::CloseOutputFilesASCII (void)
{
	fclose(pfile_DOF_1);
	fclose(pfile_DOF_2);
	fclose(pfile_DOF_3);
	fclose(pfile_DOF_4);
	fclose(pfile_DOF_5);
	fclose(pfile_DOF_6);
	fclose(pfile_HSF);
	fclose(pfile_WRF);
	fclose(pfile_BCPF);
	fclose(pfile_WEF);
}


// Escibir datos a fichero
void Body::WriteOut(double t)
{
	fprintf(pfile_DOF_1, "%f    %f    %f    %f \n",t,this->pos(0,0),this->vel(0,0),this->acc(0,0));
	fprintf(pfile_DOF_2, "%f    %f    %f    %f \n",t,this->pos(1,0),this->vel(1,0),this->acc(1,0));
	fprintf(pfile_DOF_3, "%f    %f    %f    %f \n",t,this->pos(2,0),this->vel(2,0),this->acc(2,0));
	fprintf(pfile_DOF_4, "%f    %f    %f    %f \n",t,this->pos(3,0),this->vel(3,0),this->acc(3,0));
	fprintf(pfile_DOF_5, "%f    %f    %f    %f \n",t,this->pos(4,0),this->vel(4,0),this->acc(4,0));
	fprintf(pfile_DOF_6, "%f    %f    %f    %f \n",t,this->pos(5,0),this->vel(5,0),this->acc(5,0));
	fprintf(pfile_HSF, "%f    ", t);
	for(int ii=0;ii<6;ii=ii+1) fprintf(pfile_HSF, "%f    ", hydrostaticForces(ii, 0));
	fprintf(pfile_HSF, "\n");

	fprintf(pfile_WRF, "%f    ", t);
	for(int ii=0;ii<6;ii=ii+1) fprintf(pfile_WRF, "%f    ", radiationForces(ii, 0));
	fprintf(pfile_WRF, "\n");

	fprintf(pfile_WEF, "%f    ", t);
	for(int ii=0;ii<6;ii=ii+1) fprintf(pfile_WEF, "%f    ", excitationForces_1(ii, 0));
	for(int ii=0;ii<6;ii=ii+1) fprintf(pfile_WEF, "%f    ", excitationForces_2(ii, 0));
	fprintf(pfile_WEF, "\n");

	fprintf(pfile_BCPF, "%f    ", t);
	for(int ii=0;ii<numBcps;ii=ii+1)
	{
		for(int jj=0;jj<6;jj=jj+1)
		{
			fprintf(pfile_BCPF, "%f    ", pBodyBcps[ii]->temp(jj, 0));
			//fprintf(pfile_BCPF, "%f    ", pBodyBcps[ii]->forceBcp(jj, 0));
		}
	}
	for(int jj=0;jj<6;jj=jj+1)
	{
		fprintf(pfile_BCPF, "%f    ", bcpForces(jj, 0));
	}
	fprintf(pfile_BCPF, "\n");
}
