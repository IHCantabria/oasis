#include <armadillo>
#include <string>
#include "Morison.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Simulations/Simulation.hpp"


Morison::Morison(int numBodies_inp, Simulation* pSim_inp)
{
	pi = arma::datum::pi;

	numBodies = numBodies_inp;
	pSim = pSim_inp;

	pWindFKCoeff = new arma::cube* [numBodies];
	pWindDragCoeff = new arma::cube* [numBodies];
	pCurrFKCoeff = new arma::cube* [numBodies];
	pCurrDragCoeff = new arma::cube* [numBodies];
}


void Morison::ReadMorisonData(void)
{
    std::string file_path = JoinPath(pSim->inputFolderPath, "dataMorison.dat");

    FILE* file_pointer = fopen(file_path.c_str(), "r");
	
	if (file_pointer == NULL)
	{
        std::stringstream ss;
        ss << "Not possible to open the file: dataMorison.dat\n    ->Dir: " << pSim->inputFolderPath << std::endl;
        throw IOError(ss.str());
	}
	
    char bufferLine [1000];

    // Ignore three header lines for flow definition
	for(int ii=0; ii<3; ii++)
	{
		fgets(bufferLine, sizeof(bufferLine), file_pointer);
	}

	// Get flow type line
    fscanf(file_pointer, "%d %[^\n]\n", &FlowType_flag, bufferLine);

    if (FlowType_flag==1)
    {
    	fgets(bufferLine, sizeof(bufferLine), file_pointer); // Ignore one line
    	double wspd, wdir, cspd, cdir;
    	time = arma::linspace(-1, 1e6, 2);
    	wind_spd = arma::zeros(size(time));
    	wind_dir = arma::zeros(size(time));
    	curr_spd = arma::zeros(size(time));
    	curr_dir = arma::zeros(size(time));

    	fscanf(file_pointer, "%lf %[^\n]\n", &wspd, bufferLine); wind_spd = wind_spd + wspd;
    	fscanf(file_pointer, "%lf %[^\n]\n", &wdir, bufferLine); wind_dir = wind_dir + wdir;
    	fscanf(file_pointer, "%lf %[^\n]\n", &cspd, bufferLine); curr_spd = curr_spd + cspd;
    	fscanf(file_pointer, "%lf %[^\n]\n", &cdir, bufferLine); curr_dir = curr_dir + cdir;

    	// Ignore two lines for variable flow
    	for(int ii=0; ii<2; ii++)
		{
			fgets(bufferLine, sizeof(bufferLine), file_pointer);
		}
    }
    else if (FlowType_flag==2)
    {
    	// Ignore six lines for constant flow
		for(int ii=0; ii<6; ii++)
		{
			fgets(bufferLine, sizeof(bufferLine), file_pointer);
		}

		char cFlowDataFile [1000];
		fscanf(file_pointer, "%s %[^\n]\n", &cFlowDataFile, bufferLine);
		FlowDataFile = cFlowDataFile;
		ReadFlowData_HDF5();
    }
    else
	{
	    std::stringstream ss;
	    ss << "Error while parsing file: dataMorison.dat; Unexpected flow type. \n";
	    throw ValueError(ss.str());
	}


	// Ignore four header lines for coefficients definition
	for(int ii=0; ii<4; ii++)
	{
		fgets(bufferLine, sizeof(bufferLine), file_pointer);
	}

    // Get symetry line
    int Sym_flag;
    fscanf(file_pointer, "%d %[^\n]\n", &Sym_flag, bufferLine);

    if (Sym_flag==1)
    {
    	fgets(bufferLine, sizeof(bufferLine), file_pointer); // Ignore one line

    	char cMorCoeffDataFile [1000];
		fscanf(file_pointer, "%s %[^\n]\n", &cMorCoeffDataFile, bufferLine);
		MorCoeffDataFile = cMorCoeffDataFile;
		ReadMorCoeffData_HDF5();
    }
    else if (Sym_flag==2)
    {
    	// Ignore three lines for non symmetric coefficients definition
		for(int ii=0; ii<3; ii++)
		{
			fgets(bufferLine, sizeof(bufferLine), file_pointer);
		}

		// Get flow type line
    	int SymOrder;
    	fscanf(file_pointer, "%d %[^\n]\n", &SymOrder, bufferLine);
    	// Close file
    	fclose(file_pointer);

    	arma::mat temp = arma::zeros(6,2);
    	arma::mat windFKCoef_X = temp, windDragCoef_X = temp, windFKCoef_Y = temp, windDragCoef_Y = temp;
    	arma::mat currFKCoef_X = temp, currDragCoef_X = temp, currFKCoef_Y = temp, currDragCoef_Y = temp;

    	//Abro el fichero de nuevo
		std::ifstream dataMorCoeff(file_path);
		std::string Dummy; int ii;

		for(ii=1; ii<=20; ii=ii+1){
			dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}

		// Read WIND Froude–Krylov matrix for symmetry axis heading (X)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(int ii=0;ii<6;ii++){
			dataMorCoeff >> windFKCoef_X(ii,0); dataMorCoeff >> windFKCoef_X(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read WIND Drag matrix for symmetry axis heading (X)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> windDragCoef_X(ii,0); dataMorCoeff >> windDragCoef_X(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read WIND Froude–Krylov matrix for heading tangent to symmetry axis (Y) (only for case with '2' symmetry)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> windFKCoef_Y(ii,0); dataMorCoeff >> windFKCoef_Y(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read WIND Drag matrix for heading tangent to symmetry axis (Y) (only for case with '2' symmetry)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> windDragCoef_Y(ii,0); dataMorCoeff >> windDragCoef_Y(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read CURRENTS Froude–Krylov matrix for symmetry axis heading (X)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> currFKCoef_X(ii,0); dataMorCoeff >> currFKCoef_X(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read CURRENTS Drag matrix for symmetry axis heading (X)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> currDragCoef_X(ii,0); dataMorCoeff >> currDragCoef_X(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read CURRENTS Froude–Krylov matrix for heading tangent to symmetry axis (Y) (only for case with '2' symmetry)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> currFKCoef_Y(ii,0); dataMorCoeff >> currFKCoef_Y(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}
		// Read CURRENTS Drag matrix for heading tangent to symmetry axis (Y) (only for case with '2' symmetry)
		dataMorCoeff >> Dummy; dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		for(ii=0;ii<6;ii=ii+1){
			dataMorCoeff >> currDragCoef_Y(ii,0); dataMorCoeff >> currDragCoef_Y(ii,1);  dataMorCoeff.ignore(std::numeric_limits<int>::max(), '\n');
		}

		//Cierro el fichero
		dataMorCoeff.close();

		// Postproceso de los datos leidos para generar las matrices usadas por la clase...
		if (SymOrder>2)
		{
			headings = arma::linspace(0,360,SymOrder+1);
			std::stringstream ss;
		    ss << "FK forces not implemented yet. \n";
		    throw NotImplementedError(ss.str());
		}
		if (SymOrder<2)
		{
			headings = arma::linspace(0,360,SymOrder+1);
			std::stringstream ss;
		    ss << "FK forces not implemented yet. \n";
		    throw NotImplementedError(ss.str());
		}
		if (SymOrder==2)
		{
			headings = arma::linspace(0,360,5);
			arma::cube WindDragCoeff = arma::zeros(5,6,2);
			arma::cube CurrDragCoeff = arma::zeros(5,6,2);
			arma::cube WindFKCoeff = arma::zeros(5,6,2);
			arma::cube CurrFKCoeff = arma::zeros(5,6,2);
			for(ii=0;ii<5;ii=ii+2)
			{
				WindFKCoeff(arma::span(ii),arma::span::all,arma::span::all) = windFKCoef_X;
				WindDragCoeff(arma::span(ii),arma::span::all,arma::span::all) = windDragCoef_X;
				CurrFKCoeff(arma::span(ii),arma::span::all,arma::span::all) = currFKCoef_X;
				CurrDragCoeff(arma::span(ii),arma::span::all,arma::span::all) = currDragCoef_X;
			}
			for(ii=1;ii<5;ii=ii+2)
			{
				WindFKCoeff(arma::span(ii),arma::span::all,arma::span::all) = windFKCoef_Y;
				WindDragCoeff(arma::span(ii),arma::span::all,arma::span::all) = windDragCoef_Y;
				CurrFKCoeff(arma::span(ii),arma::span::all,arma::span::all) = currFKCoef_Y;
				CurrDragCoeff(arma::span(ii),arma::span::all,arma::span::all) = currDragCoef_Y;
			}
			for(ii=0;ii<numBodies;ii=ii+1)
			{
				pWindFKCoeff[ii] = new arma::cube; *pWindFKCoeff[ii] = WindFKCoeff;
				pWindDragCoeff[ii] = new arma::cube; *pWindDragCoeff[ii] = WindDragCoeff;
				pCurrFKCoeff[ii] = new arma::cube; *pCurrFKCoeff[ii] = CurrFKCoeff;
				pCurrDragCoeff[ii] = new arma::cube; *pCurrDragCoeff[ii] = CurrDragCoeff;
			}
		}

    }
    else
	{
	    std::stringstream ss;
	    ss << "Error while parsing file: dataMorison.dat; Unexpected symmetry flag. \n";
	    throw ValueError(ss.str());
	}

	if(arma::accu(arma::abs(wind_spd))>0)
	{
		std::cout << "    --> Morison forces for wind are activated." << std::endl;
		flag_wind = true;
	}

	if(arma::accu(arma::abs(curr_spd))>0)
	{
		std::cout << "    --> Morison forces for currents are activated." << std::endl;
		flag_curr = true;
	}

}


void Morison::ReadFlowData_HDF5(void)
{
	std::cout << "--> Reading Flow (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadFlowData_HDF5 in class Morison not implemented yet. \n";
    throw NotImplementedError(ss.str());
    std::cout << "----> Wave Spectrum Read" << std::endl;
}


void Morison::ReadMorCoeffData_HDF5(void)
{
	std::cout << "--> Reading Morison Coefficients (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadMorCoeffData_HDF5 in class Morison not implemented yet. \n";
    throw NotImplementedError(ss.str());
    std::cout << "----> Wave Spectrum Read" << std::endl;
}


arma::mat Morison::ComputeWindForce(int idBody, double yaw, double t)
{
	arma::mat F = arma::zeros(6,1);

	double wspd, wdir;

	if (FlowType_flag == 1)
	{
		wspd = arma::as_scalar(wind_spd(0,0));
		wdir = arma::as_scalar(wind_dir(0,0));
	}
	if (FlowType_flag == 2)
	{
		//arma::mat tt = arma::zeros(1,1) + t;
		//wspd = arma::as_scalar(arma::interp1(time,wind_spd,tt));
		//wdir = arma::as_scalar(arma::interp1(time,wind_dir,tt));
	    std::stringstream ss;
	    ss << "FK forces not implemented yet. \n";
	    throw NotImplementedError(ss.str());
	}	

	arma::mat vel = arma::zeros(2,1);
	vel(0,0) = wspd*cos(pi*wdir/180);
	vel(1,0) = wspd*sin(pi*wdir/180); 

	arma::mat h = arma::zeros(1,1) + yaw + wdir; // REVISAR SIGNOS

	arma::cube temp_B = interp1(headings,*(pWindDragCoeff[idBody]),h);
	arma::mat B = temp_B(arma::span(0),arma::span::all,arma::span::all);

	F = F + B*(vel%arma::abs(vel));

	return F;
}


arma::mat Morison::ComputeCurrForce(int idBody, double yaw, double t)
{
	arma::mat F = arma::zeros(6,1);

	double cspd, cdir;
	if (FlowType_flag == 1)
	{
		cspd = arma::as_scalar(curr_spd(0,0));
		cdir = arma::as_scalar(curr_dir(0,0));
	}
	if (FlowType_flag == 2)
	{
		//arma::mat tt = arma::zeros(1,1) + t;
		//cspd = arma::as_scalar(arma::interp1(time,curr_spd,tt));
		//cdir = arma::as_scalar(arma::interp1(time,curr_dir,tt));
	    std::stringstream ss;
	    ss << "FK forces not implemented yet. \n";
	    throw NotImplementedError(ss.str());
	}	

	arma::mat vel = arma::zeros(2,1);
	vel(0,0) = cspd*cos(pi*cdir/180);
	vel(1,0) = cspd*sin(pi*cdir/180); 

	arma::mat h = arma::zeros(1,1) + yaw*0 + cdir; // REVISAR SIGNOS

	arma::cube temp_B = interp1(headings,*(pCurrDragCoeff[idBody]),h);
	arma::mat B = temp_B(arma::span(0),arma::span::all,arma::span::all);

	F = F + B*(vel%arma::abs(vel));

	return F;
}