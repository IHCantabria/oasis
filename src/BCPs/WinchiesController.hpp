
#ifndef WINCHIE_CONTROLER_FLAG
#define WINCHIE_CONTROLER_FLAG
#include <armadillo>
#include <string>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <ctime>
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "Winchies.hpp"

class Simulation;

class WinchieController
{
public:
	int nWinchies;		// Number of winchies
	Winchie **Winchies; // Pointers a los winchies
	Simulation *pSim;	// Pointer to simulation instance
	arma::mat T;		// Tension vector [N]

	// Output files (common to all controllers)
	FILE *pfile_TW; // WinchesTensions.txt
	FILE *pfile_LL; // WinchedLinesLengths.txt

	WinchieController(void);
	WinchieController(int n, Winchie **Ws, Simulation *pIncSim);
	virtual ~WinchieController();

	virtual void ReadPropertiesASCII(FILE *pFile) = 0;
	virtual void SetUpWinchiesController(void);
	virtual void controlWinchies(double time) = 0;
	virtual void OpenOutputFilesASCII(std::string path);
	virtual void CloseOutputFilesASCII(void);
	virtual void WriteOut(double t);

	// Factory method to create the appropriate controller type
	static WinchieController *Create(int controllerType, int n, Winchie **Ws, Simulation *pIncSim);

protected:
	void applyTensions(void);
};

#endif
