
#include <armadillo>
#include <string>
#include "WinchiesController.hpp"
#include "WinchiesControllerConstantTension.hpp"
#include "WinchiesControllerHorizontal.hpp"
#include "../Simulations/Simulation.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"

WinchieController::WinchieController(void)
{
    nWinchies = 0;
    pfile_TW = NULL;
    pfile_LL = NULL;
}

WinchieController::WinchieController(int n, Winchie** Ws, Simulation* pIncSim)
{
    nWinchies = n;
    Winchies = Ws;
    pSim = pIncSim;
    pfile_TW = NULL;
    pfile_LL = NULL;
}

WinchieController::~WinchieController()
{
}

void WinchieController::SetUpWinchiesController(void)
{
}

void WinchieController::applyTensions(void)
{
    for (int ii = 0; ii < nWinchies; ii++)
        Winchies[ii]->tau = T(ii, 0) * Winchies[ii]->radius;
}

void WinchieController::OpenOutputFilesASCII(std::string path)
{
    std::string file_path;

    file_path = JoinPath(path, "WinchesTensions.txt");
    pfile_TW = fopen(file_path.c_str(), "w");
    if (pfile_TW == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: WinchesTensions.txt\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }

    file_path = JoinPath(path, "WinchedLinesLengths.txt");
    pfile_LL = fopen(file_path.c_str(), "w");
    if (pfile_LL == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: WinchedLinesLengths.txt\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
}

void WinchieController::CloseOutputFilesASCII(void)
{
    if (pfile_TW)
        fclose(pfile_TW);
    if (pfile_LL)
        fclose(pfile_LL);
}

void WinchieController::OpenOutputFilesCSV(std::string path)
{
    std::string file_path;

    file_path = JoinPath(path, "winches_tensions.csv");
    pfile_TW_csv = fopen(file_path.c_str(), "w");
    if (pfile_TW_csv == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: winches_tensions.csv\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    // Write header
    fprintf(pfile_TW_csv, "time");
    for (int ii = 0; ii < nWinchies; ii++)
        fprintf(pfile_TW_csv, ",w%d_ten", ii);
    fprintf(pfile_TW_csv, "\n");

    file_path = JoinPath(path, "winches_lengths.csv");
    pfile_LL_csv = fopen(file_path.c_str(), "w");
    if (pfile_LL_csv == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: winches_lengths.csv\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    // Write header
    fprintf(pfile_LL_csv, "time");
    for (int ii = 0; ii < nWinchies; ii++)
        fprintf(pfile_LL_csv, ",w%d_len", ii);
    fprintf(pfile_LL_csv, "\n");
}

void WinchieController::CloseOutputFilesCSV(void)
{
    if (pfile_TW_csv)
        fclose(pfile_TW_csv);
    if (pfile_LL_csv)
        fclose(pfile_LL_csv);
}

void WinchieController::OpenOutputFiles(std::string path)
{
    if (pSim->outputFormat == 1)
        OpenOutputFilesCSV(path);
    else
        OpenOutputFilesASCII(path);
}

void WinchieController::CloseOutputFiles(void)
{
    if (pSim->outputFormat == 1)
        CloseOutputFilesCSV();
    else
        CloseOutputFilesASCII();
}

void WinchieController::WriteOut(double t)
{
    if (pSim->outputFormat == 1)
    {
        fprintf(pfile_TW_csv, "%f", t);
        for (int ii = 0; ii < nWinchies; ii++)
            fprintf(pfile_TW_csv, ",%f", T(ii, 0));
        fprintf(pfile_TW_csv, "\n");

        fprintf(pfile_LL_csv, "%f", t);
        for (int ii = 0; ii < nWinchies; ii++)
            fprintf(pfile_LL_csv, ",%f", Winchies[ii]->LineW->L * Winchies[ii]->LineW->dL / Winchies[ii]->LineW->dL0);
        fprintf(pfile_LL_csv, "\n");
    }
    else
    {
        fprintf(pfile_TW, "%f    ", t);
        for (int ii = 0; ii < nWinchies; ii++)
            fprintf(pfile_TW, "%f    ", T(ii, 0));
        fprintf(pfile_TW, "\n");

        fprintf(pfile_LL, "%f    ", t);
        for (int ii = 0; ii < nWinchies; ii++)
            fprintf(pfile_LL, "%f    ", Winchies[ii]->LineW->L * Winchies[ii]->LineW->dL / Winchies[ii]->LineW->dL0);
        fprintf(pfile_LL, "\n");
    }
}

WinchieController* WinchieController::Create(int controllerType, int n, Winchie** Ws, Simulation* pIncSim)
{
    switch (controllerType)
    {
    case 1:
        return new WinchieControllerConstantTension(n, Ws, pIncSim);
    case 2:
        return new WinchieControllerHorizontal(n, Ws, pIncSim);
    default:
    {
        std::stringstream ss;
        ss << "Winch controller type " << controllerType << " is not implemented.\n"
           << "Available types: 1 (Constant tension), 2 (Horizontal control).\n";
        throw ValueError(ss.str());
    }
    }
}
