// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#endif

#ifndef os_tools_hpp___
#define os_tools_hpp___
#include "os_tools.hpp"
#include <cstdio>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <armadillo>
#include "./Exceptions/Exception.hpp"

// Filesystem: prefer C++17, fall back to experimental, then boost
#if __cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)
#include <filesystem>
namespace fs = std::filesystem;
#elif defined(__has_include) && __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#elif defined(__has_include) && __has_include(<boost/filesystem.hpp>)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#include <windows.h>

std::string GetWorkingDir()
{
    char buffer[MAX_PATH];
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::string::size_type pos = std::string(buffer).find_last_of("\\/");
    return std::string(buffer).substr(0, pos);
}
#elif defined(__unix__)
#include <unistd.h>

std::string GetWorkingDir()
{
    char buffer[FILENAME_MAX];
    getcwd(buffer, FILENAME_MAX);
    std::string::size_type pos = std::string(buffer).find_last_of("\\/");
    return std::string(buffer).substr(0, pos);
}
#endif

bool CheckDirExits(std::string folderPath)
{

    bool exists = false;
    struct stat info;

    if (stat(folderPath.c_str(), &info) != 0) // S_ISDIR() doesn't exist on my windows
        exists = false;
    else if (info.st_mode & S_IFDIR)
        exists = true;
    else
        exists = false;

    return exists;
}

inline bool CheckFileExists(const std::string& name)
{
    if (FILE* file = fopen(name.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    else
    {
        return false;
    }
}

void CheckInputFile(std::string inputFile, std::string specsStr)
{
    if (inputFile.length() > 4)
    {
        if (inputFile.substr(inputFile.length() - 4, 4).compare(".dat") != 0)
        {
            std::stringstream ss;
            ss << "Invalid file data format for: " << specsStr.c_str();
            throw ValueError(ss.str());
        }
        else if (!CheckFileExists(inputFile))
        {
            std::stringstream ss;
            ss << "Actuator file does not exits for: " << specsStr.c_str();
            throw ValueError(ss.str());
        }
    }
    else
    {
        std::stringstream ss;
        ss << "No input file name for: " << specsStr.c_str();
        throw ValueError(ss.str());
    }
}

std::string CorrectBackSlashes(std::string path)
{
    int pos_backslash = path.find("\\", 0);
    while (pos_backslash > -1)
    {
        path.replace(pos_backslash, 1, "/");
        pos_backslash = path.find("\\", pos_backslash);
    }

    return path;
}

std::string JoinPath(std::string basePath, std::string subdirName)
{
    // Declare variables
    fs::path p1 = basePath.c_str();

    // Join paths
    p1 /= subdirName.c_str();
    std::string full_dir = CorrectBackSlashes(p1.u8string());

    return full_dir;
}

std::string JoinPaths(std::string basePath, std::vector<std::string> subDirs)
{
    for (int i = 0; i < subDirs.size(); i++)
    {
        basePath = JoinPath(basePath, subDirs[i]);
    }

    return basePath;
}

void WriteASCII(std::string filePath, arma::mat outVec, bool showColsNum)
{
    std::ofstream xpos;
    xpos.open(filePath);
    if (showColsNum)
    {
        for (int ii = 0; ii < outVec.n_rows; ii++)
            xpos << std::setw(15) << ii;
        xpos << std::endl;
    }
    for (int ii = 0; ii < outVec.n_rows; ii++)
        xpos << std::setw(15) << outVec(ii, 0);
    xpos << std::endl;
    xpos.close();
}

#endif
