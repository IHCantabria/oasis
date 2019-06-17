

#ifndef os_tools_hpp___
#define os_tools_hpp___
#include "os_tools.hpp"
#include <cstdio>
#include <string>
#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
namespace fs = std::experimental::filesystem;

#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
	#include <windows.h>
	
	std::string GetWorkingDir() {
		char buffer[MAX_PATH];
		GetModuleFileName( NULL, buffer, MAX_PATH );
		std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
		return std::string( buffer ).substr( 0, pos);
	}
#elif defined(__unix__)
	#include <unistd.h>
	
	std::string GetWorkingDir() {
		char buffer[FILENAME_MAX];
		getcwd(buffer, FILENAME_MAX);
		std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
		return std::string( buffer ).substr( 0, pos);
	}
#endif


bool CheckDirExits(std::string folderPath)
{
	
	bool exists = false;
	struct stat info;

	if( stat( folderPath.c_str(), &info ) != 0 ) // S_ISDIR() doesn't exist on my windows 
		exists = false;
	else if( info.st_mode & S_IFDIR )  
		exists = true;
	else
		exists = false;

	return exists;
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
	
	printf("Resulting path is: %s\n", full_dir.c_str());
	
	return full_dir;
}


std::string JoinPaths(std::string basePath, std::vector<std::string> subDirs)
{
	for (int i=0; i<subDirs.size(); i++)
	{
		basePath = JoinPath(basePath, subDirs[i]);
	}
	
	return basePath;
}

#endif