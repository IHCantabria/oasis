

#ifndef os_tools_hpp___
#define os_tools_hpp___

#include <string>
#include <vector>

bool CheckDirExits(std::string folderPath);
std::string CorrectBackSlashes(std::string path);
void CheckInputFile(std::string inputFile, std::string specsStr);
std::string GetWorkingDir(void);
std::string JoinPath(std::string basePath, std::string subdirName);
std::string JoinPaths(std::string basePath, std::vector<std::string> subDirs);

#endif