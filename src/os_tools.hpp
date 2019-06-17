

#ifndef os_tools_hpp___
#define os_tools_hpp___

#include <string>
#include <vector>

std::string GetWorkingDir(void);
std::string CorrectBackSlashes(std::string path);
std::string JoinPath(std::string basePath, std::string subdirName);
std::string JoinPaths(std::string basePath, std::vector<std::string> subDirs);

#endif