
#include <cstdio>
#include <string>
#include <sstream>
#include "string.h"
#include "Exceptions/Exception.hpp"

int parse_file(std::string filePath)
{
    // Declare local variables
    char buffer[1000];
    int count_line = 0;
    int num_items = 0;
    int num_head_lines = 0;
    FILE *pFile;

    // Open file
    pFile = fopen(filePath.c_str(), "r");
    if (pFile == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosBodies.dat\n    ->Full path: " << filePath.c_str() << std::endl;
        throw IOError(ss.str());
    }

    // Parse input file
    while (fgets(buffer, sizeof(buffer), pFile) != NULL)
    {
        // Advnace line count
        count_line++;

        // std::cout << buffer << strncmp(buffer, "/", 1) << std::endl;
        if (strncmp(buffer, "///", 3) == 0)
        {
            num_head_lines++;
        }

        if (num_head_lines == 3)
        {
            count_line++;

            if (fgets(buffer, sizeof(buffer), pFile) == NULL)
            {
                std::stringstream ss;
                ss << "Error while parsing file: " << filePath.c_str() << "\n --> Expected item definition at file line:" << count_line << "\n";
                throw ValueError(ss.str());
            }

            num_items++;
            num_head_lines = 0;
        }
    }
    fclose(pFile);

    return num_items;
}