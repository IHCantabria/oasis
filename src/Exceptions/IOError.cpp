
#include <iostream>
#include <string>
#include <sstream>
#include "Exception.hpp"


// Attribute Error constructor
IOError::IOError(std::string probString) : Exception("IOError", probString)
{
	mProblem = probString;
}