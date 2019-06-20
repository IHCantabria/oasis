
#include <iostream>
#include <string>
#include <sstream>
#include "Exception.hpp"

// Attribute Error constructor
NotImplementedError::NotImplementedError(std::string probString) : Exception("NotImplementedError", probString)
{
	mProblem = probString;
}