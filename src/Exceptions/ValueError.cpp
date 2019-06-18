
#include <iostream>
#include <string>
#include <sstream>
#include "Exception.hpp"

// Attribute Error constructor
ValueError::ValueError(std::string probString) : Exception("ValueError", probString)
{
	mProblem = probString;
}