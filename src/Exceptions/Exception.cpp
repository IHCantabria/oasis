
#include <iostream>
#include <string>
#include <sstream>
#include "Exception.hpp"

std::stringstream ss;

// Exception Constructor
Exception::Exception(std::string tagString, std::string probString)
{
    mTag = tagString;
    mProblem = probString;
}

void Exception::PrintDebug() const
{
    std::cerr << "** Error (" << mTag << ") **\n";
    std::cerr << "Problem: " << mProblem << "\n\n";
}

// Attribute Error constructor
ValueError::ValueError(std::string probString)
    : Exception("ValueError", probString)
{
    mProblem = probString;
}

// Attribute Error constructor
IOError::IOError(std::string probString)
    : Exception("IOError", probString)
{
    mProblem = probString;
}

// Attribute Error constructor
NotImplementedError::NotImplementedError(std::string probString)
    : Exception("NotImplementedError", probString)
{
    mProblem = probString;
}
