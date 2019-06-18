
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
	std::cerr << "** Error ("<<mTag<<") **\n";
	std::cerr << "Problem: " << mProblem << "\n\n";
}
