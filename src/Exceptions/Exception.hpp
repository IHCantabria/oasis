
#ifndef exceptiondef_hpp__
#define exceptiondef_hpp__

#include <string>

class Exception
{
public:
	std::string mTag, mProblem;

	Exception(std::string tagString, std::string probString);
	void PrintDebug() const;
};

class ValueError :  public Exception
{
public:
	ValueError(std::string probString);
};


class IOError :  public Exception
{
public:
	IOError(std::string probString);
};

#endif //exceptiondef_hpp__