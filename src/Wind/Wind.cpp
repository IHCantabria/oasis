
#include <armadillo>
#include <tuple>
#include "Wind.hpp"
#include "../Exceptions/Exception.hpp"
#include "../MathTools.hpp"
#include "../os_tools.hpp"

Wind::Wind(double S, double D)
{
    speed = S;
    heading = D*pi/180.0;
}