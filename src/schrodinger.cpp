#include <functional>
#include "schrodinger.h"
#include <cmath>

Schrodinger::Schrodinger(double E, std::function<double(double)> potential)
{
    // Constructur. Inputs are the eigen-energy: "E", and the potential function: "potential"
 
    this->E = E;
    this->potential = potential;
}
double Schrodinger::query(double x_query, double y_query) 
{
    return 2 * ( potential(x_query) - E ) * y_query;
}

double Schrodinger::initial_conditions(double x_query)
{
    return exp( -sqrt( 2 * ( potential(x_query) - E ) ) * fabs(x_query) );
}


