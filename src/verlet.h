#ifndef V_H
#define V_H

#include <functional>


double verlet_init(double dx, double x_init, double y_init, double y_prime_init, std::function<double(double, double)> func);

double verlet_proc(double dx, double x, double psi_curr, double psi_prev, std::function<double(double, double)> func);

double verlet(double x_first,
	    double x_last,	
	    double dx, 
	    double y_init, 
	    double y_prime_init, 
	    std::function<double(double, double)> func);


double verlet(double x_first,
	    double x_last,	
	    double dx, 
	    Schrodinger schro,
	    bool log);
#endif
