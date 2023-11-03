#include <iostream>
#include <iomanip>
#include <functional>
#include "schrodinger.h"
#include "bisection.h"
#include <cmath>

double verlet_init(double dx, double x_init, double y_init, double y_prime_init, std::function<double(double, double)> func)
{
    return y_init + y_prime_init*dx + 0.5*func(x_init, y_init)*dx*dx;
}

double verlet_proc(double dx, double x, double psi_curr, double psi_prev, std::function<double(double, double)> func)
{
    return 2*psi_curr - psi_prev + func(x, psi_curr) * dx * dx;
}

double verlet(double x_first,
	    double x_last,	
	    double dx, 
	    double y_init, 
	    double y_prime_init, 
	    std::function<double(double, double)> func)
{
    double x_curr = x_first;
    double y_curr = y_init;

    std::cout << x_curr << "," << y_curr << std::endl;
    

    x_curr += dx;
    y_curr = verlet_init(dx, x_curr, y_init, y_prime_init, func);

    std::cout << x_curr << "," << y_curr << std::endl;

    double y_last = y_init;
    
    int num_iters = ((x_last - x_first) / dx);
    double y_next;
    
    for(int i = 0; i < num_iters; i++)
    {
	x_curr += dx;
	y_next = verlet_proc(dx, x_curr, y_curr, y_last, func);

	y_last = y_curr;
	y_curr = y_next;


        std::cout << x_curr << "," << y_curr << std::endl;

    }

    return y_curr;
}

double normalize(double x_first,
	         double x_last,	
	         double dx, 
	         Schrodinger schro)
{
    double norm_factor = 0;

    auto func = [&] (double x, double y) -> double{return schro.query(x, y);};
    double x_curr = x_first;
    double y_curr = schro.initialConditions(x_first);
    double y_last = y_curr;

    norm_factor += fabs(y_curr) * dx;

    x_curr += dx;
    y_curr = schro.initialConditions(x_curr);
    norm_factor += y_curr*y_curr * dx;

    int num_iters = ((x_last - x_first) / dx);
    double y_next;

    
    for(int i = 0; i < num_iters; i++)
    {
	x_curr += dx;
	y_next = verlet_proc(dx, x_curr, y_curr, y_last, func);

	y_last = y_curr;
	y_curr = y_next;
        norm_factor += fabs(y_curr) * dx;
    }

    return norm_factor;
}



double verlet(double x_first,
	      double x_last,	
	      double dx, 
	      Schrodinger schro,
	      bool log)
{
    // Verlet integration to solve second order differential equation with good accuracy
    // Applied specifically for Schrodinger equation
    //
    // Data and equation for schrodinger encapsulated with Schrodinger object "schro"
    
    double norm_factor = normalize(x_first, x_last, dx, schro);
    auto func = [&] (double x, double y) -> double{return schro.query(x, y);};
    double x_curr = x_first;
    double y_curr = schro.initialConditions(x_first);
    double y_last = y_curr;
    double e = schro.getEnergy();

    if(log) std::cout << x_curr << "," << e + y_curr / norm_factor << std::endl;
    
    x_curr += dx;
    y_curr = schro.initialConditions(x_curr);

    if(log) std::cout << x_curr << "," << e + y_curr / norm_factor << std::endl;

    int num_iters = ((x_last - x_first) / dx);
    double y_next;
    
    for(int i = 0; i < num_iters; i++)
    {
	x_curr += dx;
	y_next = verlet_proc(dx, x_curr, y_curr, y_last, func);

	y_last = y_curr;
	y_curr = y_next;


        if(log) std::cout << x_curr << "," << e + y_curr / norm_factor << std::endl;

    }

    return y_curr;

}


double quad_potl(double x)
{
   return 0.5*x*x;
}

double quart_potl(double x)
{
   return 0.25*x*x*x*x;
}

int main()
{
    double dx = 0.05;
    double x_first = -5;
    double x_last = 5;
    double y_init = 0;
    double y_prime_init = 10;

    double e = 0.475;

    auto bisect = Bisection(dx, x_first, x_last, y_init, 0);
    e = bisect.bisection(6.5, 9, quart_potl, false);
    auto schro = Schrodinger(e, quart_potl);

    verlet(x_first, x_last, dx, schro, true);

}
