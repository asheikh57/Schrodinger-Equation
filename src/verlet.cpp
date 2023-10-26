#include <iostream>
#include <iomanip>
#include <functional>
#include "schrodinger.h"
#include <cmath>

double schrodinger(double psi, double E, std::function<double(double)> potential)
{
   return 2 * (potential(psi) - E) * psi; 
}

double sho(double x, double y)
{
    return -y;
}

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


double verlet(double x_first,
	    double x_last,	
	    double dx, 
	    Schrodinger schro,
	    std::function<double(double, double)> func,
	    bool log)
{
    double x_curr = x_first;
    double y_curr = schro.initialConditions(x_first);
    double y_last = y_curr;
    if(log) std::cout << x_curr << "," << y_curr << std::endl;
    
    x_curr += dx;
    y_curr = schro.initialConditions(x_curr);

    if(log) std::cout << x_curr << "," << y_curr << std::endl;

    int num_iters = ((x_last - x_first) / dx);
    double y_next;
    
    for(int i = 0; i < num_iters; i++)
    {
	x_curr += dx;
	y_next = verlet_proc(dx, x_curr, y_curr, y_last, func);

	y_last = y_curr;
	y_curr = y_next;


        if(log) std::cout << x_curr << "," << y_curr << std::endl;

    }

    return y_curr;

 }

double bisection(double lower, double upper, std::function<double(double)> potential, bool debug) 
{

   double e_guess = (upper + lower) / 2;
   Schrodinger schro = Schrodinger(e_guess, potential);

   double dx = 0.001;
   double x_first = -5;
   double x_last = 5;
   double y_init = 0;

   double bound = 1;

   auto fxn = [&] (double x, double y) -> double{return schro.query(x, y);};
   
   double accuracy = 0.00000001;

   while(fabs(bound) > accuracy)
   {
       bound = verlet(x_first, x_last, dx, schro, fxn, false);

       if(debug) std::cout << e_guess << std::endl;


       if(bound > 0)
       {
	   lower = e_guess;
	   e_guess = (lower + upper) / 2;
	   schro.setEnergy(e_guess);
       }

       else if(bound < 0) 
       { 
	   upper = e_guess;
	   e_guess = (lower + upper) / 2;
	   schro.setEnergy(e_guess);
       }

   }
   
   return e_guess;
}

double quad_potl(double x)
{
   return 0.5*x*x;
}

int main()
{
    //std::cout << "Hello World \n";
    double dx = 0.001;
    double x_first = -5;
    double x_last = 5;
    double y_init = 0;
    double y_prime_init = 10;

    double e = 0.475;

    auto schro = Schrodinger(e, quad_potl);
    auto fxn = [&] (double x, double y) -> double{return schro.query(x, y);};

    //verlet(x_first, x_last, dx, schro, fxn, true);

    bisection(4.2, 4.9, quad_potl, true);
    //verlet(x_first, x_last, dx, y_init, y_prime_init, sho);
}
