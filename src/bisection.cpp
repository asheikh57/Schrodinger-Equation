#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

#include "bisection.h"
#include "schrodinger.h"
#include "verlet.h"


Bisection::Bisection(double dx, double x_first, double x_last, double y_init, double accuracy)
{

    //std::cout.unsetf(std::ios::floatfield);
    //std::cout << std::setprecision(128);
    this->dx = dx;
    this->x_first = x_first;
    this->x_last = x_last;
    this->y_init = y_init;
    this->accuracy = accuracy;
}

// Setters
void Bisection::setDX(double dx)
{
    this->dx = dx;
}

void Bisection::setXFirst(double x_first)
{
    this->x_first = x_first;
}

void Bisection::setXLast(double x_last)
{
    this->x_last = x_last;
}

void Bisection::setYInit(double y_init)
{
    this->y_init = y_init;
}

void Bisection::setAccuracy(double accuracy)
{
    this->accuracy = accuracy;
}

// Getters
double Bisection::getDX()
{
    return dx;
}

double Bisection::getXFirst()
{
    return x_first;
}

double Bisection::getXLast()
{
    return x_last;
}

double Bisection::getYInit()
{
    return y_init;
}

double Bisection::getAccuracy()
{
    return accuracy;
}


double Bisection::bisectionEven(double lower, double upper, std::function<double(double)> potential, bool debug)
{

   double e_guess = (upper + lower) / 2;
   double e_guess_prior = e_guess - accuracy - 1;
   Schrodinger schro = Schrodinger(e_guess, potential);

   double bound = 1;

   auto fxn = [&] (double x, double y) -> double{return schro.query(x, y);};
   
   while(fabs(e_guess - e_guess_prior) > accuracy / 2)
   {
       bound = verlet(x_first, x_last, dx, schro, fxn, false);

       if(debug) std::cout << e_guess << std::endl;


       if(bound > 0)
       {
	   lower = e_guess;
	   e_guess_prior = e_guess;
	   e_guess = (lower + upper) / 2;
	   schro.setEnergy(e_guess);
       }

       else if(bound < 0) 
       { 
	   upper = e_guess;
	   e_guess_prior = e_guess;
	   e_guess = (lower + upper) / 2;
	   schro.setEnergy(e_guess);
       }

   }
   
   return e_guess;
}

double Bisection::bisectionOdd(double lower, double upper, std::function<double(double)> potential, bool debug)
{

   double e_guess = (upper + lower) / 2;
   double e_guess_prior = e_guess - accuracy - 1;
   Schrodinger schro = Schrodinger(e_guess, potential);

   double bound = 1;

   auto fxn = [&] (double x, double y) -> double{return schro.query(x, y);};
   

   while(fabs(e_guess - e_guess_prior) > accuracy / 2)
   {
       bound = verlet(x_first, x_last, dx, schro, fxn, false);

       if(debug) std::cout << e_guess << std::endl;


       if(bound < 0)
       {
	   lower = e_guess;
	   e_guess_prior = e_guess;
	   e_guess = (lower + upper) / 2;
	   schro.setEnergy(e_guess);
       }

       else if(bound > 0) 
       { 
	   upper = e_guess;
	   e_guess_prior = e_guess;
	   e_guess = (lower + upper) / 2;
	   schro.setEnergy(e_guess);
       }

   }
   
   return e_guess;
}

double Bisection::bisection(double lower, double upper, std::function<double(double)> potential, int parity, bool debug)
{
    if(parity) return bisectionOdd(lower, upper, potential, debug);
    else{ return bisectionEven(lower, upper, potential, debug); }
}

double Bisection::bisection(double lower, double upper, std::function<double(double)> potential, bool debug)
{
   double upper_copy = upper;
   double lower_copy = lower;

   double e_guess = (upper_copy + lower_copy) / 2;
   Schrodinger schro = Schrodinger(e_guess, potential);

   double bound_odd, bound_even;

   auto fxn = [&] (double x, double y) -> double{return schro.query(x, y);};
   

   if(debug) std::cout << "\n=== ODD BISECTION ===\n" << std::endl;

   for(int i = 0; i < 50; i++)
   {
       bound_odd = verlet(x_first, x_last, dx, schro, fxn, false);

       if(debug) std::cout << e_guess << std::endl;


       if(bound_odd < 0)
       {
	   lower_copy = e_guess;
	   e_guess = (lower_copy + upper_copy) / 2;
	   schro.setEnergy(e_guess);
       }

       else if(bound_odd > 0) 
       { 
	   upper_copy = e_guess;
	   e_guess = (lower_copy + upper_copy) / 2;
	   schro.setEnergy(e_guess);
       }

   }


   upper_copy = upper;
   lower_copy = lower;

   if(debug) std::cout << "\n=== EVEN BISECTION ===\n" << std::endl;
   
   for(int i = 0; i < 50; i++)
   {
       bound_even = verlet(x_first, x_last, dx, schro, fxn, false);

       if(debug) std::cout << e_guess << std::endl;


       if(bound_even > 0)
       {
	   lower_copy = e_guess;
	   e_guess = (lower_copy + upper_copy) / 2;
	   schro.setEnergy(e_guess);
       }

       else if(bound_even < 0) 
       { 
	   upper_copy = e_guess;
	   e_guess = (lower_copy + upper_copy) / 2;
	   schro.setEnergy(e_guess);
       }

   }

   if(bound_odd < bound_even)
   {
       if(debug) std::cout << "\n=== CHOSEN ODD BISECTION ===\n" << std::endl;
       return bisectionOdd(lower, upper, potential, debug);
   }

   else
   {
       if(debug) std::cout << "\n=== CHOSEN EVEN BISECTION ===\n" << std::endl; 
       return bisectionEven(lower, upper, potential, debug);
   }
   
}

Bisection::~Bisection() {};
