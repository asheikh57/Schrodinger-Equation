#include "schrodinger.h"

#ifndef B_H
#define B_H
class Bisection
{
    private:
        double dx;         // For verlet; specficies dx to be used in verlet
	double x_first;    // For verlet; specifies where verlet begins
	double x_last;     // For verlet; specifies where verlet ends
	double y_init;     // For verlet; initial condition
	double accuracy;   // For bisection; specifies termination condition

    public:

	// Constructor
	Bisection(double dx, double x_first, double x_last, double y_init, double accuracy);

	// Setters
	void setDX(double dx);
	void setXFirst(double x_first);
	void setXLast(double x_last);
	void setYInit(double y_init);
	void setAccuracy(double accuracy);

	// Getters
	double getDX();
	double getXFirst();
	double getXLast();
	double getYInit();
	double getAccuracy();

        /* FUNCTION NAME: bisectionEven
	 * 
	 * DESCRIPTION: Bisection for Schrodinger equation for "even" eigen-energes (i.e. ground state, state "0" would be an "even" eigen-energy)
	 *
	 * @PARAMS:
	 *     @double lower - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @double upper - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @std::function<double(double)> potential - potential function to be used for Schrodinger equation
	 *     @bool debug - boolean allowing function to output values to see how algorithm is performing
	 * 
	 * @Return double e_guess - the eigen-energy that bisection converged to given the parameters above 
	 */
        double bisectionEven(double lower, double upper, std::function<double(double)> potential, bool debug);

        /* FUNCTION NAME: bisectionOdd
	 * 
	 * DESCRIPTION: Bisection for Schrodinger equation for "odd" eigen-energes (i.e. first excited state, state "1", would be an "odd" eigen-energy)
	 *
	 * @PARAMS:
	 *     @double lower - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @double upper - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @std::function<double(double)> potential - potential function to be used for Schrodinger equation
	 *     @bool debug - boolean allowing function to output values to see how algorithm is performing
	 * 
	 * @Return double e_guess - the eigen-energy that bisection converged to given the parameters above 
	 */
        double bisectionOdd(double lower, double upper, std::function<double(double)> potential, bool debug);

        /* FUNCTION NAME: bisection
	 * 
	 * DESCRIPTION: Bisection for Schrodinger equation 
	 *
	 * @PARAMS:
	 *     @double lower - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @double upper - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @std::function<double(double)> potential - potential function to be used for Schrodinger equation
	 *     @int parity - selects between bisectionOdd and bisectionEven
	 *     @bool debug - boolean allowing function to output values to see how algorithm is performing
	 * 
	 * @Return double e_guess - the eigen-energy that bisection converged to given the parameters above 
	 */
        double bisection(double lower, double upper, std::function<double(double)> potential, int parity, bool debug);

        /* FUNCTION NAME: bisection
	 * 
	 * DESCRIPTION: Bisection for Schrodinger equation, automatically determines whether to use bisectionOdd or bisectionEven
	 *
	 * @PARAMS:
	 *     @double lower - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @double upper - lower part of the interval of guesses, i.e. guess the actual value is between lower and upper
	 *     @std::function<double(double)> potential - potential function to be used for Schrodinger equation
	 *     @bool debug - boolean allowing function to output values to see how algorithm is performing
	 * 
	 * @Return double e_guess - the eigen-energy that bisection converged to given the parameters above 
	 */
        double bisection(double lower, double upper, std::function<double(double)> potential, bool debug);

	~Bisection();
};


#endif
