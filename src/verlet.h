#ifndef V_H
#define V_H

#include <functional>

/* FUNCTION NAME: verlet_init
 *
 * DESCRIPTION: intializes verlet algorithm (not velocity verlet)
 * @PARAMS:
 *    @double param dx - delta x difference between each x "query" point for the differential equation being solved
 *    @double x_init - x-coordinate of initial condition, t0 for time dependent differential equation
 *    @double y_init - initial condition value, or y-value at x_init for solution function
 *    @y_prime_init - initial condition value, velocity initial condition or initial condition on derivative of solution at x_init
 *    @func - function of "accleration on position", i.e. represents second order differential equation depending only on x (in general time) and y 
 *
 * @RETURN next y value at x_init + dx
 *
 *
 */
double verlet_init(double dx, double x_init, double y_init, double y_prime_init, std::function<double(double, double)> func);



/* FUNCTION NAME: verlet_proc
 *
 * DESCRIPTION: performs verlet algorithm iteration given
 * @PARAMS:
 *    @double param dx - delta x difference between each x "query" point for the differential equation being solved
 *    @double x - the current value of that the algorithm has gotten to, the query point
 *    @double psi_curr - the current value of the approximated numerical solution at "x"
 *    @double psi_prev - the last value of the approximated numerical solution at "x - dx"
 *    @func - function of "accleration on position", i.e. represents second order differential equation depending only on x (in general time) and y 
 *
 * @RETURN next y value at x + dx
 *
 *
 */
double verlet_proc(double dx, double x, double psi_curr, double psi_prev, std::function<double(double, double)> func);

/* FUNCTION NAME: verlet
 *
 * DESCRIPTION: intializes verlet algorithm (not velocity verlet)
 * @PARAMS:
 *    @double x_first - First x query value -- x value of initial condition (e.g. y(x_first) = y_init)
 *    @double x_last - Last x query value, x_last > x_first,
 *    @double param dx - delta x difference between each x "query" point for the differential equation being solved
 *    @double y_init - initial condition value, or y-value at x_init for solution function
 *    @y_prime_init - initial condition value, velocity initial condition or initial condition on derivative of solution at x_init
 *    @func - function of "accleration on position", i.e. represents second order differential equation depending only on x (in general time) and y 
 *
 * @RETURN y value at x_last
 *
 *
 */
double verlet(double x_first,
	    double x_last,	
	    double dx, 
	    double y_init, 
	    double y_prime_init, 
	    std::function<double(double, double)> func);


/* FUNCTION NAME: verlet
 *
 * DESCRIPTION: intializes verlet algorithm (not velocity verlet)
 * @PARAMS:
 *    @double x_first - First x query value -- x value of initial condition (e.g. y(x_first) = y_init)
 *    @double x_last - Last x query value, x_last > x_first,
 *    @double param dx - delta x difference between each x "query" point for the differential equation being solved
 *    @double y_init - initial condition value, or y-value at x_init for solution function
 *    @y_prime_init - initial condition value, velocity initial condition or initial condition on derivative of solution at x_init
 *    @func - function of "accleration on position", i.e. represents second order differential equation depending only on x (in general time) and y 
 *
 * @RETURN y value at x_last
 *
 *
 */
double verlet(double x_first,
	    double x_last,	
	    double dx, 
	    Schrodinger schro,
	    bool log);
#endif
