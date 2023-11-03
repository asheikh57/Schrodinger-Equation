#include<functional>
#ifndef SD_H
#define SD_H

class Schrodinger
{
   private:
      double E;                                   // Eigen-energy for the Schrodinger equation
      std::function<double(double)> potential;    // Potential function for the Schrodinger equation

   public:

      // Constructor
      Schrodinger(double energy, std::function<double(double)> potential_fxn);

      // Setters
      void setEnergy(double energy);
      void setPotential(std::function<double(double)> potential);

      // Getter
      double getEnergy();

      /* FUNCTION NAME: query
       * 
       * DESCRIPTION: Queries Schrodinger equation for "acceleration" value (second x derivative of solution y) from a value of x and value of y at x
       * @PARAMS
       *     @double x_query - value of x for query
       *     @double y_query - value of y for query, schrodinger depends on x due to potential and y due to the form of the equation itself
       *
       * @RETURN acceleration; the second derivative of the solution y according to the Schrodinger equation
       */ 
      double query(double x_query, double y_query);

      /* FUNCTION NAME: initialConditions
       * 
       * DESCRIPTION: Queries for approximate initial conditions using eigen-energy, potential and an x value
       * @PARAMS
       *     @double x_query - value of x in 
       *
       * @RETURN initial condition (approximate y-value, solution value) at x_query to specify solution and start verlet
       */
      double initialConditions(double x_query);

      ~Schrodinger();
};


#endif
