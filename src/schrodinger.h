#include<functional>
#ifndef SD_H
#define SD_H

class Schrodinger
{
   private:
      double E;
      std::function<double(double)> potential;

   public:
      Schrodinger(double energy, std::function<double(double)> potential_fxn);
      double query(double x_query, double y_query);
      double initial_conditions(double x_query);
};


#endif
