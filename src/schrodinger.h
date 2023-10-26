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

      void setEnergy(double energy);
      void setPotential(std::function<double(double)> potential);

      double getEnergy();

      double query(double x_query, double y_query);
      double initialConditions(double x_query);

      ~Schrodinger();
};


#endif
