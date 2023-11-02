#include "schrodinger.h"

#ifndef B_H
#define B_H
class Bisection
{
    private:
        double dx;
	double x_first;
	double x_last;
	double y_init;
	double accuracy;

    public:
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


        double bisectionEven(double lower, double upper, std::function<double(double)> potential, bool debug);

        double bisectionOdd(double lower, double upper, std::function<double(double)> potential, bool debug);

        double bisection(double lower, double upper, std::function<double(double)> potential, int parity, bool debug);

        double bisection(double lower, double upper, std::function<double(double)> potential, bool debug);

	~Bisection();
};


#endif
