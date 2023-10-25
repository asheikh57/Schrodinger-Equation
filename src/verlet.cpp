#include<iostream>
#include<iomanip>
#include<functional>

double schrodinger(double psi, double E, std::function<double(double)> potential)
{
   return 2 * (potential(psi) - E) * psi; 
}

double sho(double x)
{
    return -1*x;
}

double verlet_init(double dx, double y_init, double y_prime_init, std::function<double(double)> func)
{
    return y_init + y_prime_init*dx + 0.5*func(y_init)*dx*dx;
}

double verlet_proc(double dx, double psi_curr, double psi_prev, std::function<double(double)> func)
{
    return 2*psi_curr - psi_prev + func(psi_curr)*dx*dx;
}

void verlet(double x_first,
	    double x_last,	
	    double dx, 
	    double y_init, 
	    double y_prime_init, 
	    std::function<double(double)> func)
{
    double x_curr = x_first;
    double y_curr = y_init;

    std::cout << x_curr << "," << y_curr << std::endl;
    

    x_curr += dx;
    y_curr = verlet_init(dx, y_init, y_prime_init, func);

    std::cout << x_curr << "," << y_curr << std::endl;

    double y_last = y_init;
    
    int num_iters = ((x_last - x_first) / dx);
    double y_next;
    
    for(int i = 0; i < num_iters; i++)
    {
	x_curr += dx;
	y_next = verlet_proc(dx, y_curr, y_last, func);

	y_last = y_curr;
	y_curr = y_next;


        std::cout << x_curr << "," << y_curr << std::endl;

    }

    
}

int main()
{
    //std::cout << "Hello World \n";
    double dx = 0.1;
    double x_first = 0;
    double x_last = 100;
    double y_init = 0;
    double y_prime_init = 1;

    verlet(x_first, x_last, dx, y_init, y_prime_init, sho);
}
