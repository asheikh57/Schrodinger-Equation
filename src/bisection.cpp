
BisectionConfig::BisectionConfig(double dx, double x_first, double x_last, double y_init, double accuracy);

	// Setters
void Bisection::setDX(double dx);
void Bisection::setXFirst(double x_first);
void Bisection::setXLast(double x_last);
void Bisection::setYInit(double y_init);
void Bisection::setAccuracy(double accuracy);

	// Getters
double Bisection::getDX();
double Bisection::getXFirst();
double Bisection::getXLast();
double Bisection::getYInit();
double Bisection::getAccuracy();


double Bisection::bisection_even(double lower, double upper, std::function<double(double)> potential, bool debug);

double Bisection::bisection_odd(double lower, double upper, std::function<double(double)> potential, bool debug);

double Bisection::bisection(double lower, double upper, std::function<double(double)> potential, int parity, bool debug);

double Bisection::bisection(double lower, double upper, std::function<double(double)> potential, bool debug);

BisectionConfig::~BisectionConfig();
