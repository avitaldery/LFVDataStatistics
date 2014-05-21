#include "LFVDataStatistics/minim.h"
#include <math.h>


const double GOLD2 = (sqrt(5.0) - 1.0) / 2.0;
const double GOLD1 = 1.0 - GOLD2;

double Minimizer::minimize(IObjective& obj, double& amin, double& amax, double tolerance)
{
	double a0 = amin;
	obj.function(a0);
	double a3 = amax;
	obj.function(a3);
	double a1 = GOLD2 * a0 + GOLD1 * a3, f1 = obj.function(a1);
	double a2 = GOLD1 * a0 + GOLD2 * a3, f2 = obj.function(a2);
	double tol = tolerance * (a3 - a0);

	while (a3 - a0 > tol)
	{
//		printf("(a0,a1,a2,a3)=(%8.7f,%8.7f,%8.7f,%8.7f)    ",a0,a1,a2,a3);
//		printf("(f1,f2)=(%8.7f,%8.7f)\n",f1,f2);

		if (f1 > f2)
		{
			a0 = a1;
			a1 = a2;
			f1 = f2;
			a2 = GOLD1 * a3 + GOLD2 * a1;
			f2 = obj.function(a2);
		}
		else
		{
			a3 = a2;
			a2 = a1;
			f2 = f1;
			a1 = GOLD2 * a2 + GOLD1 * a0;
			f1 = obj.function(a1);
		}
	}
	a1 = (a1 + a2) / 2.0; minF = obj.function(a1);
	//     amin = a0;
	//     amax = a3;

//cout << endl;
	return a1;
}


// To use the following test code, compile this module with -DTESTING

#ifdef TESTING
#include <iostream>

class parabola : public IObjective
{
	public:
	virtual double function(double a)
	{
	return 1.0 + 2.0 * a + 3.0 * a * a;
	}
};

int main(int argc, char* argv[])
{
	parabola p;
	Minimizer minim;
	double mina = minim.minimize(p, -2.0, 2.0, 0.0001);

//std::cout << "Minimum is f=" << minim.minF << " at a=" << mina << std::endl;
}
#endif // TESTING
