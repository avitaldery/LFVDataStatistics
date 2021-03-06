#ifndef MINIM_H
#define MINIM_H

/**
* Interface to the objective function
*/
class IObjective
{
public:
/**
* A function of one parameter
* @param a The value of the parameter
* @return The objective function evaluated at the given paramenter
*/
virtual double function(double a) = 0;
}; // class IObjective

/**
* A minimizer of a function of one variable, using the golden section line
* search method
*/
class Minimizer
{
	public:
	/**
	* The main minimization routine
	* @param obj       An instance of the objective function object
	* @param amin      The lower bound of the search interval
	* @param amax      The upper bound of the search interval
	* @param tolerance The fraction of the origina interval where the
	*                  minimization stops. The default is 0.0001.
	* @return The value of the parameter at the minimum
	*/
	double minimize(IObjective& obj, double& amin, double& amax, double tolerance = 0.000001);
	double minF;    /**< The last value of the objective function */
}; // class Minimizer
#endif // MINIM_H
