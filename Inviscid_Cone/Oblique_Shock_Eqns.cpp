#include "stdafx.h"
#include "Oblique_Shock_Eqns.h"
#include "Constants.h"

double Calculate2DShockAngle(const double M, const double theta)
{
	int n = 1;									//Current iteration
	double Eps = 1;								//Relative difference in successive guesses
	double f_1, f_2;							//Function values
	double temp;								//Temporary value for swapping
	double beta_1 = theta + (2.0*deg2rad);		//Previous guess for 2D wedge shock angle (radians)
	double beta_2 = theta + (10.0*deg2rad);		//Current guess for shock angle
	double beta_3;								//Next guess for shock angle
	
	while ((n <= EXIT_IT) && (Eps >= CONV_CRIT))
	{
		//Evaluate function at two current points to simplify Secant expression below
		f_1 = Calculate2DWedgeSecantFunction(M, theta, beta_1);
		f_2 = Calculate2DWedgeSecantFunction(M, theta, beta_2);

		//Calculate guess for root based off slope of endpoints
		beta_3 = beta_2 - (f_2 / ((f_2 - f_1) / (beta_2 - beta_1)));

		//Calculate relative difference in successive guesses for root
		Eps = abs(beta_3 - beta_2) / abs(beta_2);

		//Update variables
		temp = beta_2;
		beta_2 = beta_3;
		beta_1 = temp;

		//Increment iteration value
		n++;
	}

	return beta_2;
}

double Calculate2DWedgeSecantFunction(const double M, const double theta, const double beta)
{
	double funct_evaluation;

	funct_evaluation = (2 * (1.0 / tan(beta))*((pow(M, 2)*pow(sin(beta), 2) - 1) / (pow(M, 2)*(gamma + cos(2 * beta)) + 2))) - tan(theta);

	return funct_evaluation;
}

double CalculateDeflectionAngle(const double M, const double B)
{
	double num, den, delta;		//numerator, denominator, shock deflection angle (radians)

	//Split up for readability
	num = (pow(M, 2)*pow(sin(B), 2)) - 1;
	den = (pow(M, 2)*(gamma + cos(2 * B))) + 2;

	delta = atan(2 * (1 / tan(B))*(num / den));
	return delta;
}

double CalculateM2(const double M, const double B, const double theta)
{
	double M_1n, M_2n, M_2;		//Normal component upstream, normal component downstream, Mach downstream
	double num, den;			//numerator, denominator 

	//Normal component upstream
	M_1n = M * sin(B);

	//Split up normal shock relation for readability
	num = 1+((gamma-1)/2)*pow(M_1n,2);
	den = (gamma*pow(M_1n,2))-((gamma-1)/2);
	M_2n = sqrt(num/den);

	//Mach downstream
	M_2 = M_2n/sin(B-theta);

	return M_2;
}

void CalculateComponentsBehindShock(double M, vector<double> &V_r_values, vector<double> &V_t_values, const double B, const double theta)
{
	double V;				//Nondimensional velocity magnitude behind shock
	double inner_power;		//interior power variable, used for readability

	//Split up for readability
	inner_power = (2/((gamma-1)*pow(M,2)))+1;
	V = pow(inner_power,-0.5);

	//Compute components based on geometry, used as initial values for ODE's
	V_r_values[0] = V*cos(B-theta);
	V_t_values[0] = -V*sin(B-theta);

	return;
}