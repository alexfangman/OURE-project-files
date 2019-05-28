#include "stdafx.h"
#include "Solver.h"
#include "Constants.h"

void Solve(const double h, vector<double> &V_r_values, vector<double> &V_t_values, const vector<double> &theta_values)
{
	int EXIT_IT = V_r_values.size();	//Max allowable iterations
	int i = 0;							//Iteration/index value
	double slope_r, slope_t;			//Weighted slope for radial and theta components of velocity (k1+2k2+2k3+k4)

	while (i <= (EXIT_IT-2))
	{
		//Calculate weighted slopes 
		CalculateSlope(h, V_r_values[i], V_t_values[i], slope_r, slope_t, theta_values[i]);

		//Iterate to next value for both components of velocity (soln of ODE's)
		V_r_values[i+1] = V_r_values[i] + ((1.0/6)*slope_r*h);
		V_t_values[i+1] = V_t_values[i] + ((1.0/6)*slope_t*h);

		//Update index/iteration value
		i++;
	}

	return;
}

void CalculateSlope(const double h, const double x1, const double x2, double &slope_r, double &slope_t, const double theta)
{
	double k1, k2, k3, k4;		//k values for slope_r
	double k11, k22, k33, k44;	//k values for slope_t
	double x1_temp, x2_temp;	//Variables used as midpt predicitions

	//Calculate k1, slope at beginning of interval
	k1 = CalculateFunction1Value(x2);
	k11 = CalculateFunction2Value(x1, x2, theta);

	//Calculate first values of x1 and x2 at midpt
	x1_temp = x1+(k1*(h/2));
	x2_temp = x2+(k11*(h/2));

	//Calculate first set of midpt slopes, k2's
	k2 = CalculateFunction1Value(x2_temp);
	k22 = CalculateFunction2Value(x1_temp, x2_temp, (theta+(h/2)));

	//Use these to calculate 2nd set of midpt predictions
	x1_temp = x1 + (k2*(h/2));
	x2_temp = x2 + (k22*(h/2));

	//Calculate second set of midpt slopes, k3's
	k3 = CalculateFunction1Value(x2_temp);
	k33 = CalculateFunction2Value(x1_temp, x2_temp, (theta+(h/2)));

	//Use these to determine predictions at end of interval
	x1_temp = x1 + (k3*h);
	x2_temp = x2 + (k33*h);

	//Use these to compute endpt slopes
	k4 = CalculateFunction1Value(x2_temp);
	k44 = CalculateFunction2Value(x1_temp, x2_temp, (theta+h));

	//Compute weight slopes
	slope_r = k1 + (2*k2) + (2*k3) + k4;
	slope_t = k11 + (2*k22) + (2*k33) + k44;

	return;
}

double CalculateFunction1Value(const double x2)
{
	//Returns input b/c 1st ODE is x1' = x2
	return x2;
}

double CalculateFunction2Value(const double x1, const double x2, const double theta)
{
	double function_value;
	double a = (gamma-1)/2.0;								    //Defined constant to simplify expression
	double num_inner1, num_inner2, num, den;					//Numerator inner two bracketed terms, num, and denominator

	//Split up eqn for readability
	num_inner1 = (1.0/tan(theta))*(x2-(pow(x1,2)*x2)-pow(x2,3));
	num_inner2 = 2*x1*(1-pow(x1,2)-pow(x2,2));
	num = (a*(num_inner1+num_inner2))-(x1*pow(x2,2));
	den = (a*(pow(x1,2)-1))+(pow(x2,2)*(a+1));

	function_value = num/den;

	return function_value;
}

void CalculateFlowField(const double M_inf, const double T_inf, const double p_inf, const double beta, const vector<double> & V_r_values, const vector<double> & V_t_values, vector<double> &V_nondim, vector<double> &Mach, vector<double> &T, vector<double> &p, vector<double> &rho)
{
	double M_inf_n;					//Normal component upstream
	double T_0;						//Total temperature (constant)
	double p_0, p_02;				//Total pressure before and after shock

	//Intermediate terms when calculating total pressure behind shock	
	double term1_num, term1_den, term1, term2_num, term2_den, term2, term3;


	//Calculate normal component of freestream mach
	M_inf_n = M_inf * sin(beta);

	//Calculate nondimensional velocity magnitude and Mach number at along each ray
	for (int i = 0; i <= (V_nondim.size() - 1); i++)
	{
		V_nondim[i] = sqrt(pow(V_r_values[i], 2) + pow(V_t_values[i], 2));

		//Mach calculated from eqn (13.81)
		Mach[i] = sqrt(2.0 / ((gamma - 1)*(pow(V_nondim[i], -2) - 1)));
	}

	//Calculate total temperature using isentropic relations
	T_0 = T_inf * (1 + (((gamma - 1) / 2.0)*pow(M_inf, 2)));

	//Calculate total pressure behind shock using eqn from NASA, split up for readability
	term1_num = (gamma + 1)*pow(M_inf_n, 2);
	term1_den = ((gamma - 1)*(pow(M_inf_n, 2))) + 2;
	term1 = pow((term1_num / term1_den), (gamma / (gamma - 1)));
	term2_num = gamma + 1;
	term2_den = (2 * gamma*pow(M_inf_n, 2)) - (gamma - 1);
	term2 = pow((term2_num / term2_den), (1 / (gamma - 1)));

	//Total pressure before shock
	term3 = 1 + (((gamma - 1) / 2)*pow(M_inf, 2));
	p_0 = p_inf * pow(term3, (gamma / (gamma - 1)));

	//Total pressure behind shock
	p_02 = p_0 * term1*term2;

	//Use isentropic relations to solve for rest of flow field along rays
	for (int i = 0; i <= (V_nondim.size() - 1); i++)
	{
		T[i] = T_0 / (1 + (((gamma - 1) / 2)*pow(Mach[i], 2)));
		term1 = 1 + (((gamma - 1) / 2)*pow(Mach[i], 2));
		p[i] = p_02 / pow(term1, (gamma / (gamma - 1)));
		//Eqn of state to find density
		rho[i] = p[i] / (R*T[i]);
	}

	return;
}