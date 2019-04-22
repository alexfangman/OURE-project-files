#include "stdafx.h"
#include "Secant.h"

double Calc_Prantl_Meyer_Function(const double M, const double function_value, const int option)
{
	double f;						//P-M function value
	double term1, term2, term3;		//Intermediate terms

	term1 = sqrt( (gamma + 1) / (gamma - 1) );
	term2 = atan( sqrt( ( (gamma - 1) / (gamma + 1) )*( pow(M, 2) - 1 ) ) );
	term3 = atan( sqrt( pow(M, 2) - 1) );
	if (option == 1)
	{
		f = (term1*term2) - term3;
	}
	else
	{
		f = (term1*term2) - term3 - function_value;
	}

	return f;
}

double Secant_Method(const double function_value, const double M_i, const double M_i1)
{
	int n = 1;			//Current iteration
	double Eps = 1;		//Relative difference in successive guesses
	double f_1, f_2;	//P-M function values
	double temp;		//Temporary value for swapping

	double M_1 = M_i;	//Previous guess for Mach
	double M_2 = M_i1;	//Current guess for Mach
	double M_3;			//Next guess for Mach

	while ((n <= MAX_IT_SECANT) && (Eps >= CONV_CRIT_SECANT))
	{
		//Evaluate function at two current points to simplify Secant expression below
		f_1 = Calc_Prantl_Meyer_Function(M_1, function_value, 2);
		f_2 = Calc_Prantl_Meyer_Function(M_2, function_value, 2);

		//Calculate guess for root based off slope of endpoints
		M_3 = M_2 - (f_2 / ((f_2 - f_1) / (M_2 - M_1)));

		//Calculate relative difference in successive guesses for root
		Eps = abs(M_3 - M_2) / abs(M_2);

		//Update variables
		temp = M_2;
		M_2 = M_3;
		M_1 = temp;

		//Increment iteration value
		n++;
	}

	return M_2;
}