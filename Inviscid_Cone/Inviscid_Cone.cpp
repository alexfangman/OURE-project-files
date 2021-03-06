// Inviscid_Cone.cpp : 

#include "stdafx.h"
#include "Oblique_Shock_Eqns.h"
#include <iostream>
#include "Solver.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include "Constants.h"
#include "InputOutput.h"

using namespace std;

int main()
{
	const double h = -.001*deg2rad; //Step size in theta (degrees)
	double M_inf;					//Freestream Mach
	double M_2;						//Mach Number behind shock
	double cone_angle;				//Half-vertex angle of cone
	double beta_1, beta_2, beta_3;	//Previous, current, & next cone shock angle for secant method 
	int n = 1;						//Iteration number
	double delta;					//Flow deflection angle
	double f_1 = 1, f_2 = 1;		//Variables that store function evaluations for secant method (function is theta-comp. of velocity on surface of cone)
	int ARRAY_SIZE;					//Size of flow variable vectors
	double temp;					//Temporary variable used for swapping values when updating

	vector<double> V_r_values;		//Radial components of velocity (nondim)
	vector<double> V_t_values;		//Theta components of velocity (nondim)
	vector<double> theta_values;	//Vector of shock region theta values (rays along which flow field is constant)

	//Freestream values
	const double T_inf = 227.13;	//Freestream temperature (K)
	const double p_inf = 1090.16;	//Freestream pressure (Pa)

	//Prompt for Mach and shock angle
	PromptMachConeAngle(M_inf, cone_angle);

	//------------------------------------This should be changed---------------------------------------------------//
	//Initial guesses for beta are: (cone half-vertex angle + 2.0 degrees and 2D wedge o.s.w. soln - 5 deg.)
	beta_1 = cone_angle+(2.0*deg2rad);
	
	//Use Mach-Beta-Theta relation to calculate second guess for shock angle
	beta_2 = Calculate2DShockAngle(M_inf, cone_angle) - 5.0*deg2rad;
	//-------------------------------------------------------------------------------------------------------------//

	while ((n <= EXIT_IT) && (abs(f_2) >= CONV_CRIT))
	{
		//Evaluate beta_1 on first iteration for secant method, beta_2 will always be evaluated and updated to be beta_1 before next iteration
		if (n == 1)
		{
			//Calculate flow deflection angle and Mach behind shock using o.s.w. analysis
			delta = CalculateDeflectionAngle(M_inf, beta_1);
			M_2 = CalculateM2(M_inf, beta_1, delta);

			//Appropriate size for allocating array, rounding up
			ARRAY_SIZE = 1.0 + abs((ceil((beta_1-cone_angle) / h)));

			//Vectors to store radial and theta components (nondim) along rays, initialized to zero
			V_r_values.resize(ARRAY_SIZE, 0.0);
			V_t_values.resize(ARRAY_SIZE, 0.0);

			//Calculate nondimensional velocity magnitude behind shock from Anderson eqn (13.81) to be used as initial values for ODE's
			CalculateComponentsBehindShock(M_2, V_r_values, V_t_values, beta_1, delta);

			//Vector of theta data values
			theta_values.resize(ARRAY_SIZE, 0.0);
			for (int i = 0; i < ARRAY_SIZE; i++)
			{
				theta_values[i] = beta_1 + (i*h);
			}

			//Solve the 2 first-order ODE's until cone_angle with initial conditions V_r and V_theta from initial shock guess
			Solve(h, V_r_values, V_t_values, theta_values);

			//Find V_t on surface of cone (last entry in vector) 
			f_1 = V_t_values[ARRAY_SIZE-1];
		}

		//Calculate flow deflection angle and Mach behind shock using o.s.w. analysis
		delta = CalculateDeflectionAngle(M_inf, beta_2);
		M_2 = CalculateM2(M_inf, beta_2, delta);

		//Appropriate size for allocating array, rounding up
		ARRAY_SIZE = 1.0 + abs((ceil((beta_2-cone_angle) / h)));

		//Vectors to store radial and theta components (nondim) along rays, initialized to zero
		V_r_values.resize(ARRAY_SIZE, 0.0);
		V_t_values.resize(ARRAY_SIZE, 0.0);

		//Calculate nondimensional velocity magnitude behind shock from Anderson eqn (13.81) to be used as initial values for ODE's
		CalculateComponentsBehindShock(M_2, V_r_values, V_t_values, beta_2, delta);

		//Vector of theta data values
		theta_values.resize(ARRAY_SIZE, 0.0);
		for (int i = 0; i < ARRAY_SIZE; i++)
		{
			theta_values[i] = beta_2 + (i*h);
		}

		//Solve the 2 first-order ODE's until cone_angle with initial conditions V_r and V_theta from initial shock guess
		Solve(h, V_r_values, V_t_values, theta_values);

		//Find V_t on surface of cone (last entry in vector) 
		f_2 = V_t_values[ARRAY_SIZE - 1];

		//Calculate guess for root based off slope of endpoints
		beta_3 = beta_2 - (f_2 / ((f_2 - f_1) / (beta_2 - beta_1)));

//-------------------------------------------------------------------------------------------------//
cout << "Iteration: " << n << "     " << "f_1: " << f_1 << "     f_2: " << f_2 << endl;
cout << "Shock angle after step: " << beta_3 * rad2deg << endl << endl;
//-------------------------------------------------------------------------------------------------//		

		//Update variables
		temp = beta_2;
		beta_2 = beta_3;
		beta_1 = temp;
		f_1 = f_2;

		//Increment iteration number
		n++;
	}

	//Vectors to store nondim velocity magnitude and Mach number along rays
	vector<double> V_nondim(ARRAY_SIZE, 0.0);
	vector<double> Mach(ARRAY_SIZE, 0.0);

	//Vectors to store thermodynamic quantities
	vector<double> T(ARRAY_SIZE, 0.0);
	vector<double> p(ARRAY_SIZE, 0.0);
	vector<double> rho(ARRAY_SIZE, 0.0);

	//Calculate rest of flow field using isentropic relations
	CalculateFlowField(M_inf, T_inf, p_inf, beta_3, V_r_values, V_t_values, V_nondim, Mach, T, p, rho);

	//Output results to text file for importing into MATLAB for plotting [theta|pressure|temperature|density]
	ofstream myfile;
	myfile.open("cone_flow_field.txt");
	for (int i = 0; i <= (V_nondim.size() - 1); i++)
	{
		myfile << setprecision(10) << theta_values[i]*rad2deg << setw(15) << p[i] << setw(15) << T[i] << setw(15) << rho[i] << "\n";
	}
	myfile.close();

	//Output shock angle and important surface quantities
	OutputValues(cone_angle, T[ARRAY_SIZE - 1], p[ARRAY_SIZE - 1], rho[ARRAY_SIZE - 1], Mach[ARRAY_SIZE - 1], beta_3);
		
cin >> beta_1;

    return 0;
}