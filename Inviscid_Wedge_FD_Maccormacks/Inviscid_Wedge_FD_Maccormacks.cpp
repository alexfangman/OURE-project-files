//Main.cpp

#include "stdafx.h"
#include <iostream>
#include "Maccormacks.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include "Constants.h"
#include "InputOutput.h"
#include "NewtonsNonlinearSolver.h"
 
using namespace std;

int main()
{
	const int NUM_EQNS = 4;					//Number of governing eqn's
	int j_max, i_max;						//Nodes in the y and x-direction
	double M_inf;							//Freestream Mach
	double C, C_visc;						//Courant number and coefficient of artificial viscosity
	double delta_y, delta_xi, delta_eta;    //Step sizes in physical and computational domains (d_x = d_xi)
	double xi = 0.0;						//Current xi value
	int iteration = 0;		     			//Current spacial iteration step
	double local_ht;						//Local height of domain
	ofstream txt_file;						//Stream object for writing soln values to txt file

	//Write header in soln file
	Write_Header(txt_file);

	//Prompt for size of mesh, Courant number, and coefficient of viscosity
	Prompt_Problem_Setup(j_max, C, C_visc);

	//Vector for primitive variables at current step and predictor step [u|v|p|T|rho]
	vector< vector<double> > prim_vars(j_max, vector<double>(5, 0));
	vector< vector<double> > prim_vars_p(j_max, vector<double>(5, 0));

	//Prompt for freestream conditions
	Prompt_Initial_Values(M_inf, prim_vars);

	vector< vector<double> > F_i(j_max, vector<double>(NUM_EQNS, 0));		//F vector for current nodes (i) 
	vector< vector<double> > F_i1(j_max, vector<double>(NUM_EQNS, 0));		//F vector for next nodes (i+1)
	vector< vector<double> > G_i(j_max, vector<double>(NUM_EQNS, 0));		//G vector for current nodes(i)
	vector< vector<double> > dF_p(j_max, vector<double>(NUM_EQNS, 0));		//F derivative vector for predictor step
	vector< vector<double> > SF_p(j_max, vector<double>(NUM_EQNS, 0));		//Artificial viscosity vector for predictor step
	vector< vector<double> > F_p(j_max, vector<double>(NUM_EQNS, 0));		//Predicted F vector
	vector< vector<double> > G_p(j_max, vector<double>(NUM_EQNS, 0));		//Predicted G vector
	vector< vector<double> > dF_c(j_max, vector<double>(NUM_EQNS, 0));		//F derivative vector for corrector step
	vector< vector<double> > dF_avg(j_max, vector<double>(NUM_EQNS, 0));	//Average F derivative vector
	vector< vector<double> > SF_c(j_max, vector<double>(NUM_EQNS, 0));		//Artificial viscosity vector for corrector step
	vector<double> d_eta_d_x(j_max, 0.0);		//Vector for d_eta/d_x metric for each node in y-direction at current xi-location

	//Calculate delta_y and delta_eta
	delta_y = ht_domain / (j_max - 1.0);
	delta_eta = 1.0 / (j_max - 1.0);

	//Place eta values into a vector 
	vector<double> eta_values(j_max, 0.0);
	for (int i = 0; i <= (j_max-1); i++)
	{
		eta_values[i] = delta_eta * i;
	}

	//Calculate step size, delta_xi using idea of characteristics (assuming all flow initially 0 degrees w.r.t. horizontal)
	delta_xi = C * ( delta_y / ( tan( asin( 1.0/M_inf ) ) ) );

	//Adjust step size and round conservatively so it stops marching at end of l_domain
	i_max = ceil(l_domain / delta_xi);
	delta_xi = l_domain / i_max;

	//Calculate initial F vector from freestream
	Calc_F_Vector(F_i, prim_vars);

	//Output initial values to soln file
	WriteSoln(prim_vars, txt_file, xi, j_max, eta_values, iteration);

	//Calculate initial G vector from freestream
	Calc_G_Vector(G_i, prim_vars);

	//Update xi value for first iteration
	xi = xi + delta_xi;
	iteration++;

	//March forward in x-direction until soln reaches l_domain
	while (iteration <= i_max)
	{
		//Calculate local height of physical domain
		local_ht = Calc_Local_Ht(xi);

		//Calculate d_eta/d_x metric for each node in y-direction at current xi-location
		Calc_d_eta_d_x_Metrics(d_eta_d_x, eta_values, xi, j_max, local_ht);

		//Boundary condition at top wall nodes (constant flow field)
		ConstBoundaryCond(F_i, F_i1);

		//Calculate current G vector 
		Calc_G_Vector(G_i, prim_vars);

		//-----------------Interior and bottom wall nodes: Maccormack's Method-----------------------------//
		//Predictor step (forward differences for spatial terms)
		dF_p = Calc_Predictor_F_Derivative(F_i, G_i, d_eta_d_x, delta_eta, local_ht, j_max, NUM_EQNS);
		SF_p = Calc_Artificial_Visc(F_i, prim_vars, C_visc, j_max, NUM_EQNS);
		F_p = Calc_Predictor_F(F_i, dF_p, SF_p, delta_xi, j_max, NUM_EQNS);

		//Solve non-linear for predicted prim vars
		Calc_Primitive_Vars(F_p, prim_vars_p, xi);	

		//Calculate G vector using primitive variables
		Calc_G_Vector(G_p, prim_vars_p);

		//Corrector step (rearward spacial differences): 1-sided difference for bottom wall nodes
		dF_c = Calc_Corrector_F_Derivative(F_p, G_p, d_eta_d_x, delta_eta, local_ht, j_max, NUM_EQNS);
		dF_avg = Calc_Avg_Derivative(dF_p, dF_c, j_max, NUM_EQNS);
		SF_c = Calc_Artificial_Visc(F_i, prim_vars_p, C_visc, j_max, NUM_EQNS);
		Advance(F_i, F_i1, dF_avg, SF_c, delta_xi);

		//Calculate primitive variables at next xi-location: change to Newton's Method for nonlinear eqn's 
		Calc_Primitive_Vars(F_i1, prim_vars, xi);

		//Boundary condition correction at bottom wall node: Abbett's Euler wall using Prantl-Meyer expansion
		Euler_Wall_BC(prim_vars, xi);

		//Recalculate F vector w/ corrected primitive variables for next iteration (update is already done b/c prim_vars from F_i1)
		Calc_F_Vector(F_i, prim_vars);

		//Output values to file 
		WriteSoln(prim_vars, txt_file, xi, j_max, eta_values, iteration);

		//Update x-location
		xi = xi + delta_xi;

		//Update iteration counter
		iteration++;
	}

    return 0;
}