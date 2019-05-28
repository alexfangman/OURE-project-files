#include "stdafx.h"
#include "Maccormacks.h"

void Calc_F_Vector(vector< vector<double> > &F, const vector< vector<double> > &prim_vars)
{
	
	for (int i = 0; i <= F.size() - 1; i++)
	{
		for (int j = 0; j <= F[0].size() - 1; j++)
		{
			if (j == 0)
			{
				//F1 = rho * u;
				F[i][j] = prim_vars[i][4] * prim_vars[i][0];
			}
			else if (j == 1)
			{
				//F2 = (rho * pow(u, 2) ) + p;
				F[i][j] = (prim_vars[i][4] * pow(prim_vars[i][0], 2) ) + prim_vars[i][2];
			}
			else if (j == 2)
			{
				//F3 = rho * u * v;
				F[i][j] = prim_vars[i][4] * prim_vars[i][0] * prim_vars[i][1];
			}
			else
			{
				//F4 = ( gamma / ( gamma - 1.0 ) )*(p * u) + (rho * u)*( ( pow(u,2) + pow(v, 2)) / 2 );
				F[i][j] = (gamma / (gamma - 1.0))*(prim_vars[i][2] * prim_vars[i][0]) + (prim_vars[i][4] * prim_vars[i][0])*((pow(prim_vars[i][0], 2) + pow(prim_vars[i][1], 2)) / 2);
			}
		}
	}

	return;
}

void Calc_G_Vector(vector< vector<double> > &G, const vector< vector<double> > &prim_vars)
{

	for (int i = 0; i <= G.size() - 1; i++)
	{
		for (int j = 0; j <= G[0].size() - 1; j++)
		{
			if (j == 0)
			{
				//G1 = rho * v;
				G[i][j] = prim_vars[i][4] * prim_vars[i][1];
			}
			else if (j == 1)
			{
				//G2 = rho * u * v;
				G[i][j] = prim_vars[i][4] * prim_vars[i][0] * prim_vars[i][1];
			}
			else if (j == 2)
			{
				//G3 = (rho * pow(v, 2)) + p;
				G[i][j] = (prim_vars[i][4] * pow(prim_vars[i][1], 2) ) + prim_vars[i][2];
			}
			else
			{
				//G4 = (gamma / (gamma - 1.0))*(p * v) + (rho * v)*((pow(u, 2) + pow(v, 2)) / 2);
				G[i][j] = (gamma / (gamma - 1.0))*(prim_vars[i][2] * prim_vars[i][1]) + (prim_vars[i][4] * prim_vars[i][1])*((pow(prim_vars[i][0], 2) + pow(prim_vars[i][1], 2)) / 2);
			}
		}
	}

	return;
}

void ConstBoundaryCond(vector< vector<double> > &F_i, vector< vector<double> > &F_i1)
{
	int i = (F_i.size() - 1);				//Top wall index value

	//Keep values constant at top boundary
	for (int j = 0; j <= F_i[0].size() - 1; j++)
	{
		F_i1[i][j] = F_i[i][j];
	}
	
	return;
}

void Update_F_Vector(vector< vector<double> > &F_i, vector< vector<double> > &F_i1)
{
	for (int i = 0; i <= F_i.size() - 1; i++)
	{
		for (int j = 0; j <= F_i[0].size() - 1; j++)
		{
			F_i[i][j] = F_i1[i][j];
		}
	}
	return;
}

double Calc_Local_Ht(const double xi)
{
	double local_ht;

	if (xi <= wedge_location)
	{
		local_ht = ht_domain;
	}
	else
	{
		local_ht = ht_domain - ((xi - wedge_location)*tan(wedge_angle));
	}

	return local_ht;
}

double Calc_Local_y_s(const double xi)
{
	double local_y_s;

	if (xi <= wedge_location)
	{
		local_y_s = 0.0;
	}
	else
	{
		local_y_s = (xi - wedge_location)*tan(wedge_angle);
	}

	return local_y_s;
}

void Calc_d_eta_d_x_Metrics(vector<double> &d_eta_d_x, const vector<double> &eta_values, const double xi, const double j_max, const double local_ht)
{
	if (xi <= wedge_location)
	{
		for (int i = 0; i <= (j_max - 1); i++)
		{
			d_eta_d_x[i] = 0;
		}
	}
	else
	{
		for (int i = 0; i <= (j_max - 1); i++)
		{
			d_eta_d_x[i] = (eta_values[i]-1)*(tan(wedge_angle)/local_ht);
		}
	}

	return;
}

vector< vector<double> > Calc_Predictor_F_Derivative(const vector< vector<double> > &F_i, const vector< vector<double> > &G_i, const vector<double> &d_eta_d_x, const double delta_eta, const double local_ht, const int j_max, const int NUM_EQNS)
{
	vector< vector<double> > dF_p(j_max, vector<double>(NUM_EQNS, 0));		//F derivative vector for predictor step

	//Spatial derivatives using forward differences: only calculating interior nodes (outside for loop begins at 1)
	for (int i = 0; i <= (dF_p.size() - 2); i++)
	{
		for (int j = 0; j <= (dF_p[0].size() - 1); j++)
		{
			dF_p[i][j] = ( d_eta_d_x[i] ) * ( ( F_i[i][j] - F_i[i+1][j] ) / delta_eta) + (1.0 / local_ht) * ( ( G_i[i][j] - G_i[i+1][j] ) / delta_eta );
		}
	}

	return dF_p;
}

vector< vector<double> > Calc_Artificial_Visc(const vector< vector<double> > &F_i, const vector< vector<double> > &prim_vars, const double C_visc, const int j_max, const int NUM_EQNS)
{
	vector< vector<double> > SF(j_max, vector<double>(NUM_EQNS, 0));		//Artificial viscosity vector for predictor step

	for (int i = 1; i <= (SF.size() - 2); i++)
	{
		for (int j = 0; j <= (SF[0].size() - 1); j++)
		{
			SF[i][j] = C_visc * ( (prim_vars[i+1][2] - 2*prim_vars[i][2] + prim_vars[i-1][2]) / (prim_vars[i+1][2] + 2*prim_vars[i][2] + prim_vars[i - 1][2]) ) * (F_i[i+1][j] - 2*F_i[i][j] + F_i[i-1][j]);
		}
	}

	return SF;
}

vector< vector<double> > Calc_Predictor_F(const vector< vector<double> > &F_i, const vector< vector<double> > &dF_p, const vector< vector<double> > &SF_p, const double delta_xi, const int j_max, const int NUM_EQNS)
{
	vector< vector<double> > F_p(j_max, vector<double>(NUM_EQNS, 0));		//Predicted F vector

	for (int i = 0; i <= (F_p.size() - 1); i++)
	{
		for (int j = 0; j <= (F_p[0].size() - 1); j++)
		{
			F_p[i][j] =  F_i[i][j] + ( dF_p[i][j] * delta_xi ) + SF_p[i][j];
		}
	}

	return F_p;
}

vector< vector<double> > Calc_Corrector_F_Derivative(const vector< vector<double> > &F_p, const vector< vector<double> > &G_p, const vector<double> &d_eta_d_x, const double delta_eta, const double local_ht, const int j_max, const int NUM_EQNS)
{
	int wall_node = 0;														//Index for wall node
	vector< vector<double> > dF_c(j_max, vector<double>(NUM_EQNS, 0));		//F derivative vector for predictor step

	//Spatial derivatives using rearward differences: only calculating interior nodes (outside for loop begins at 1)
	for (int i = 1; i <= (dF_c.size() - 2); i++)
	{
		for (int j = 0; j <= (dF_c[0].size() - 1); j++)
		{
			dF_c[i][j] = (d_eta_d_x[i]) * ((F_p[i-1][j] - F_p[i][j]) / delta_eta) + (1.0 / local_ht) * ((G_p[i-1][j] - G_p[i][j]) / delta_eta);
		}
	}

	//Wall node using 1-sided differences
	for (int j = 0; j <= (dF_c[0].size() - 1); j++)
	{
		dF_c[wall_node][j] = (d_eta_d_x[wall_node]) * ((F_p[wall_node][j] - F_p[wall_node + 1][j]) / delta_eta) + (1.0 / local_ht) * ((G_p[wall_node][j] - G_p[wall_node + 1][j]) / delta_eta);
	}

	return dF_c;
}

vector< vector<double> > Calc_Avg_Derivative(const vector< vector<double> > &dF_p, const vector< vector<double> > &dF_c, const int j_max, const int NUM_EQNS)
{
	vector< vector<double> > dF_avg(j_max, vector<double>(NUM_EQNS, 0));		//Average F derivative 

	//Average derivative for all nodes except top boundary b/c constant B.C.
	for (int i = 0; i <= (dF_avg.size() - 2); i++)
	{
		for (int j = 0; j <= (dF_avg[0].size() - 1); j++)
		{
			dF_avg[i][j] = 0.5*( dF_p[i][j] + dF_c[i][j] );
		}
	}

	return dF_avg;
}

void Advance(const vector< vector<double> > &F_i, vector< vector<double> > &F_i1, const vector< vector<double> > &dF_avg, const vector< vector<double> > &SF_c, const double delta_xi)
{
	for (int i = 0; i <= (F_i1.size() - 2); i++)
	{
		for (int j = 0; j <= (F_i1[0].size() - 1); j++)
		{
			F_i1[i][j] = F_i[i][j] + ( dF_avg[i][j] * delta_xi ) + SF_c[i][j];
		}
	}
}

void Euler_Wall_BC(vector< vector<double> > &prim_vars, const double xi)
{
	double phi;						//Expansion angle correction (positive or negative depending on into/out of surface)
	double psi;						//Initially calculated angle on ramp
	double M_cal, M_act;			//Initially calculated and corrected Mach
	double f_cal, f_act;			//Initial and corrected P-M f(M)
	double p_act, T_act, rho_act;	//Corrected pressure, temperature, density
	double term1, term2, term3;		//Intermediate terms for P-M function
	double num, den;				//Intermediate  numerator and denominator
	double v_act;					//Corrected y-velocity

	//Calculate expansion angle correction (must determine if along flat edge or ramp for local horizontal)
	if (xi < wedge_location)
	{
		//tan^-1(v/u)
		phi = atan(prim_vars[0][1] / prim_vars[0][0]);
	}
	else
	{
		psi = atan( abs(prim_vars[0][1]) / prim_vars[0][0]);
		phi = psi - wedge_angle;
	}

	//Calculate initial Mach at wall
	M_cal = sqrt( pow(prim_vars[0][0], 2) + pow(prim_vars[0][1], 2) ) / sqrt (gamma * R * prim_vars[0][3]);

	//Calculate initial Prantl-Meyer
	f_cal = Calc_Prantl_Meyer_Function(M_cal, 0, 1);

	//Calculate actual Prantl-Meyer
	f_act = f_cal + phi;

	//Solve for actual Mach using Secant method
	M_act = Secant_Method(f_act, M_cal, (M_cal+.5));

	//Calculate corrected pressure
	num = 1 + (((gamma - 1) / 2) * pow(M_cal, 2) );
	den = 1 + (((gamma - 1) / 2) * pow(M_act, 2));
	p_act = prim_vars[0][2]*pow((num / den), (gamma / (gamma - 1)));
	//Calculate corrected Temperature
	T_act = prim_vars[0][3] * ( num / den );
	//Calculate corrected density
	rho_act = p_act / (R * T_act);
	//Calculate corrected y-velocity, v, so tangent to wall
	if (xi < wedge_location)
	{
		v_act = 0;
	}
	else
	{
		v_act = prim_vars[0][0]*tan(wedge_angle);
	}

	//Store corrected values in prim_vars vector
	prim_vars[0][1] = v_act;
	prim_vars[0][2] = p_act;
	prim_vars[0][3] = T_act;
	prim_vars[0][4] = rho_act;

	return;
}