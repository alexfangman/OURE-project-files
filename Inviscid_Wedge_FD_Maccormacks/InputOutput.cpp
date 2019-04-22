#include "stdafx.h"
#include "InputOutput.h"

void Prompt_Problem_Setup(int &j_max, double &C, double &C_visc)
{
	cout << "Nodes in y-direction: ";
	cin >> j_max;
	cout << "Courant Number (0 < C < 1.0): ";
	cin >> C;
	cout << "Coefficient of Viscosity: ";
	cin >> C_visc;

	return;
}

void Prompt_Initial_Values(double &M_inf, vector< vector<double> > &prim_vars)
{
	double T_inf, p_inf, rho_inf, u_inf;
	double v_inf = 0;

	cout << "Freestream Mach: ";
	cin >> M_inf;
	cout << "Freestream Temperature (K): ";
	cin >> T_inf;
	cout << "Freestream Pressure (Pa): ";
	cin >> p_inf;
	cout << "Freestream Density (kg/m^3): ";
	cin >> rho_inf;

	//Calculate velocity in x-direction
	u_inf = M_inf * sqrt(gamma*R*T_inf);

	//Fill in primitive variable vector w/ freestream values
	for (int i = 0; i <= (prim_vars.size() - 1); i++)
	{
		prim_vars[i][0] = u_inf;
		prim_vars[i][1] = v_inf;
		prim_vars[i][2] = p_inf;
		prim_vars[i][3] = T_inf;
		prim_vars[i][4] = rho_inf;
	}

	return;
}

void Write_Header(ofstream& txt_file)
{
	//Open file
	txt_file.open("flow_field.txt");

	//Write headers
	txt_file << setw(15) << "Iteration" << setw(15)<< "x-coord" << setw(15) << "y-coord" << setw(15) << "u" << setw(15) << "v" << setw(15) << "M" << setw(15) << "p" << setw(15) << "T" << setw(15) << "rho" << "\n";

	//Close file
	txt_file.close();

	return;
}

void WriteSoln(const vector< vector<double> > &prim_vars, ofstream& txt_file, const double xi, const double j_max, const vector<double> &eta_values, const int iteration)
{
	double local_ht;					//Local height of domain at current xi-location
	double local_y_s;					//Local height of bottom surface in physical domain

	vector<double> y(j_max, 0.0);		//y-coordinates
	vector<double> M(j_max, 0.0);		//Mach number

	//Calculate local height of domain
	local_ht = Calc_Local_Ht(xi);

	//Calculate local height of bottom surface
	local_y_s = Calc_Local_y_s(xi);

	//Calculate y-coordinates by transforming them from computational domain eta values
	for (int i = 0; i <= (j_max-1); i++)
	{
		y[i] = (eta_values[i]*local_ht) + local_y_s;
	}

	//Calculate Mach number (V/a)
	for (int i = 0; i <= (j_max - 1); i++)
	{
		M[i] = pow(( pow(prim_vars[i][0], 2) + pow(prim_vars[i][1], 2) ), 0.5) / sqrt(gamma * R * prim_vars[i][3]);
	}

	//Write variables to file: |x-coord|y-coord|u|v|M|p|T|rho|   
	txt_file.open("flow_field.txt", fstream::app);
	for (int i = 0; i <= (j_max - 1); i++)
	{
		txt_file << setprecision(10) << setw(15) << iteration << setw(15) << xi << setw(15) << y[i] << setw(15) << prim_vars[i][0] << setw(15) << prim_vars[i][1] << setw(15) << M[i] << setw(15) << prim_vars[i][2] << setw(15) << prim_vars[i][3] << setw(15) << prim_vars[i][4] << "\n";
	}
	txt_file.close();

	return;
}

void Print_1D_Vector(const vector<double> &V, const int iteration)
{
	cout << "iteration: " << iteration << endl;

	for (int i = 0; i <= V.size() - 1; i++)
	{
		cout << V[i] << "                      ";
	}
	cout << endl;

	return;
}

void Print_2D_Vector(const vector< vector<double> > &V, const int iteration)
{
	cout << "iteration: " << iteration << endl;

	for (int i = 0; i <= V.size() - 1; i++)
	{
		for (int j = 0; j <= V[0].size() - 1; j++)
		{
			cout << V[i][j] << "                      ";
		}
		cout << endl;
	}
	cout << endl;

	return;
}