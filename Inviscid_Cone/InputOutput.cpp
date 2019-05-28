#include "stdafx.h"
#include "InputOutput.h"


void PromptMachConeAngle(double &M_inf, double &cone_angle)
{
	//Prompt for Mach and shock angle
	cout << "Enter freestream Mach Number: ";
	cin >> M_inf;
	cout << "Enter cone half-vertex angle (degrees): ";
	cin >> cone_angle;
	cout << endl;

	//Convert cone angle to radians
	cone_angle = cone_angle * deg2rad;

	return;
}

void OutputValues(const double cone_angle, const double T, const double p, const double rho, const double M, const double beta)
{
	cout << endl << "Cone Angle: " << cone_angle*rad2deg << " degrees" << endl;
	cout << "Shock Angle: " << beta * rad2deg << " degrees" << endl;
	cout << "T_cone: " << T << ", p_cone: " << p << ", rho_cone: " << rho << ", Mach along cone: " << M << endl;

	return;
}