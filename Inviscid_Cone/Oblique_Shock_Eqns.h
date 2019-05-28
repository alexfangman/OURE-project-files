#pragma once

#ifndef OBLIQUE_SHOCK_EQNS_H
#define OBLIQUE_SHOCK_EQNS_H

#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include "Constants.h"

using namespace std;

//Function to calculate shock angle for 2D wedge using Mach-Beta-Theta relation (Secant method)
double Calculate2DShockAngle(const double M, const double theta);

//Function to calculate F(beta,theta,M) for finding root F(beta) = 0
double Calculate2DWedgeSecantFunction(const double M, const double theta, const double beta);

//Function to calculate flow deflection angle using Mach-Beta-Theta relation
double CalculateDeflectionAngle(const double M, const double B);

//Function to calculate Mach behind o.s.w.
double CalculateM2(const double M, const double B, const double theta);

//Function to calculate nondimensional velocity and components behind o.s.w.
void CalculateComponentsBehindShock(double M, vector<double> &V_r_values, vector<double> &V_t_values, const double B, const double theta);

#endif 