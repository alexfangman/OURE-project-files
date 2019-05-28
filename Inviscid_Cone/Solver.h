#pragma once

#ifndef SOLVER_H
#define SOLVER_H

#include "stdafx.h"
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "Constants.h"

using namespace std;

//Function to solve using RK4 method
void Solve(const double h, vector<double> &V_r_values, vector<double> &V_t_values, const vector<double> &theta_values);

//Function to calculate weighted slope for RK4 method (k1+2k2+2k3+k4)
void CalculateSlope(const double h, const double x1, const double x2, double &slope_r, double &slope_t, const double theta);

//Function to evaluate 1st ODE function x_1' = x_2 where V_r = x_1 and dV_r/dtheta = x_2
double CalculateFunction1Value(const double x2);

//Function to evaluate 2nd ODE function x_2' = "eqn"
double CalculateFunction2Value(const double x1, const double x2, const double theta);

//Function to solve for rest of inviscid shock region using isentropic relations
void CalculateFlowField(const double M_inf, const double T_inf, const double p_inf, const double beta, const vector<double> & V_r_values, const vector<double> & V_t_values, vector<double> &V_nondim, vector<double> &Mach, vector<double> &T, vector<double> &p, vector<double> &rho);

#endif 