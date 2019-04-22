#pragma once

#ifndef MACCORMACKS_H
#define MACCORMACKS_H

#include "stdafx.h"
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "InputOutput.h"
#include "Constants.h"
#include "Secant.h"

using namespace std;

//Function that calculates components of F flux vector given primative variables
void Calc_F_Vector(vector< vector<double> > &F, const vector< vector<double> > &prim_vars);

//Function that calculates components of G flux vector given primative variables
void Calc_G_Vector(vector< vector<double> > &G, const vector< vector<double> > &prim_vars);

//Function that keeps values constant along top boundary
void ConstBoundaryCond(vector< vector<double> > &F_i, vector< vector<double> > &F_i1);

//Function that updates vectors by swapping the old with the new
void Update_F_Vector(vector< vector<double> > &F_i, vector< vector<double> > &F_i1);

//Function that returns the local height of the domain at given xi-location
double Calc_Local_Ht(const double xi);

//Function that returns the local height of the bottom surface in physical domain at given xi-location
double Calc_Local_y_s(const double xi);

//Function to calculate d_eta_d_x metrics
void Calc_d_eta_d_x_Metrics(vector<double> &d_eta_d_x, const vector<double> &eta_values, const double xi, const double j_max, const double local_ht);

//Function to calculate F derivative from current F values
vector< vector<double> > Calc_Predictor_F_Derivative(const vector< vector<double> > &F_i, const vector< vector<double> > &G_i, const vector<double> &d_eta_d_x, const double delta_eta, const double local_ht, const int j_max, const int NUM_EQNS);

//Function to calculate artificial viscosity, no viscosity for wall nodes
vector< vector<double> > Calc_Artificial_Visc(const vector< vector<double> > &F_i, const vector< vector<double> > &prim_vars, const double C_visc, const int j_max, const int NUM_EQNS);

//Function to calculate predictor F vector
vector< vector<double> > Calc_Predictor_F(const vector< vector<double> > &F_i, const vector< vector<double> > &dF_p, const vector< vector<double> > &SF_p, const double delta_xi, const int j_max, const int NUM_EQNS);

//Function to calculate F derivative from predicted values, using 1-sided difference for wall nodes
vector< vector<double> > Calc_Corrector_F_Derivative(const vector< vector<double> > &F_p, const vector< vector<double> > &G_p, const vector<double> &d_eta_d_x, const double delta_eta, const double local_ht, const int j_max, const int NUM_EQNS);

//Function to calculate average F derivative
vector< vector<double> > Calc_Avg_Derivative(const vector< vector<double> > &dF_p, const vector< vector<double> > &dF_c, const int j_max, const int NUM_EQNS);

//Function to calculate F values at next xi-location
void Advance(const vector< vector<double> > &F_i, vector< vector<double> > &F_i1, const vector< vector<double> > &dF_avg, const vector< vector<double> > &SF_c, const double delta_xi);

//Function to apply Abbett's Euler wall B.C.
void Euler_Wall_BC(vector< vector<double> > &prim_vars, const double xi);

//Function that primitive variables from flux vector, F


#endif 
