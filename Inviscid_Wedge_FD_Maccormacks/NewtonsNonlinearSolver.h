#pragma once

#ifndef NEWTONS_NONLINEAR_SOLVER_H
#define NEWTONS_NONLINEAR_SOLVER_H

#include "stdafx.h"
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "Constants.h"
#include "InputOutput.h"

using namespace std;

//Function that uses Newton's Method to solve nonlinear set of eqn's for primitive variables
//Gaussian elimination with partial pivoting is used to solve the resulting set of linear eqn's
void Calc_Primitive_Vars(const vector< vector<double> > &F, vector< vector<double> > &prim_vars, const double xi);

#endif