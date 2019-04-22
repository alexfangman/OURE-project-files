#pragma once

#ifndef SECANT_H
#define SECANT_H

#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include "Maccormacks.h"
#include <cmath>
#include <math.h>
#include "Constants.h"

const double CONV_CRIT_SECANT = 0.000001;
const double MAX_IT_SECANT = 100;

//Function to evaluate and return P-M function value
double Calc_Prantl_Meyer_Function(const double M, const double function_value, const int option);

//Function to calculate Mach from P-M function from 2 initial guesses using Secant method
double Secant_Method(const double function_value, const double M_i, const double M_i1);

#endif