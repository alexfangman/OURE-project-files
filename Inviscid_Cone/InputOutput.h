#pragma once

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include <iostream>
#include "Constants.h"

using namespace std;

//Function to prompt for and read in freestream Mach and cone half-vertex angle
void PromptMachConeAngle(double &M_inf, double &cone_angle);

//Function that outputs quantities of interest
void OutputValues(const double cone_angle, const double T, const double p, const double rho, const double M, const double beta);

#endif