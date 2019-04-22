#pragma once

#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "Maccormacks.h"
#include <iomanip>
#include <cmath>
#include <math.h>
#include "Constants.h"

using namespace std;

//Function to prompt for number of nodes in y-direction, Courant number, and coefficient of viscosity
void Prompt_Problem_Setup(int &j_max, double &C, double &C_visc);

//Function to prompt for freestream values
void Prompt_Initial_Values(double &M_inf, vector< vector<double> > &prim_vars);

//Function to write header in soln txt file
void Write_Header(ofstream& txt_file);

//Function to write values to .txt file
void WriteSoln(const vector< vector<double> > &prim_vars, ofstream& txt_file, const double xi, const double j_max, const vector<double> &eta_values, const int iteration);

//Function to output 1D vector to console window for checking values
void Print_1D_Vector(const vector<double> &V, const int iteration);

//Function to output 2D vector to console window for checking values
void Print_2D_Vector(const vector< vector<double> > &V, const int iteration);

#endif 