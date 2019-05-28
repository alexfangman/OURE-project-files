#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>
#include <cmath>

constexpr double gamma = 1.4;
constexpr double PI = 3.1415926535897932384626433;
constexpr double deg2rad = 3.1415926535897932384626433 / 180;
constexpr double rad2deg = 180 / 3.1415926535897932384626433;
constexpr double R = 287.058;
constexpr double CONV_CRIT = 0.00001;			//Convergence criteria for secant method
constexpr int EXIT_IT = 100;		            //Max allowable iterations for secant method

#endif