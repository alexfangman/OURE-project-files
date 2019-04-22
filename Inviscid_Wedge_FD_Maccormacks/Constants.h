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

constexpr double ZERO_CRIT = 0.001;

//Domain
constexpr double wedge_angle = 5 * (3.1415926535897932384626433 / 180);  //Half-wedge angle (radians)
constexpr double ht_domain = 1;										  //Height of physical domain (m)
constexpr double l_domain = 1.3;										  //Length of physical domain, wedge chord is 1 (m)
constexpr double wedge_location = 0.3;									  //x-location of ramp vertex: corresponds to (0.3,0) coordinate in physical domain

#endif