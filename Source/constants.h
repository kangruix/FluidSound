#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#define M_PI 3.141592653589

#include <cmath>
static const double RHO_WATER = 998.; // density
static const double SIGMA = 0.0726; // surface tension
static const double GAMMA = 1.4; // specific heat ratio
static const double MU = 8.9e-4; // dynamic viscosity of water
static const double GTH = 1.6e6; // thermal damping constant
static const double CF = 1497; // sound speed in water
static const double CA = 343; // sound speed in air
static const double G = 1.0; // TODO: ?
static const double ATM = 101325; // atmospheric pressure
static const double ETA = 0.84; // tan(40 degrees)
static const double ETA_SPLIT = 0.364; // tan(20 degrees)
static const double ETA_ENTRAIN = 1.732; // tan(60 degrees)

#endif // _CONSTANTS_H