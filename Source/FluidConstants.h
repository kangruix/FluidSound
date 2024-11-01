/** (c) 2024 Kangrui Xue
 *
 *  FluidConstants.h
 * 
 */

#ifndef _FLUID_CONSTANTS_H
#define _FLUID_CONSTANTS_H


#define _USE_MATH_DEFINES  // needed for M_PI (in Visual Studio)
#include "math.h"

namespace FluidSound {

#ifdef USE_FLOAT64 
	typedef double REAL;
#else
	typedef float REAL;
#endif

static const REAL RHO_WATER = 998.;		// density of water
static const REAL SIGMA = 0.0726;		// surface tension
static const REAL GAMMA = 1.4;			// gas heat capacity ratio
static const REAL MU = 8.9e-4;			// dynamic viscosity of water
static const REAL GTH = 1.6e6;			// thermal damping constant
static const REAL CF = 1497;			// speed of sound in water
static const REAL CA = 343;				// speed of sound in air
static const REAL G = 1.0;				// 
static const REAL ATM = 101325;			// atmospheric pressure
static const REAL ETA = 0.84;			// tan(40 degrees)
static const REAL ETA_SPLIT = 0.364;	// tan(20 degrees)
static const REAL ETA_ENTRAIN = 1.732;	// tan(60 degrees)

} // namespace FluidSound

#endif // _FLUID_CONSTANTS_H