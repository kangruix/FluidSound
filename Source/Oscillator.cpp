/** (c) 2024 Kangrui Xue
 *
 * \file Oscillator.cpp
 * \brief Implements Oscillator class forcing and damping functions
 * 
 * References:
 *   [Langlois et al. 2016] Toward Animating Water with Complex Acoustic Bubbles
 */

#include <random>

#include "Oscillator.h"


namespace FluidSound {

static const double RHO_WATER = 998.;	// density of water
static const double GAMMA = 1.4;		// gas heat capacity ratio
static const double SIGMA = 0.0726;		// surface tension
static const double MU = 8.9e-4;		// dynamic viscosity of water
static const double GTH = 1.6e6;		// thermal damping constant
static const double G = 1.0;			// TODO: ?
static const double CF = 1497;			// speed of sound in water
static const double ATM = 101325;		// atmospheric pressure


static const double MAX_CUTOFF = 0.0006;

static std::default_random_engine s_forcingRnd;
static std::uniform_real_distribution<double> s_eta(0.4, 1.5);
static std::uniform_real_distribution<double> s_frac(0.4, 0.8);


/** Eq. 14 from [Langlois et al. 2016] */
template <typename T>
std::pair<T, T> Oscillator<T>::CzerskiJetForcing(T radius)
{
    T eta = 0.95;  // TODO

    T cutoff = std::min(MAX_CUTOFF, 0.5 / (3. / radius));   // 1/2 minnaert period

    T pressure_in0 = (ATM + 2. * SIGMA / radius);
    T weight = -9. * GAMMA * SIGMA * eta * pressure_in0 * std::sqrt(1. + eta * eta) / (4. * RHO_WATER * radius * radius * radius);

    T mass = (RHO_WATER / (4. * M_PI * radius));

    return std::pair<T, T>(cutoff, weight / mass);
}

/** Eq. 15 from [Langlois et al. 2016] */
template <typename T>
std::pair<T, T> Oscillator<T>::MergeForcing(T radius, T r1, T r2)
{
    T frac = s_frac(s_forcingRnd);
    T factor = std::pow(2. * SIGMA * r1 * r2 / (RHO_WATER * (r1 + r2)), 0.25);

    T cutoff = std::min(MAX_CUTOFF, 0.5 / (3. / radius));   // 1/2 minnaert period
    T tmp = std::pow(frac * std::min(r1, r2) / 2. / factor, 2);     // TODO: cleanup
    cutoff = std::min(cutoff, tmp);

    T pressure_in0 = (ATM + 2. * SIGMA / radius);
    T weight = 6. * SIGMA * GAMMA * pressure_in0 / (RHO_WATER * radius * radius * radius);

    T mass = (RHO_WATER / (4. * M_PI * radius));

    return std::pair<T, T>(cutoff, weight / mass);
}

/**  */
template <typename T>
T Oscillator<T>::calcBeta(T radius, T w0)
{
    T dr = w0 * radius / CF;
    T dvis = 4 * MU / (RHO_WATER * w0 * radius * radius);
    T phi = 16. * GTH * G / (9 * (GAMMA - 1) * (GAMMA - 1) * w0 / 2. / M_PI);
    T dth = 2 * (std::sqrt(phi - 3) - (3 * GAMMA - 1) / (3 * (GAMMA - 1))) / (phi - 4);

    T dtotal = dr + dvis + dth;

    return w0 * dtotal / std::sqrt(dtotal * dtotal + 4);
}

template class Oscillator<float>;
template struct Oscillator<double>;

}