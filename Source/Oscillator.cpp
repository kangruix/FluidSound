/** (c) 2024 Kangrui Xue
 *
 * \file Oscillator.cpp
 */

#include "Oscillator.h"



namespace FluidSound {

static const double fTimeCutoff = 0.0006;  // TODO

template <typename T>
std::pair<T, T> Oscillator<T>::_CzerskiJetForcing(T radius)
{
    T m_r = radius;  // TODO
    T eta = 0.95;  // TODO

    T cutoff = std::min(fTimeCutoff, std::min(5000., 0.5 / (3.0 / radius)));

    T weight = -9 * GAMMA * SIGMA * m_eta * (ATM + 2 * SIGMA / m_r) * sqrt(1 + m_eta * m_eta) / (4 * RHO_WATER * RHO_WATER * m_r * m_r * m_r * m_r * m_r);
    weight *= m_r;

    // Convert to pressure
    double mrp = RHO_WATER * m_r;
    weight *= mrp;

    double mass = (RHO_WATER / (4. * M_PI * m_r));
    weight /= mass;

    return std::pair<T, T>(cutoff, weight);
}

template <typename T>
std::pair<T, T> Oscillator<T>::_MergeForcing(T radius, T r1, T r2)
{
    T m_r = radius;

    double frac = 0.64;
    frac = s_frac(s_forcingRnd);
    double factor = std::pow(2. * SIGMA * r1 * r2 / (RHO_WATER * (r1 + r2)), 0.25);

    T cutoff = std::min(fTimeCutoff, 0.5 / (3.0 / radius)); // 1/2 minnaert period
    cutoff = std::min(cutoff, std::pow(frac * std::min(r1, r2) / 2. / factor, 2));

    T weight = 6 * SIGMA * GAMMA * (ATM + 2 * SIGMA / m_r) *
        (RHO_WATER * RHO_WATER * m_r * m_r * m_r * m_r * m_r);

    // Convert to radius (instead of fractional radius)
    weight *= m_r;

    // Convert to pressure
    double mrp = RHO_WATER * m_r;
    weight *= mrp;

    double mass = (RHO_WATER / (4. * M_PI * m_r));
    weight /= mass;

    return std::pair<T, T>(cutoff, weight);
}

template <typename T>
T Oscillator<T>::calcBeta(T radius, T w0)
{
    double dr = w0 * radius / CF;
    double dvis = 4 * MU / (RHO_WATER * w0 * radius * radius);
    double phi = 16. * GTH * G / (9 * (GAMMA - 1) * (GAMMA - 1) * w0 / 2. / M_PI);
    double dth = 2 * (std::sqrt(phi - 3) - (3 * GAMMA - 1) / (3 * (GAMMA - 1))) / (phi - 4);

    double dtotal = dr + dvis + dth;

    return w0 * dtotal / std::sqrt(dtotal * dtotal + 4);
}

//template class Oscillator<float>;
template struct Oscillator<double>;

}