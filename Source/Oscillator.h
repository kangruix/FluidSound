/** (c) 2024 Kangrui Xue
 *
 * \file Oscillator.h
 * \brief
 */

#ifndef _FS_OSCILLATOR_H
#define _FS_OSCILLATOR_H


#include <random>
#include <Eigen/Dense>

#include "BubbleUtils.h"


namespace FluidSound {

typedef double REAL;

static const REAL RHO_WATER = 998.;		// density of water
static const REAL SIGMA = 0.0726;		// surface tension
static const REAL GAMMA = 1.4;			// gas heat capacity ratio
static const REAL MU = 8.9e-4;			// dynamic viscosity of water
static const REAL GTH = 1.6e6;			// thermal damping constant
static const REAL CF = 1497;			// speed of sound in water
static const REAL G = 1.0;				// 
static const REAL ATM = 101325;			// atmospheric pressure
static const REAL ETA = 0.84;			// tan(40 degrees)

/**
 * \class Oscillator
 * \brief
 */
template <typename T>
struct Oscillator
{
    /** \brief vector of IDs of Bubbles belonging to this Oscillator, sorted by increasing start time */
    std::vector<int> bubIDs;

    double startTime = -1.;
    double endTime = -1.;

    /** \brief current volume displacement and volume velocity state vector: [v v'] */
    Eigen::Vector2<T> state = { 0., 0. };
    T accel = 0.;   //<! \brief current volume acceleration: v''

    std::vector<double> solveTimes;
    Eigen::Array<T, 7, Eigen::Dynamic> solveData;
    /* For times (0, ..., N), solveData is given by:
     *  [ radius(0) ... radius(N) ]
     *  [ wfreqs(0) ... wfreqs(N) ]
     *  [ x(0)      ... x(N)      ]
     *  [ y(0)      ... y(N)      ]
     *  [ z(0)      ... z(N)      ]
     */

    /** \brief Returns array of linearly interpolated solve data at specified time */
    Eigen::Array<T, 7, 1> interp(double time)
    {
        if (time >= solveTimes.back()) { return solveData.col(solveTimes.size() - 1); }
        else if (time <= solveTimes[0]) { return solveData.col(0); }

        while (time < solveTimes[_idx]) { _idx--; }
        while (_idx < solveTimes.size() - 1 && time > solveTimes[_idx + 1]) { _idx++; }

        double alpha = (time - solveTimes[_idx]) / (solveTimes[_idx + 1] - solveTimes[_idx]);
        return (1. - alpha) * solveData.col(_idx) + alpha * solveData.col(_idx + 1);
    }

    /** \brief Returns true if this Oscillator has decayed sufficiently */
    bool is_dead() const { return state.norm() < 1e-10; }
    
    bool operator < (const Oscillator& osc) const { return startTime < osc.startTime; }


    //std::vector< std::pair<double, std::shared_ptr<ForcingFunction>> > m_forcing;
    std::vector<std::pair<double, std::pair<T, T>>> m_forcing;


    static std::pair<T, T> _CzerskiJetForcing(T radius);
    static std::pair<T, T> _MergeForcing(T radius, T r1, T r2);
    static T calcBeta(T radius, T w0);

    

private:
    int _idx = 0;
};

static const double fTimeCutoff = 0.0006;

static std::default_random_engine s_forcingRnd;
static std::uniform_real_distribution<double> s_eta(0.4, 1.5);
static std::uniform_real_distribution<double> s_frac(0.4, 0.8);

} // namespace FluidSound

#endif // _FS_OSCILLATOR_H