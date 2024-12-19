/** (c) 2024 Kangrui Xue
 *
 * \file Oscillator.h
 * \brief Declares Oscillator struct, including computation of forcing and damping terms
 */

#ifndef _FS_OSCILLATOR_H
#define _FS_OSCILLATOR_H

#include <Eigen/Dense>

#include "BubbleUtils.h"


namespace FluidSound {

/**
 * \class Oscillator
 * \brief Represents a single oscillator, \f$ \ddot{v}(t) + 2\beta \dot{v}(t) + \omega_0^2 v(t) = p(t) / m \f$
 * 
 * Whereas the Bubble struct corresponds to physical bubbles, the Oscillator struct is more of a
 *   mathematical abstraction, meant to interface efficiently with the Integrator class
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
    T accel = 0.;   /**< \brief current volume acceleration : v'' */

    std::vector<double> solveTimes;  //!< times [ t(0) ... t(N) ] corresponding to solveData
    Eigen::Array<T, 6, Eigen::Dynamic> solveData;
    /**< \brief For indices (0, ..., N), solveData is given by:
     * [[ r(0)  ... r(N)  ]
     *  [ ω0(0) ... ω0(N) ]
     *  [ x(0)  ... x(N)  ]
     *  [ y(0)  ... y(N)  ]
     *  [ z(0)  ... z(N)  ]
     *  [ 2β(0) ... 2β(N) ]]
     */

    /** \brief Returns array of linearly interpolated solve data at specified time */
    Eigen::Array<T, 6, 1> interp(double time)
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
    
    bool operator < (const Oscillator<T>& osc) const { return startTime < osc.startTime; }

    Eigen::Array<T, 3, Eigen::Dynamic> forceData;
    /**< \brief For indices (0, ..., F), forceData is given by:
     * [[ forceTime(0) ... forceTime(F) ]
     *  [ cutoff(0)    ... cutoff(F)    ]
     *  [ weight(0)    ... weight(F)    ]]
     *
     * All forcing functions have the form F(t) = (t < cutoff) * weight * t * t
     *  (where t is relative to the force start time : t = time - forceTime) 
     */

    /** 
     * \brief Neck collapse forcing model from Czerksi/Deane [2008; 2010]
     * \param[in]  radius  bubble equilibrium radius
     * \return  (cutoff, weight) pair
     */
    static std::pair<T, T> CzerskiJetForcing(T radius);
    
    /** 
     * \brief Neck expansion forcing model from [Czerski 2011]
     * \param[in]  radius  bubble equilibrium radius
     * \param[in]  r1, r2  radii of parent bubbles
     * \return  (cutoff, weight) pair
     */
    static std::pair<T, T> MergeForcing(T radius, T r1, T r2);
    
    /** \brief Damping via radiative, viscous, and thermal effects */
    static T calcBeta(T radius, T w0);
    
private:
    int _idx = 0;
};

} // namespace FluidSound

#endif // _FS_OSCILLATOR_H