/** (c) 2024 Kangrui Xue
 *
 * \file FluidSound.h
 * \brief Public interface for FluidSound library; declares Solver class
 * 
 * Based on code from Timothy Langlois and Ryan Aronson. Thanks to Zhehao Li for reviewing!
 */

#ifndef FLUID_SOUND_H
#define FLUID_SOUND_H

#include <set>

#include "Integrators.h"


namespace FluidSound {

/**
 * \class Solver
 * \brief High-level manager for bubble-based water sound synthesis
 */
template <typename T>
class Solver
{
public:
    /**
     * \brief Constructor: reads Bubble data from file and initializes Oscillators
     * \param[in]  bubFile  path to bubble tracking file
     * \param[in]  dt       timestep size
     * \param[in]  scheme   coupling scheme (0 - uncoupled, 1 - coupled)
     * \param[in]  ts       simulation start time (default 0.)
     */
    Solver(const std::string& bubFile, double dt, int scheme, double ts = 0.);

    /** \brief Timesteps Oscillator vibrations */
    T step();

    //void loadState(const std::string &stateFile);
    //void saveState(const std::string &stateFile);

    /** \brief Returns vector of ALL Oscillators, sorted by start time */
    std::vector<Oscillator<T>>& oscillators() { return _oscillators; }
    
    /** \brief Returns vector of sorted event times (i.e., when bubbles are added or removed) */
    const std::vector<double>& eventTimes() { return _eventTimes; }

    /** \brief Prints timings from Integrator */
    void printTimings()
    {
        std::cout << "K,C,F time: " << _integrator->coeff_time.count() << std::endl;
        std::cout << "M^-1 time:  " << _integrator->mass_time.count() << std::endl;
        std::cout << "Solve time: " << _integrator->solve_time.count() << std::endl;
    }

    ~Solver() { delete _integrator; }

private:
    double _dt = 0.;    //!< timestep size
    double _ts = 0.;    //!< simulation start time
    int _step = 0;      //!< current time step

    Integrator<T>* _integrator;

    std::vector<Oscillator<T>*> _coupled_osc;
    std::vector<Oscillator<T>*> _uncoupled_osc;

    std::vector<Oscillator<T>> _oscillators;   //!< vector of ALL Oscillators, sorted by start time
    int _osID = 0;  //!< current _oscillators index

    std::vector<double> _eventTimes;    //!< vector of sorted event times (i.e., when to refactor the mass matrix)
    int _evID = 0;  //!< current _eventTimes index

    /**
     * \private Given bubble data, chains Bubbles together to form Oscillators
     * \param[in]  bubMap  map from Bubble IDs to Bubble objects
     */
    void _makeOscillators(const std::map<int, Bubble<T>>& bubMap);
};

} // namespace FluidSound

#endif // #ifndef FLUID_SOUND_H