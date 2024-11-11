/** (c) 2024 Kangrui Xue
 *
 * \file FluidSound.h
 * \brief
 */

#ifndef FLUID_SOUND_H
#define FLUID_SOUND_H

#include <set>

#include "Integrators.h"


namespace FluidSound {

/**
 * \class Solver
 * \brief
 */
template <typename T>
class Solver
{
public:
    /**
     * \brief Constructor: reads Bubble data from file and initializes _oscillators
     * \param[in]  bubFile  path to bubble tracking file
     * \param[in]  dt       timestep size
     * \param[in]  scheme   coupling scheme (0 - uncoupled, 1 - coupled)
     */
    Solver(const std::string& bubFile, double dt, int scheme);

    /** \brief Timesteps Oscillator vibrations */
    T step();

    //void loadState(const std::string &stateFile);
    //void saveState(const std::string &stateFile);

    std::vector<Oscillator<T>>& oscillators() { return _oscillators; }
    const std::vector<double>& eventTimes() { return _eventTimes; }

    /** \brief Prints timings from _integrator */
    void printTimings()
    {
        std::cout << "K,C,F time: " << _integrator->coeff_time.count() << std::endl;
        std::cout << "M^-1 time:  " << _integrator->mass_time.count() << std::endl;
        std::cout << "Solve time: " << _integrator->solve_time.count() << std::endl;
    }

    ~Solver() { delete _integrator; }

private:
    double _dt = 0.;    //!< timestep size
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