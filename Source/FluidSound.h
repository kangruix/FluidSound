/*
 *  FluidSound.h/cpp
 */

#ifndef FLUID_SOUND_H
#define FLUID_SOUND_H

#include "Integrator.h"
#include <set>

#define MAX_OSC 1024


namespace FluidSound {

class Solver
{
public:
    // Constructor
    Solver(const std::string& bubFile, double dt, int scheme);

    // Timesteps oscillators
    double step();

    //void loadState(const std::string &stateFile);
    //void saveState(const std::string &stateFile);

    std::vector<Oscillator>& oscillators() { return _oscillators; }
    const std::vector<double>& eventTimes() { return _eventTimes; }

    void printTimings();

private:
    double _dt;
    int _step = 0;
    Integrator* _integrator;

    std::vector<Oscillator*> coupled_osc, uncoupled_osc;

    std::vector<Oscillator> _oscillators;
    int _osID = 0;

    std::vector<double> _eventTimes;
    int _evID = 0;

    void _makeOscillators(const std::map<int, Bubble>& bubMap);
};

}

#endif // #ifndef FLUID_SOUND_H