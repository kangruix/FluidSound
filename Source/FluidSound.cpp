#include "FluidSound.h"
#include <iostream>

namespace FluidSound {

//##############################################################################
Solver::Solver(const std::string &bubFile, double dt, int scheme)
/*
 *  IN : bubFile: path to bubble tracking file
 *  IN : dt     : timestep size
 *  IN : scheme : coupling scheme (0 - uncoupled, 1 - coupled)
 */
{
    _dt = dt;
    std::map<int, Bubble> allBubbles;

    std::cout << "Reading bubble data from \"" << bubFile << "\"" << std::endl;
    loadBubbleFile(allBubbles, bubFile);
     
    _makeOscillators(allBubbles);
    std::cout << "Total number of oscillators = " << _oscillators.size() << std::endl;
    
    coupled_osc.reserve(MAX_OSC); uncoupled_osc.reserve(MAX_OSC);

    switch (scheme)
    {
        case 1: _integrator = new Coupled_Direct(uncoupled_osc, coupled_osc); break;
        default: _integrator = new Uncoupled(uncoupled_osc, coupled_osc); break;
    }
}

//##############################################################################
void Solver::printTimings()
{
    std::cout << "K,C,F time: " << _integrator->coeff_time.count() << std::endl;
    std::cout << "M^-1 time: "  << _integrator->mass_time.count() << std::endl;
    std::cout << "Solve time: " << _integrator->solve_time.count() << std::endl;
}

//##############################################################################
double Solver::step()
/* Timesteps Oscillator vibrations.
 *    IN : time    : current sim time
 *    OUT: (return): sum of volume acceleration of all Oscillators
 */
{
    double time = _dt * _step + 0.2;
    while (time >= _eventTimes[_evID])
    {
        if (time < _eventTimes[_evID + 1])
        {
            double time1 = _eventTimes[_evID]; double time2 = _eventTimes[_evID + 1];

            // 1. Check if any Oscillators have ended by 'time1'. We then uncouple them and
            //    continue timestepping to avoid discontinuities.
            int coupled_idx = 0;
            for (Oscillator* osc : coupled_osc)
            {
                if (time1 >= osc->m_endTime)
                {
                    uncoupled_osc.push_back(osc);
                    continue;
                }
                coupled_osc[coupled_idx] = osc;
                ++coupled_idx;
            }
            coupled_osc.resize(coupled_idx);

            // 2. Of uncoupled Oscillators, check if any have decayed sufficiently by 'time1'.
            //    These Oscillators are finally removed.
            int uncoupled_idx = 0;
            for (Oscillator* osc : uncoupled_osc)
            {
                if (osc->is_dead())
                {
                    osc->clear();
                    continue;
                }
                uncoupled_osc[uncoupled_idx] = osc;
                ++uncoupled_idx;
            }
            uncoupled_osc.resize(uncoupled_idx);

            // 3. Check if any Oscillators will start between 'time1' and 'time2'. NOTE: '_oscillators'
            //    sorted by increasing 'm_startTime'.
            while (_osID < _oscillators.size() && time >= _oscillators[_osID].m_startTime && 
                   _oscillators[_osID].m_startTime < time2)
            {
                Oscillator* osc = &(_oscillators[_osID]);
                if (time1 < osc->m_endTime)
                {
                    coupled_osc.push_back(osc);
                }
                _osID += 1;
            }

            _integrator->updateData(time1, time2);
            _integrator->refactor();
        }
        _evID += 1;
    }
    _step += 1;

    std::vector<Oscillator*> total_osc(coupled_osc.begin(), coupled_osc.end());
    total_osc.insert(total_osc.end(), uncoupled_osc.begin(), uncoupled_osc.end());
    size_t num_total = total_osc.size();

    if (num_total == 0) { return 0.; }
    if (num_total > MAX_OSC) { throw std::runtime_error("Too many bubbles for coupling!"); }

    // Package all oscillator states into a single vector...
    Eigen::ArrayXd curStates(2 * num_total);
    for (int i = 0; i < num_total; i++)
    {
        curStates(i) = total_osc[i]->m_state(0);           // v - volume displacement
        curStates(i + num_total) = total_osc[i]->m_state(1); // v'- volume velocity
    }
    curStates = _integrator->step(curStates, time, _dt);

    // ...and unpackage vector to update Oscillator states accordingly
    double total_response = 0.0;
    for (int i = 0; i < num_total; i++)
    {
        total_osc[i]->m_state(0) = curStates(i);
        total_osc[i]->m_state(1) = curStates(i + num_total);

        total_osc[i]->m_accel = _integrator->curDeriv(i + num_total);
        total_response += _integrator->curDeriv(i + num_total);
    }
    if (std::abs(total_response) > 100.) { throw std::runtime_error("Instability detected!"); }

    return total_response;
}


//##############################################################################
void Solver::_makeOscillators(const std::map<int, Bubble> &bubMap)
/*
 *  IN : bubMap:
 */
{
    _oscillators.clear();
    std::set<double> eventTimesSet;
    
    std::set<int> used;
    for (const std::pair<int, Bubble>& bubPair : bubMap)
    {
        // Skip if bubble has already been used
        if (used.count(bubPair.first)) continue;
        
        int curBubID = bubPair.first;
        const Bubble * curBub = &bubPair.second;

        Oscillator osc;
        osc.m_startTime = curBub->m_startTime;
        std::vector<double> times, r, wfreqs, x, y, z, pressures, Cvals;

        double prevStartTime = -1.;
        while (true)
        {
            bool lastBub = true;
            used.insert(curBubID);
            
            // Add the solve data for this bubble (if there is any)
            times.insert(times.end(), curBub->m_times.begin(), curBub->m_times.end());
            r.insert(r.end(), curBub->m_times.size(), curBub->m_radius);
            wfreqs.insert(wfreqs.end(), curBub->m_wfreqs.begin(), curBub->m_wfreqs.end());
            x.insert(x.end(), curBub->m_x.begin(), curBub->m_x.end());
            y.insert(y.end(), curBub->m_y.begin(), curBub->m_y.end());
            z.insert(z.end(), curBub->m_z.begin(), curBub->m_z.end());
            pressures.insert(pressures.end(), curBub->m_pressures.begin(), curBub->m_pressures.end());

            osc.m_bubIDs.push_back(curBubID);
            
            // Filter forcing events: only use this one if it's as least 1ms since the last one
            //if (curBub->m_startTime - prevStartTime > 0.001)
            //if (curBub->m_startTime <= curBub->m_endTime)
            {
                std::shared_ptr<ForcingFunction> force = makeForcingFunc(curBubID, *curBub, bubMap);
                osc.m_forcing.push_back(std::make_pair(curBub->m_startTime, force));
                prevStartTime = curBub->m_startTime;
            }

            if (curBub->m_endType == Bubble::MERGE)
            {
                const Bubble& nextBub = bubMap.at(curBub->m_nextBubIDs.at(0));
                int largestParent = largestBubble(nextBub.m_prevBubIDs, bubMap);
                
                if (largestParent == curBubID)
                {
                    lastBub = false;
                    curBubID = curBub->m_nextBubIDs.at(0);
                    curBub = &bubMap.at(curBubID);
                }
            }
            else if (curBub->m_endType == Bubble::SPLIT)
            {
                // Sort children in order of size
                std::multimap<double, int, std::greater<double> > children;
                for (int childID : curBub->m_nextBubIDs)
                {
                    children.insert(std::make_pair(bubMap.at(childID).m_radius, childID));
                }
                // Continue this bubble to the largest child bubble where this is the largest parent
                for (auto &child : children)
                {
                    if (used.count(child.second)) continue;
                    
                    int largestParent = largestBubble(bubMap.at(child.second).m_prevBubIDs, bubMap);
                    if (largestParent == curBubID)
                    {
                        lastBub = false;
                        curBubID = child.second;
                        curBub = &bubMap.at(curBubID);
                        break;
                    }
                }
            }
            if (lastBub) { break; }
        }

        for (int i = 0; i < times.size(); i++)
        {
            Cvals.push_back(2. * calcBeta(r[i], wfreqs[i]));
            pressures[i] *= GAMMA / (4./3. * M_PI * std::pow(r[i], 3));  // TODO: why divide by volume?
        }
        osc.m_times = Eigen::Map<Eigen::VectorXd> (times.data(), times.size());
        osc.m_data.resize(7, times.size());
        
        osc.m_data.row(0) = Eigen::Map<Eigen::VectorXd> (r.data(), r.size());
        osc.m_data.row(1) = Eigen::Map<Eigen::VectorXd> (wfreqs.data(), wfreqs.size());
        osc.m_data.row(2) = Eigen::Map<Eigen::VectorXd> (x.data(), x.size());
        osc.m_data.row(3) = Eigen::Map<Eigen::VectorXd> (y.data(), y.size());
        osc.m_data.row(4) = Eigen::Map<Eigen::VectorXd> (z.data(), z.size());
        osc.m_data.row(5) = Eigen::Map<Eigen::VectorXd> (pressures.data(), pressures.size());
        osc.m_data.row(6) = Eigen::Map<Eigen::VectorXd> (Cvals.data(), Cvals.size());

        // End this oscillator
        osc.m_endTime = curBub->m_endTime;
        
        // Skip if there are not enough frequency solves
        if (times.size() < 1) { continue; }
        // Skip high frequency bubbles
        if (osc.m_data.row(1).maxCoeff() > 18000.0 * 2 * M_PI) { continue; }
        // Skip short blips
        if (osc.m_endTime - osc.m_startTime < 3 * 2 * M_PI / wfreqs[0]) { continue; }

        _oscillators.push_back(osc);
        eventTimesSet.insert(osc.m_startTime);
        eventTimesSet.insert(osc.m_endTime);
        
    } // END for (const std::pair<int, Bubble>& bubPair : bubMap)

    std::sort(_oscillators.begin(), _oscillators.end());
    _eventTimes.assign(eventTimesSet.begin(), eventTimesSet.end());
}

} // namespace FluidSound