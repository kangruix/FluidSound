/** (c) 2024 Kangrui Xue
 *
 * @file FluidSound.cpp
 */

#include "FluidSound.h"

namespace FluidSound {

/** */
Solver::Solver(const std::string &bubFile, double dt, int scheme) : _dt(dt)
{
    std::map<int, Bubble<REAL>> allBubbles;

    std::cout << "Reading bubble data from \"" << bubFile << "\"" << std::endl;
    BubbleUtils<REAL>::loadBubbleFile(allBubbles, bubFile);
     
    _makeOscillators(allBubbles);
    std::cout << "Total number of oscillators = " << _oscillators.size() << std::endl;
    
    switch (scheme)
    {
        case 1: _integrator = new Coupled_Direct(_uncoupled_osc, _coupled_osc); break;
        default: _integrator = new Uncoupled(_uncoupled_osc, _coupled_osc); break;
    }
}

/** */
double Solver::step()
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
            for (Oscillator* osc : _coupled_osc)
            {
                if (time1 >= osc->m_endTime)
                {
                    _uncoupled_osc.push_back(osc);
                    continue;
                }
                _coupled_osc[coupled_idx] = osc;
                ++coupled_idx;
            }
            _coupled_osc.resize(coupled_idx);

            // 2. Of uncoupled Oscillators, check if any have decayed sufficiently by 'time1'.
            //    These Oscillators are finally removed.
            int uncoupled_idx = 0;
            for (Oscillator* osc : _uncoupled_osc)
            {
                if (osc->is_dead())
                {
                    osc->clear();
                    continue;
                }
                _uncoupled_osc[uncoupled_idx] = osc;
                ++uncoupled_idx;
            }
            _uncoupled_osc.resize(uncoupled_idx);

            // 3. Check if any Oscillators will start between 'time1' and 'time2'. NOTE: '_oscillators'
            //    sorted by increasing 'm_startTime'.
            while (_osID < _oscillators.size() && time >= _oscillators[_osID].m_startTime && 
                   _oscillators[_osID].m_startTime < time2)
            {
                Oscillator* osc = &(_oscillators[_osID]);
                if (time1 < osc->m_endTime)
                {
                    _coupled_osc.push_back(osc);
                }
                _osID += 1;
            }

            _integrator->updateData(time1, time2);
            _integrator->refactor();
        }
        _evID += 1;
    }
    _step += 1;

    std::vector<Oscillator*> total_osc(_coupled_osc.begin(), _coupled_osc.end());
    total_osc.insert(total_osc.end(), _uncoupled_osc.begin(), _uncoupled_osc.end());
    size_t num_total = total_osc.size();

    if (num_total == 0) { return 0.; }
    if (num_total > 1024) { throw std::runtime_error("Too many bubbles for coupling!"); }

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
void Solver::_makeOscillators(const std::map<int, Bubble<double>> &bubMap)
/*
 *  IN : bubMap:
 */
{
    _oscillators.clear();
    std::set<double> eventTimesSet;
    
    std::set<int> used;
    for (const std::pair<int, Bubble<REAL>>& bubPair : bubMap)
    {
        // Skip if bubble has already been used
        if (used.count(bubPair.first)) continue;
        
        int curBubID = bubPair.first;
        const Bubble<REAL> * curBub = &bubPair.second;

        Oscillator osc;
        osc.m_startTime = curBub->startTime;
        std::vector<double> solveTimes;
        std::vector<REAL> radii, wfreqs, x, y, z, pressures, Cvals;

        double prevStartTime = -1.;
        while (true)
        {
            bool lastBub = true;
            used.insert(curBubID);
            
            // Add the solve data for this bubble (if there is any)
            solveTimes.insert(solveTimes.end(), curBub->solveTimes.begin(), curBub->solveTimes.end());
            radii.insert(radii.end(), curBub->solveTimes.size(), curBub->radius);
            wfreqs.insert(wfreqs.end(), curBub->wfreqs.begin(), curBub->wfreqs.end());
            x.insert(x.end(), curBub->x.begin(), curBub->x.end());
            y.insert(y.end(), curBub->y.begin(), curBub->y.end());
            z.insert(z.end(), curBub->z.begin(), curBub->z.end());
            pressures.insert(pressures.end(), curBub->pressures.begin(), curBub->pressures.end());

            osc.m_bubIDs.push_back(curBubID);
            
            // Filter forcing events: only use this one if it's as least 1ms since the last one
            //if (curBub->m_startTime - prevStartTime > 0.001)
            //if (curBub->m_startTime <= curBub->m_endTime)
            {
                std::shared_ptr<ForcingFunction> force = makeForcingFunc(curBubID, *curBub, bubMap);
                osc.m_forcing.push_back(std::make_pair(curBub->startTime, force));
                prevStartTime = curBub->startTime;
            }

            if (curBub->endType == EventType::MERGE)
            {
                const Bubble<REAL>& nextBub = bubMap.at(curBub->nextBubIDs.at(0));
                int largestParent = BubbleUtils<REAL>::largestBubbleID(nextBub.prevBubIDs, bubMap);
                
                if (largestParent == curBubID)
                {
                    lastBub = false;
                    curBubID = curBub->nextBubIDs.at(0);
                    curBub = &bubMap.at(curBubID);
                }
            }
            else if (curBub->endType == EventType::SPLIT)
            {
                // Sort children in order of size
                std::multimap<double, int, std::greater<double> > children;
                for (int childID : curBub->nextBubIDs)
                {
                    children.insert(std::make_pair(bubMap.at(childID).radius, childID));
                }
                // Continue this bubble to the largest child bubble where this is the largest parent
                for (auto &child : children)
                {
                    if (used.count(child.second)) continue;
                    
                    int largestParent = BubbleUtils<REAL>::largestBubbleID(bubMap.at(child.second).prevBubIDs, bubMap);
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

        for (int i = 0; i < solveTimes.size(); i++)
        {
            Cvals.push_back(2. * calcBeta(radii[i], wfreqs[i]));
            //pressures[i] *= GAMMA / (4./3. * M_PI * std::pow(r[i], 3));  // TODO: why divide by volume?
        }
        osc.m_times = Eigen::Map<Eigen::VectorXd> (solveTimes.data(), solveTimes.size());
        osc.m_data.resize(7, solveTimes.size());
        
        osc.m_data.row(0) = Eigen::Map<Eigen::VectorXd> (radii.data(), radii.size());
        osc.m_data.row(1) = Eigen::Map<Eigen::VectorXd> (wfreqs.data(), wfreqs.size());
        osc.m_data.row(2) = Eigen::Map<Eigen::VectorXd> (x.data(), x.size());
        osc.m_data.row(3) = Eigen::Map<Eigen::VectorXd> (y.data(), y.size());
        osc.m_data.row(4) = Eigen::Map<Eigen::VectorXd> (z.data(), z.size());
        osc.m_data.row(5) = Eigen::Map<Eigen::VectorXd> (pressures.data(), pressures.size());
        osc.m_data.row(6) = Eigen::Map<Eigen::VectorXd> (Cvals.data(), Cvals.size());

        // End this oscillator
        osc.m_endTime = curBub->endTime;
        
        // Skip if there are not enough frequency solves
        if (solveTimes.size() < 1) { continue; }
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