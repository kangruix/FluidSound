/** (c) 2024 Kangrui Xue
 *
 * \file FluidSound.cpp
 */

#include "FluidSound.h"

namespace FluidSound {

/** */
template <typename T>
Solver<T>::Solver(const std::string &bubFile, double dt, int scheme) : _dt(dt)
{
    std::map<int, Bubble<T>> allBubbles;

    std::cout << "Reading bubble data from \"" << bubFile << "\"" << std::endl;
    BubbleUtils<T>::loadBubbleFile(allBubbles, bubFile);
     
    _makeOscillators(allBubbles);
    std::cout << "Total number of oscillators = " << _oscillators.size() << std::endl;
    
    switch (scheme)
    {
        case 1: _integrator = new Coupled_Direct<T>(dt); break;
        default: _integrator = new Uncoupled<T>(dt); break;
    }
}

/** */
template <typename T>
T Solver<T>::step()
{
    double time = _dt * _step + 0.2;
    while (_evID < _eventTimes.size() && time >= _eventTimes[_evID])
    {
        if (time < _eventTimes[_evID + 1])
        {
            double time1 = _eventTimes[_evID]; double time2 = _eventTimes[_evID + 1];

            // 1. Check if any Oscillators have ended by time1. We then uncouple them and
            //    continue timestepping to avoid discontinuities.
            int coupled_idx = 0;
            for (Oscillator<T>* osc : _coupled_osc)
            {
                if (time1 >= osc->endTime)
                {
                    _uncoupled_osc.push_back(osc);
                    continue;
                }
                _coupled_osc[coupled_idx] = osc;
                coupled_idx++;
            }
            _coupled_osc.resize(coupled_idx);

            // 2. Of uncoupled Oscillators, check if any have decayed sufficiently by time1.
            //    These Oscillators are finally removed.
            int uncoupled_idx = 0;
            for (Oscillator<T>* osc : _uncoupled_osc)
            {
                if (osc->is_dead())
                {
                    continue;
                }
                _uncoupled_osc[uncoupled_idx] = osc;
                uncoupled_idx++;
            }
            _uncoupled_osc.resize(uncoupled_idx);

            // 3. Check if any Oscillators will start between time1 and time2. NOTE: _oscillators
            //    sorted by increasing startTime.
            while (_osID < _oscillators.size() && time >= _oscillators[_osID].startTime && 
                   _oscillators[_osID].startTime < time2)
            {
                Oscillator<T>* osc = &(_oscillators[_osID]);
                if (time1 < osc->endTime)
                {
                    _coupled_osc.push_back(osc);
                }
                _osID++;
            }

            _integrator->updateData(_coupled_osc, _uncoupled_osc, time1, time2);
            _integrator->refactor();
        }
        _evID++;
    }
    _step++;

    std::vector<Oscillator<T>*> total_osc(_coupled_osc.begin(), _coupled_osc.end());
    total_osc.insert(total_osc.end(), _uncoupled_osc.begin(), _uncoupled_osc.end());
    size_t N_total = total_osc.size();

    if (N_total == 0) { return 0.; }
    //if (N_total > 1024) { throw std::runtime_error("Too many bubbles for coupling!"); }


    _integrator->step(time);

    // Unpack _integrator->States() to update Oscillator states accordingly
    T total_response = 0.0;
    for (int i = 0; i < N_total; i++)
    {
        total_osc[i]->state(0) = _integrator->States()(i);
        total_osc[i]->state(1) = _integrator->States()(i + N_total);

        total_osc[i]->accel = _integrator->Derivs()(i + N_total);
        total_response += total_osc[i]->accel;
    }
    if (std::abs(total_response) > 100.) { throw std::runtime_error("Instability detected!"); }

    return total_response;
}


template <typename T>
void Solver<T>::_makeOscillators(const std::map<int, Bubble<T>> &bubMap)
{
    _oscillators.clear();
    std::set<double> eventTimesSet;
    
    std::set<int> usedBubIDs;
    for (const std::pair<int, Bubble<T>>& bubPair : bubMap)
    {
        // Skip if bubble has already been used
        if (usedBubIDs.count(bubPair.first)) continue;
        
        int curBubID = bubPair.first;
        const Bubble<T>* curBub = &bubPair.second;

        // Initialize Oscillator -- we will set its data as we go 
        Oscillator<T> osc;
        osc.startTime = curBub->startTime;

        // Temporary buffers for solve data and force data
        std::vector<double> solveTimes;
        std::vector<T> s_radii, s_w0, s_x, s_y, s_z;
        std::vector<T> Cvals;   // precomputed damping coefficient
        std::vector<T> forceTimes, f_cutoffs, f_weights;

        while (true)
        {
            bool lastBub = true;
            usedBubIDs.insert(curBubID);
            
            // Add the solve data for this bubble (if there is any)
            solveTimes.insert(solveTimes.end(), curBub->solveTimes.begin(), curBub->solveTimes.end());
            s_radii.insert(s_radii.end(), curBub->solveTimes.size(), curBub->radius);
            s_w0.insert(s_w0.end(), curBub->w0.begin(), curBub->w0.end());
            s_x.insert(s_x.end(), curBub->x.begin(), curBub->x.end());
            s_y.insert(s_y.end(), curBub->y.begin(), curBub->y.end());
            s_z.insert(s_z.end(), curBub->z.begin(), curBub->z.end());

            osc.bubIDs.push_back(curBubID);
            
            // ----- Handle bubble start event: forcing logic -----
            std::pair<T, T> force(0., 0.);

            if (curBub->startType == EventType::SPLIT)
            {
                int parentBubID = curBub->prevBubIDs.at(0);
                if (bubMap.at(parentBubID).radius >= curBub->radius)
                {
                    force = Oscillator<T>::CzerskiJetForcing(curBub->radius);
                }
            }
            else if (curBub->startType == EventType::MERGE)
            {
                if (curBub->prevBubIDs.size() == 2)
                {
                    bool allMerge = true;
                    for (int parentBubID : curBub->prevBubIDs)
                    {
                        allMerge = allMerge && (bubMap.at(parentBubID).endType == EventType::MERGE);
                    }

                    if (allMerge)
                    {
                        int p1 = curBub->prevBubIDs.at(0);
                        int p2 = curBub->prevBubIDs.at(1);

                        T r1 = bubMap.at(p1).radius;
                        T r2 = bubMap.at(p2).radius;

                        if (r1 + r2 > curBub->radius)
                        {
                            T v1 = 4. / 3. * M_PI * r1 * r1 * r1;
                            T v2 = 4. / 3. * M_PI * r2 * r2 * r2;
                            T vn = 4. / 3. * M_PI * curBub->radius * curBub->radius * curBub->radius;

                            T diff = v1 + v2 - vn;
                            if (diff <= std::max(v1, v2))
                            {
                                if (v1 > v2) { v1 -= diff; }
                                else { v2 -= diff; }

                                r1 = std::pow(3. / 4. / M_PI * v1, 1. / 3.);
                                r2 = std::pow(3. / 4. / M_PI * v2, 1. / 3.);

                                force = Oscillator<T>::MergeForcing(curBub->radius, r1, r2);
                            }
                        }
                        else
                        {
                            force = Oscillator<T>::MergeForcing(curBub->radius, r1, r2);
                        }
                    }
                }
            }
            else if (curBub->startType == EventType::ENTRAIN)
            {
                force = Oscillator<T>::CzerskiJetForcing(curBub->radius);
            }
            forceTimes.push_back(curBub->startTime);
            f_cutoffs.push_back(force.first);
            f_weights.push_back(force.second);


            // ----- Handle bubble end event: chaining logic -----
            if (curBub->endType == EventType::MERGE)
            {
                const Bubble<T>& nextBub = bubMap.at(curBub->nextBubIDs.at(0));
                int largestParent = BubbleUtils<T>::largestBubbleID(nextBub.prevBubIDs, bubMap);
                
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
                std::multimap<T, int, std::greater<T>> children;
                for (int childID : curBub->nextBubIDs)
                {
                    children.insert(std::make_pair(bubMap.at(childID).radius, childID));
                }
                // Continue this bubble to the largest child bubble where this is the largest parent
                for (auto &child : children)
                {
                    if (usedBubIDs.count(child.second)) continue;
                    
                    int largestParent = BubbleUtils<T>::largestBubbleID(bubMap.at(child.second).prevBubIDs, bubMap);
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
        // End this Oscillator; next, we decide whether or not to keep it
        osc.endTime = curBub->endTime;

        // Filter out if there is not enough solve data
        if (solveTimes.size() < 1) { continue; }
        // Filter out high freqency Oscillators
        if (*std::max_element(s_w0.begin(), s_w0.end()) > 2 * M_PI * 18000.) { continue; }
        // Filter out short blips
        if (osc.endTime - osc.startTime < 3 * 2 * M_PI / s_w0[0]) { continue; }


        // Transfer solve data from temporary buffers to this Oscillator
        for (int i = 0; i < solveTimes.size(); i++)
        {
            Cvals.push_back(2. * Oscillator<T>::calcBeta(s_radii[i], s_w0[i]));
        }
        osc.solveTimes = solveTimes;
        osc.solveData.resize(6, solveTimes.size());

        osc.solveData.row(0) = Eigen::Map<Eigen::VectorX<T>>(s_radii.data(), s_radii.size());
        osc.solveData.row(1) = Eigen::Map<Eigen::VectorX<T>>(s_w0.data(), s_w0.size());
        osc.solveData.row(2) = Eigen::Map<Eigen::VectorX<T>>(s_x.data(), s_x.size());
        osc.solveData.row(3) = Eigen::Map<Eigen::VectorX<T>>(s_y.data(), s_y.size());
        osc.solveData.row(4) = Eigen::Map<Eigen::VectorX<T>>(s_z.data(), s_z.size());
        osc.solveData.row(5) = Eigen::Map<Eigen::VectorX<T>>(Cvals.data(), Cvals.size());

        osc.forceData.resize(3, forceTimes.size());
        osc.forceData.row(0) = Eigen::Map<Eigen::VectorX<T>>(forceTimes.data(), forceTimes.size());
        osc.forceData.row(1) = Eigen::Map<Eigen::VectorX<T>>(f_cutoffs.data(), f_cutoffs.size());
        osc.forceData.row(2) = Eigen::Map<Eigen::VectorX<T>>(f_weights.data(), f_weights.size());


        // Finally, add this Oscillator
        _oscillators.push_back(osc);
        eventTimesSet.insert(osc.startTime);
        eventTimesSet.insert(osc.endTime);
        
    } // END for (const std::pair<int, Bubble>& bubPair : bubMap)

    std::sort(_oscillators.begin(), _oscillators.end());
    _eventTimes.assign(eventTimesSet.begin(), eventTimesSet.end());
}

template class Solver<float>;
template class Solver<double>;

} // namespace FluidSound