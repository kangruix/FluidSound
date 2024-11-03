/** (c) 2024 Kangrui Xue
 *
 * @file BubbleUtils.cpp
 */

#include "BubbleUtils.h"


namespace FluidSound {

/** */
template <typename T>
void BubbleUtils<T>::_parseBubble(std::pair<int, Bubble<T>> &bubPair, std::ifstream& in)
{
    Bubble<T>& bub = bubPair.second;
    std::string line; std::istringstream is;

    int bubID; double time;
    T freqHz, x, y, z, pressure;

    // Line 1: 'Bub <unique bubble ID> <radius>'
    std::getline(in, line);
    if (line.empty())
    {
        bubPair.first = -1;
        return;
    }
    is = std::istringstream(line.substr(4));
    is >> bub.bubID; bubPair.first = bub.bubID;
    is >> bub.radius;

    // Line 2: '  Start: <event type> <start time> <previous bubble ID(s)>'
    std::getline(in, line);
    switch (line[9])
    {
        case 'N': bub.startType = EventType::ENTRAIN; break;
        case 'M': bub.startType = EventType::MERGE; break;
        case 'S': bub.startType = EventType::SPLIT; break;
        default: throw std::runtime_error("Invalid bubble start event type");
    }
    is = std::istringstream(line.substr(11));
    is >> bub.startTime;
    while (is >> bubID)
    {
        bub.prevBubIDs.push_back(bubID);
    }

    // Line 3-n: '  <time> <freqHz> <x> <y> <z> <pressure inside bubble>'
    std::getline(in, line);
    while (line[2] != 'E' && line[2] != 'B' && !in.eof())
    {
        is = std::istringstream(line);
        is >> time >> freqHz >> x >> y >> z >> pressure;
        std::getline(in, line);

        bub.solveTimes.push_back(time); bub.wfreqs.push_back(2 * M_PI * freqHz);
        bub.x.push_back(x); bub.y.push_back(y); bub.z.push_back(z);
        bub.pressures.push_back(pressure);
    };

    // Line n: '  End: <event type> <end time> <next bubble ID(s)>'...
    if (line[2] == 'E')
    {
        switch (line[7])
        {
            case 'M': bub.endType = EventType::MERGE; break;
            case 'S': bub.endType = EventType::SPLIT; break;
            case 'C': bub.endType = EventType::COLLAPSE; break;
            default: throw std::runtime_error("Invalid bubble end event type");
        }
        is = std::istringstream(line.substr(9));
        is >> bub.endTime;
        while (is >> bubID)
        {
            bub.nextBubIDs.push_back(bubID);
        }

        if (bub.endType == EventType::MERGE && bub.nextBubIDs.size() != 1)
        {
            throw std::runtime_error(std::to_string(bubID) + std::string(" did not merge to one bubble"));
        }
        else if (bub.endType == EventType::SPLIT && bub.nextBubIDs.size() < 2)
        {
            throw std::runtime_error(std::to_string(bubID) + std::string(" split to less than 2 bubbles"));
        }
    }
    // ...however, if end info not present, set to default values
    else if (line[2] == 'B') 
    {
        if (bub.solveTimes.empty()) { bub.endTime = bub.startTime + 0.001; }
        else { bub.endTime = bub.solveTimes.back(); }
        bub.endType = EventType::COLLAPSE;
    }
}

/** */
template <typename T>
void BubbleUtils<T>::loadBubbleFile(std::map<int, Bubble<T>>& bubMap, const std::string& bubFile)
{
    bubMap.clear();
    std::ifstream inFile(bubFile);
    while (inFile.good())
    {
        std::pair<int, Bubble<T>> bubPair;
        _parseBubble(bubPair, inFile);
        if (bubPair.first != -1) { bubMap.insert(bubPair); }
    }
}

/** */
template <typename T>
int BubbleUtils<T>::largestBubbleID(const std::vector<int>& bubIDs, const std::map<int, Bubble<T>>& bubMap)
{
    T max_radius = 0.;
    int max_bubID = 0;
    for (int id : bubIDs)
    {
        if (bubMap.at(id).radius > max_radius)
        {
            max_radius = bubMap.at(id).radius;
            max_bubID = id;
        }
    }
    return max_bubID;
}

template class BubbleUtils<float>;
template class BubbleUtils<double>;

} // namespace FluidSound