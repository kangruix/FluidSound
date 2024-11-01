/*
 *  BubbleIO.h
 */


#ifndef _BUBBLE_IO_H
#define _BUBBLE_IO_H

#include <string>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include "constants.h"


namespace FluidSound {

class Bubble
{
public:
    enum EventType
    {
        ENTRAIN, MERGE, SPLIT, COLLAPSE
    };
    int m_bubID; double m_radius;
    double m_startTime; EventType m_startType;
    double m_endTime; EventType m_endType;

    // IDs of bubbles merged/split FROM
    std::vector<int> m_prevBubIDs;

    // IDs of bubbles merged/split INTO
    std::vector<int> m_nextBubIDs;
    
    // solve data (NOTE: does not include start and end times)
    std::vector<double> m_times, m_wfreqs, m_x, m_y, m_z, m_pressures;
    bool hasData() const { return !m_times.empty(); }
    
    static EventType parseEventType(char type)
    {
        switch (type)
        {
            case 'N': return Bubble::ENTRAIN; break;
            case 'M': return Bubble::MERGE; break;
            case 'S': return Bubble::SPLIT; break;
            case 'C': return Bubble::COLLAPSE; break;
            default: throw std::runtime_error("Invalid bubble event type");
        }
    }
};


//##############################################################################
static std::pair<int, Bubble> parseBubble(std::ifstream &in)
{
/*  Parses data for a single bubble
 *    IN : in     : input file stream reference
 *    RET: bubPair: pair of bubble ID and 'Bubble' object
 */
    std::pair<int, Bubble> bubPair;
    Bubble &bub = bubPair.second;

    std::string line;
    double time, freqHz, x, y, z, pressure;
    int bubID;
    
    // 1. 'Bub <unique bubble ID> <radius>'
    std::getline(in, line);
    if (line.empty())
    {
        bubPair.first = -1;
        return bubPair;
    }
    std::istringstream is1(line.substr(4));
    is1 >> bub.m_bubID; bubPair.first = bub.m_bubID;
    is1 >> bub.m_radius;

    // 2. '  Start: <event type> <start time> <previous bubble IDs>'
    std::getline(in, line);
    char startType = line[9];
    std::istringstream is2(line.substr(11));
    is2 >> bub.m_startTime;
    while (is2 >> bubID)
    {
        bub.m_prevBubIDs.push_back(bubID);
    }
    bub.m_startType = Bubble::parseEventType(startType);

    // 3. '  <time> <freqHz> <x> <y> <z> <pressure inside bubble>'
    std::getline(in, line);
    while (line[2] != 'E' && line[2] != 'B' && !in.eof())
    {
        std::istringstream is3(line);
        is3 >> time >> freqHz >> x >> y >> z >> pressure;
        std::getline(in, line);

        bub.m_times.push_back(time); bub.m_wfreqs.push_back(2 * M_PI * freqHz);
        bub.m_x.push_back(x); bub.m_y.push_back(y); bub.m_z.push_back(z);
        bub.m_pressures.push_back(pressure);
    };

    // 4. '  End: <event type> <end time> <next bubble IDs>'...
    if (line[2] == 'E')
    {
        char endType = line[7];
        std::istringstream is4(line.substr(9));
        is4 >> bub.m_endTime;
        while (is4 >> bubID)
        {
            bub.m_nextBubIDs.push_back(bubID);
        }
        bub.m_endType = Bubble::parseEventType(endType);
        
        if (bub.m_endType == Bubble::MERGE && bub.m_nextBubIDs.size() != 1)
        {
            throw std::runtime_error(std::to_string(bubID) + std::string(" did not merge to one bubble"));
        }
        else if (bub.m_endType == Bubble::SPLIT && bub.m_nextBubIDs.size() < 2)
        {
            throw std::runtime_error(std::to_string(bubID) + std::string(" split to less than 2 bubbles"));
        }
    }
    // ...however, if end info not present, set to default values
    else if (line[2] == 'B') {
        if (bub.m_times.empty()) { bub.m_endTime = bub.m_startTime + 0.001; }
        else { bub.m_endTime = bub.m_times.back(); }
        bub.m_endType = Bubble::COLLAPSE;
    }

    return bubPair;
}


//##############################################################################
static void loadBubbleFile(std::map<int, Bubble> &bubMap, const std::string &bubFile)
{
/*  Load entire bubble file
 *    OUT: bubMap : map from bubble ID to corresponding Bubble object
 *    IN : bubFile: path to bubble file
 */
    bubMap.clear();
    std::ifstream in(bubFile.c_str());
    while (in.good())
    {
        std::pair<int, Bubble> bubPair = parseBubble(in);
        if (bubPair.first >= 0)
        {
            bubMap.insert(bubPair);
        }
    }
}


//##############################################################################
static int largestBubble(const std::vector<int> &bubIDs,
                         const std::map<int, Bubble> &allBubbles)
{
    double maxR = 0;
    int maxBubID = 0;
    for (int i : bubIDs)
    {
        if (allBubbles.at(i).m_radius > maxR)
        {
            maxBubID = i;
            maxR = allBubbles.at(i).m_radius;
        }
    }
    return maxBubID;
}

} // namespace FluidSound

#endif // #ifndef _BUBBLE_IO_H