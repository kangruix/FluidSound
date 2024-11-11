/** (c) 2024 Kangrui Xue
 *
 * \file BubbleUtils.h
 * \brief
 */

#ifndef _FS_BUBBLE_UTILS_H
#define _FS_BUBBLE_UTILS_H

#define _USE_MATH_DEFINES  // needed for M_PI in Visual Studio
#include <math.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <vector>
#include <map>


namespace FluidSound {

/** */
enum EventType { ENTRAIN, MERGE, SPLIT, COLLAPSE };

/** 
 * \struct Bubble
 * \brief Represents a single, physical bubble within a fluid simulation
 * 
 * Initialization handled by loadBubbleFile()
 */
template <typename T>
struct Bubble
{
    int bubID = -1;         //!< \brief global, unique ID
    T radius = 0.;          //!< \brief effective radius

    double startTime = -1.;
    EventType startType;    //!< \brief start event (entrain, merge, or split)
    double endTime = -1.; 
    EventType endType;      //!< \brief end event (merge, split, or collapse)

    /** \brief IDs of parent Bubbles that this Bubble merges or splits from */
    std::vector<int> prevBubIDs;

    /** \brief IDs of child Bubbles that this Bubble merges or splits into */
    std::vector<int> nextBubIDs;
    
    // solve data (NOTE: does not include Bubble start and end times)
    std::vector<double> solveTimes;
    std::vector<T> w0, x, y, z;

    bool hasSolveData() const { return !solveTimes.empty(); }
};

/**
 * \class BubbleUtils
 * \brief 
 */
template <typename T>
class BubbleUtils
{
public:
    /**
     * \brief Loads entire bubble file from disk
     * \param[out] bubMap   map from bubble ID to Bubble object
     * \param[in]  bubFile  path to bubble file
     */
    static void loadBubbleFile(std::map<int, Bubble<T>>& bubMap, const std::string& bubFile);

    /** 
     * \brief Given a list of bubble IDs, returns the ID of the largest bubble in that list
     * \param[in]  bubIDs  vector of bubble IDs to consider
     * \param[in]  bubMap  map from bubble ID to Bubble object
     * \return  ID of largest bubble
     */
    static int largestBubbleID(const std::vector<int>& bubIDs, const std::map<int, Bubble<T>>& bubMap);

private:
    /** \private Parses input file stream data for a single Bubble */
    static void _parseBubble(std::pair<int, Bubble<T>>& bubPair, std::ifstream& in);
};

} // namespace FluidSound

#endif // #ifndef _FS_BUBBLE_UTILS_H