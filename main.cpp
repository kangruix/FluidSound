/** (c) 2024 Kangrui Xue
 *
 * @file main.cpp
 */

#include "FluidSound.h"


void Run(const std::string& bubFile, int srate, int scheme)
{
    double dt = 1. / srate;
    FluidSound::Solver solver(bubFile, dt, scheme);
    std::ofstream out_file("output.txt");

    // Simulate from t = 0 to t = 'tf' seconds, logging sum of oscillators in "output.txt"
    int tdx = 0;
    for (double t = 0.; t <= solver.eventTimes().back(); t += dt, ++tdx)
    {
        out_file << solver.step() << std::endl;
        if (tdx % 9600 == 0) std::cout << "At time t = " << t << std::endl;
    }
    solver.printTimings();
}


int main(int argc, char* argv[])
{
    // Default argument values:
    std::string bubFile("trackedBubInfo.txt");  // path to bubble tracking file
    int samplerate = 48000;                     // sampling rate
    int scheme = 1;                             // coupling scheme (0 - uncoupled, 1 - coupled)

    bubFile = "C:/Users/kangruix/Desktop/Stanford/Examples/2016Pour/trackedBubInfo.txt";

    if (argc > 1) bubFile = std::string(argv[1]);
    if (argc > 2) samplerate = std::atoi(argv[2]);
    if (argc > 3) scheme = std::atoi(argv[3]);
        
    Run(bubFile, samplerate, scheme);
}