/** (c) 2024 Kangrui Xue
 *
 * @file Oscillator.h
 * @brief
 */

#ifndef _FS_OSCILLATOR_H
#define _FS_OSCILLATOR_H


#include <random>
#include <Eigen/Dense>

#include "BubbleUtils.h"


namespace FluidSound {

typedef double REAL;

static const REAL RHO_WATER = 998.;		// density of water
static const REAL SIGMA = 0.0726;		// surface tension
static const REAL GAMMA = 1.4;			// gas heat capacity ratio
static const REAL MU = 8.9e-4;			// dynamic viscosity of water
static const REAL GTH = 1.6e6;			// thermal damping constant
static const REAL CF = 1497;			// speed of sound in water
static const REAL G = 1.0;				// 
static const REAL ATM = 101325;			// atmospheric pressure
static const REAL ETA = 0.84;			// tan(40 degrees)
//static const REAL ETA_SPLIT = 0.364;	// tan(20 degrees)
//static const REAL ETA_ENTRAIN = 1.732;	// tan(60 degrees)


class ForcingFunction
{
public:
    virtual ~ForcingFunction() { }
    virtual double value(double t) = 0;

    double m_cutoff = 0.;
    double m_weight = 0.;
};


struct Oscillator
{
    double m_startTime = -1.;
    double m_endTime = -1.;
    Eigen::Vector2d m_state = { 0., 0. }; // [v v']
    double m_accel = 0.;

    Eigen::ArrayXd m_times;
    Eigen::MatrixXd m_data;

    Eigen::ArrayXd interp(double time)
    {
        size_t sz = m_times.size();
        if (time >= m_times(sz - 1)) { return m_data.col(sz - 1); }
        else if (time <= m_times(0)) { return m_data.col(0); }
        while (time < m_times(m_idx)) { --m_idx; }
        while (m_idx < sz - 1 && time > m_times(m_idx + 1)) { ++m_idx; }

        double alpha = (time - m_times(m_idx)) / (m_times(m_idx + 1) - m_times(m_idx));
        return (1. - alpha) * m_data.col(m_idx) + alpha * m_data.col(m_idx + 1);
    }
    bool is_dead() const { return m_state.norm() < 1e-10; }
    void clear()
    {
        m_state.setZero(); m_accel = 0.;
        m_times.resize(0); m_data.resize(0, 0);
        m_forcing.clear();
    }
    bool operator < (const Oscillator& osc) const { return m_startTime < osc.m_startTime; }

    std::vector<int> m_bubIDs;
    std::vector< std::pair<double, std::shared_ptr<ForcingFunction>> > m_forcing;

    

private:
    int m_idx = 0;
};




static const double fTimeCutoff = 0.0006;

static std::default_random_engine s_forcingRnd;
static std::uniform_real_distribution<double> s_eta(0.4, 1.5);
static std::uniform_real_distribution<double> s_frac(0.4, 0.8);

//##############################################################################
class CzerskiJetForcing : public ForcingFunction
{
public:
    CzerskiJetForcing(double r,
        double cutoff = std::numeric_limits<double>::infinity(),
        double eta = ETA,
        bool useLaplace = false,
        bool useModulation = false,
        double multiplier = 1.0)
        : m_r(r),
        //m_cutoff (std::min(cutoff, 0.5 / (3.0 / r))), // 1/2 minnaert period
        //m_cutoff (std::min(fTimeCutoff, std::min(cutoff, 0.5 / (3.0 / r)))), // 1/2 minnaert period
        //m_cutoff (std::min(cutoff, 0.0004)), // 1/2 minnaert period
        m_eta(0.95),
        //m_eta(s_eta(s_forcingRnd)),
        m_useLaplace(useLaplace),
        m_laplaceDone(false),
        m_useModulation(useModulation),
        m_multiplier(multiplier)
    {
        m_cutoff = std::min(fTimeCutoff, std::min(5000., 0.5 / (3.0 / r)));

        m_weight = -9 * GAMMA * SIGMA * m_eta * (ATM + 2 * SIGMA / m_r) * sqrt(1 + m_eta * m_eta) / (4 * RHO_WATER * RHO_WATER * m_r * m_r * m_r * m_r * m_r);
        m_weight *= m_r;

        // Convert to pressure
        double mrp = RHO_WATER * m_r;
        m_weight *= mrp;

        double mass = (RHO_WATER / (4. * M_PI * m_r));
        m_weight /= mass;
    }

    double modulation(double t)
    {
        return 1.;
        return 0.5 - 1.0 / M_PI * atan(5000.0 * (t - m_r));
        return 0.5 * erfc(24. * t / m_r - 5.);
    }

    double value(double t)
    {
        if (t > m_cutoff) { return 0.0; }
        return modulation(t) * m_weight * t * t;
    }

private:
    double m_r;
    //double m_cutoff;
    double m_eta;

    double m_multiplier;

    bool m_useLaplace;
    bool m_laplaceDone;
    bool m_useModulation;

    //double m_weight;
};


//##############################################################################
class MergeForcing : public ForcingFunction
{
public:
    // r: new bubble radius
    // r1: parent bubble 1 radius
    // r2: parent bubble 2 radius
    MergeForcing(double r, double r1, double r2,
        double cutoff = std::numeric_limits<double>::infinity())
        : m_r(r)
    {
        // tlim is the time taken for the expanding radius
        // to reach a fixed fraction of the smaller parent
        // radius. Values in the paper ranged from 0.64 - 0.75

        // Equation 5 from Czerski 2011
        double frac = 0.64;
        frac = s_frac(s_forcingRnd);
        double factor = std::pow(2. * SIGMA * r1 * r2 / (RHO_WATER * (r1 + r2)),
            0.25);

        cutoff = std::min(cutoff, 0.5 / (3.0 / r)); // 1/2 minnaert period
        cutoff = std::min(cutoff, fTimeCutoff); // 1/2 minnaert period
        m_cutoff = std::min(cutoff,
            std::pow(frac * std::min(r1, r2) / 2. / factor,
                2));

        m_weight = 6 * SIGMA * GAMMA * (ATM + 2 * SIGMA / m_r) *
            (RHO_WATER * RHO_WATER * m_r * m_r * m_r * m_r * m_r);

        // Convert to radius (instead of fractional radius)
        m_weight *= m_r;

        // Convert to pressure
        double mrp = RHO_WATER * m_r;
        m_weight *= mrp;

        double mass = (RHO_WATER / (4. * M_PI * m_r));
        m_weight /= mass;
    }

    double modulation(double t)
    {
        // From Czerski paper, doesn't go to zero fast enough
        return 1.;
        return 0.5 - 1 / M_PI * std::atan(3 * (t - m_cutoff) / m_cutoff);

        //return 0.5 * std::erfc(4000. * (t - m_tlim));
    }

    double value(double t)
    {
        if (t > m_cutoff) { return 0; }
        return modulation(t) * m_weight * t * t;
    }

private:
    double m_r;
};


//##############################################################################
class ZeroForcing : public ForcingFunction
{
public:
    double value(double t) {
        return 0.;
    }
};


//##############################################################################
static std::shared_ptr<ForcingFunction>
makeForcingFunc(int curBubID, const Bubble<REAL>& curBub,
    const std::map<int, Bubble<REAL>>& bubMap)
{
    /*
     *    IN : curBubID:
     *    IN : curBub  :
     *    IN : bubMap  :
     */
    std::shared_ptr<ForcingFunction> forcing;
    forcing.reset(new ZeroForcing());

    // If this bubble came from another bubble, set the state correctly
    if (curBub.startType == EventType::SPLIT)
    {
        int prevBub = curBub.prevBubIDs.at(0);
        if (bubMap.at(prevBub).radius >= curBub.radius)
        {
            forcing.reset(new CzerskiJetForcing(curBub.radius,
                5000,
                ETA,
                false,
                false));
        }
    } // END if (curBub.m_startType == Bubble::SPLIT))

    else if (curBub.startType == EventType::MERGE)
    {
        if (curBub.prevBubIDs.size() == 2)
        {
            // Make sure all parent bubbles merged completely with this one
            // Sometimes a bubble can split into two, and the smaller fragment
            // can merge with another one. This is probably just resolution and/or
            // bubble tracking issues
            bool allMerge = true;
            for (auto p : curBub.prevBubIDs)
            {
                allMerge = allMerge && bubMap.at(p).endType == EventType::MERGE;
            }

            if (allMerge)
            {
                int p1 = curBub.prevBubIDs.at(0);
                int p2 = curBub.prevBubIDs.at(1);
                double r1 = 0, r2 = 0;

                if (bubMap.count(p1) && bubMap.count(p2))
                {
                    r1 = bubMap.at(p1).radius;
                    r2 = bubMap.at(p2).radius;
                }
                else
                {
                    throw std::runtime_error("Missing parent");
                }

                if (r1 + r2 > curBub.radius)
                {
                    double v1 = 4. / 3. * M_PI * r1 * r1 * r1;
                    double v2 = 4. / 3. * M_PI * r2 * r2 * r2;
                    double vn = 4. / 3. * M_PI * curBub.radius * curBub.radius * curBub.radius;

                    double diff = v1 + v2 - vn;

                    if (diff <= std::max<double>(v1, v2))
                    {
                        if (v1 > v2)
                        {
                            v1 -= diff;
                        }
                        else
                        {
                            v2 -= diff;
                        }

                        r1 = std::pow(3. / 4. / M_PI * v1, 1. / 3.);
                        r2 = std::pow(3. / 4. / M_PI * v2, 1. / 3.);

                        forcing.reset(new MergeForcing(curBub.radius,
                            r1,
                            r2,
                            5000));
                        //curBub.m_endTime - curBub.m_startTime));
                    }
                }
                else
                {
                    forcing.reset(new MergeForcing(curBub.radius,
                        r1,
                        r2,
                        5000));
                    //curBub.m_endTime - curBub.m_startTime));
                }
            }
        }
    } // END else if (curBub.m_startType == Bubble::MERGE)

    else if (curBub.startType == EventType::ENTRAIN)
    {
        forcing.reset(new CzerskiJetForcing(curBub.radius,
            5000,
            //curBub.m_endTime - curBub.m_startTime,
            //ETA_ENTRAIN,
            ETA,
            true,
            false,
            1));
        //forcing.reset(new ZeroForcing());
    }

    return forcing;
}

//##############################################################################
static double calcBeta(double radius, double w0)
{
    double dr = w0 * radius / CF;
    double dvis = 4 * MU / (RHO_WATER * w0 * radius * radius);
    double phi = 16. * GTH * G / (9 * (GAMMA - 1) * (GAMMA - 1) * w0 / 2. / M_PI);
    double dth = 2 * (std::sqrt(phi - 3) - (3 * GAMMA - 1) / (3 * (GAMMA - 1))) / (phi - 4);

    double dtotal = dr + dvis + dth;

    return w0 * dtotal / std::sqrt(dtotal * dtotal + 4);
}

} // namespace FluidSound

#endif // _FS_OSCILLATOR_H