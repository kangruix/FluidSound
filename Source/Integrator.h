/** (c) 2024 Kangrui Xue
 *
 * @file Integrator.h
 * @brief Defines classes for numerically integrating coupled (or uncoupled) Oscillator systems
 */

#ifndef _FS_INTEGRATOR_H
#define _FS_INTEGRATOR_H


#include <chrono>

#include "Oscillator.h"


namespace FluidSound {

/**
 * @class Integrator
 * @brief Base RK4 integrator for the Oscillator system Mv''(t) + Cv'(t) + Kv(t) = F(t)
 */
class Integrator
{
public:
    Integrator(double dt) : _dt(dt)
    {        
        Kvals.resize(1024);
        Cvals.resize(1024);
        Fvals.resize(1024);
        radii.resize(1024);
    }
    
    // Compute stiffness, damping, and forcing terms
    void computeKCF(double time)
    {
        auto coeff_start = std::chrono::steady_clock::now();
        
        double alpha = (time - _t1) / (_t2 - _t1);

        radii.head(_N_total) = (1. - alpha) * _solveData1.row(0) + (alpha) * _solveData2.row(0);
        Kvals.head(_N_total) = (1. - alpha) * _solveData1.row(1) + (alpha) * _solveData2.row(1);  // w0
        Kvals.head(_N_total) *= Kvals.head(_N_total);  // K = w0*w0
        Cvals.head(_N_total) = (1. - alpha) * _solveData1.row(6) + (alpha) * _solveData2.row(6);

        Eigen::ArrayXd pressures = (1. - alpha) * _solveData1.row(5) + (alpha) * _solveData2.row(5);
        Fvals.head(_N_total) *= 0.;
        for (int i = 0; i < _N_coupled; i++)
        {
            //if (time > _ft2[i]) { Fvals[i] = _forcing2[i]->value(time - _ft2[i]); }
            //else { Fvals[i] = _forcing1[i]->value(time - _ft1[i]); }

            if (time > _ft2[i])
            {
                double cutoff = _forceData2[i].first;
                double weight = _forceData2[i].second;
                double t = time - _ft2[i];

                Fvals[i] = (t < cutoff) * weight * t * t;
            }
            else
            {
                double cutoff = _forceData1[i].first;
                double weight = _forceData1[i].second;
                double t = time - _ft1[i];

                Fvals[i] = (t < cutoff) * weight * t * t;
            }
            //Fvals[i] /= pressures[i];
        }
        //Fvals.head(_N_total) *= Kvals.head(_N_total);
        
        auto coeff_end = std::chrono::steady_clock::now();
        coeff_time += coeff_end - coeff_start;
    }

    void updateData(const std::vector<Oscillator<REAL>*>& coupled_osc, const std::vector<Oscillator<REAL>*>& uncoupled_osc,
        double time1, double time2)
    { 
        std::vector<Oscillator<REAL>*> total_osc(coupled_osc.begin(), coupled_osc.end());
        total_osc.insert(total_osc.end(), uncoupled_osc.begin(), uncoupled_osc.end());
        
        _N_coupled = coupled_osc.size(); _N_total = total_osc.size();
        _t1 = time1; _t2 = time2;

        _solveData1.resize(7, _N_total); _solveData2.resize(7, _N_total);
        for (int i = 0; i < _N_total; i++)
        {
            _solveData1.col(i) = total_osc[i]->interp(time1);
            _solveData2.col(i) = total_osc[i]->interp(time2);
        }

        _ft1.resize(_N_coupled); _ft2.resize(_N_coupled);
        _forcing1.resize(_N_coupled); _forcing2.resize(_N_coupled);
        _forceData1.resize(_N_coupled); _forceData2.resize(_N_coupled);
        for (int i = 0; i < _N_coupled; i++)
        {
            for (int forceIdx = 0; forceIdx < total_osc[i]->m_forcing.size(); forceIdx++)
            {
                if (forceIdx == total_osc[i]->m_forcing.size() - 1)
                {
                    _ft1[i] = total_osc[i]->m_forcing[forceIdx].first;
                    _ft2[i] = _ft1[i];
                    std::pair<double, double> forcePair(total_osc[i]->m_forcing[forceIdx].second->m_cutoff,
                        total_osc[i]->m_forcing[forceIdx].second->m_weight);
                    _forceData1[i] = forcePair;
                    _forceData2[i] = _forceData1[i];
                    _forcing1[i] = total_osc[i]->m_forcing[forceIdx].second;
                    _forcing2[i] = _forcing1[i];
                }
                else if (time1 < total_osc[i]->m_forcing[forceIdx + 1].first)
                {
                    _ft1[i] = total_osc[i]->m_forcing[forceIdx].first;
                    _ft2[i] = total_osc[i]->m_forcing[forceIdx + 1].first;
                    std::pair<double, double> forcePair1(total_osc[i]->m_forcing[forceIdx].second->m_cutoff,
                        total_osc[i]->m_forcing[forceIdx].second->m_weight);
                    std::pair<double, double> forcePair2(total_osc[i]->m_forcing[forceIdx + 1].second->m_cutoff,
                        total_osc[i]->m_forcing[forceIdx + 1].second->m_weight);
                    _forceData1[i] = forcePair1;
                    _forceData2[i] = forcePair2;
                    _forcing1[i] = total_osc[i]->m_forcing[forceIdx].second;
                    _forcing2[i] = total_osc[i]->m_forcing[forceIdx + 1].second;
                    break;
                }
            }
        }
    }
    
    virtual void refactor() = 0;
    virtual Eigen::ArrayXd solve(const Eigen::ArrayXd &state, double time) = 0;
    
    // Runge-Kutta 4
    const Eigen::ArrayXd step(const Eigen::ArrayXd &curState, double time, double dt)
    {
        Eigen::ArrayXd k1 = solve(curState, time);
        Eigen::ArrayXd k2 = solve(curState + dt/2. * k1, time + dt/2.);
        Eigen::ArrayXd k3 = solve(curState + dt/2. * k2, time + dt/2.);
        Eigen::ArrayXd k4 = solve(curState + dt * k3, time + dt);
        
        _Derivs = (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
        return curState + dt * _Derivs;
    }
    
    //const Eigen::ArrayX<REAL>& States() { return _States; }
    const Eigen::ArrayX<REAL>& Derivs() { return _Derivs; }

protected:
    double _dt = 0.;        //!< timestep size
    int _N_coupled = 0;     //!< number of Oscillators to treat as coupled for the batch
    int _N_total = 0;       //!< total number of Oscillators for the batch

    //Eigen::ArrayX<REAL> _States;    //!< packed state vectors [v v'] (over all active Oscillators)
    Eigen::ArrayX<REAL> _Derivs;    //!< packed derivatives [v' v''] (over all active Oscillators)

    // Batch endpoint times
    double _t1 = -1., _t2 = -1.;

    // Solve data (over all active Oscillators) at times '_t1' and '_t2'
    Eigen::Array<REAL, 7, Eigen::Dynamic> _solveData1, _solveData2;


    std::vector<double> _ft1, _ft2;
    std::vector<std::shared_ptr<ForcingFunction>> _forcing1, _forcing2;
    std::vector<std::pair<double, double>> _forceData1, _forceData2;
    
    Eigen::ArrayXd Kvals;
    Eigen::ArrayXd Cvals;
    Eigen::ArrayXd Fvals;
    Eigen::ArrayXd radii;

public:
    std::chrono::duration<double> coeff_time = std::chrono::duration<double>::zero();
    std::chrono::duration<double> mass_time = std::chrono::duration<double>::zero();
    std::chrono::duration<double> solve_time = std::chrono::duration<double>::zero();
};


class Coupled_Direct : public Integrator
{
public:
    Coupled_Direct(double dt);

    void constructMass(double time, double epsSq);
    void refactor();
    Eigen::ArrayXd solve(const Eigen::ArrayXd& state, double time);
    
private:
    const double epsSq = 4.; // regularization term

    // --- host (CPU)
    Eigen::LLT<Eigen::MatrixXd> factor1, factor2;
    Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor> centers;

    Eigen::MatrixXd M;
    Eigen::VectorXd RHS;

    // --- device (GPU)
    Eigen::MatrixXd L1, L2;
    double *d_cx, *d_cy, *d_cz, *d_r;

    double *d_M;
    double *d_RHS;
};


class Uncoupled : public Integrator
{
public:
    Uncoupled(double dt) : Integrator(dt)
    { }
    
    void refactor() { }
    Eigen::ArrayXd solve(const Eigen::ArrayXd &state, double time)
    {
        computeKCF(time);
        
        auto solve_start = std::chrono::steady_clock::now();
        
        _Derivs.resize(2 * _N_total);
        _Derivs.segment(0, _N_total) = state.segment(_N_total, _N_total);
        _Derivs.segment(_N_total, _N_total) = Fvals.head(_N_total) - Cvals.head(_N_total) * state.segment(_N_total, _N_total)
            - Kvals.head(_N_total) * state.segment(0, _N_total);
        
        auto solve_end = std::chrono::steady_clock::now();
        solve_time += solve_end - solve_start;
        
        return _Derivs;
    }
};

} // namespace FluidSound

#endif // _FS_INTEGRATOR_H