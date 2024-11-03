#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include "Oscillator.h"
#include <chrono>
#include <iostream>

namespace FluidSound {

class Integrator
{
public:
    Integrator(std::vector<Oscillator*> &uncoupled_osc, std::vector<Oscillator*> &coupled_osc)
        : _uncoupled_osc(uncoupled_osc), _coupled_osc(coupled_osc)
    {        
        Kvals.resize(1024);
        Cvals.resize(1024);
        Fvals.resize(1024);
        radii.resize(1024);

        coeff_time = std::chrono::duration<double>::zero();
        mass_time = std::chrono::duration<double>::zero();
        solve_time = std::chrono::duration<double>::zero();
    }
    
    // Compute stiffness, damping, and forcing terms
    void computeKCF(double time)
    {
        auto coeff_start = std::chrono::steady_clock::now();
        
        double alpha = (time - _time1) / (_time2 - _time1);

        radii.head(n_total) = (1. - alpha) * _Data1.row(0) + (alpha) * _Data2.row(0);
        Kvals.head(n_total) = (1. - alpha) * _Data1.row(1) + (alpha) * _Data2.row(1);  // w0
        Kvals.head(n_total) *= Kvals.head(n_total);  // K = w0*w0
        Cvals.head(n_total) = (1. - alpha) * _Data1.row(6) + (alpha) * _Data2.row(6);

        Eigen::ArrayXd pressures = (1. - alpha) * _Data1.row(5) + (alpha) * _Data2.row(5);
        Fvals.head(n_total) *= 0.;
        for (int i = 0; i < n_coupled; ++i)
        {
            if (time > _ft2[i]) { Fvals[i] = _forcing2[i]->value(time - _ft2[i]); }
            else { Fvals[i] = _forcing1[i]->value(time - _ft1[i]); }

            /*if (time > _ft2[i])
            {
                double cutoff = _forceData2[i].first;
                double weight = _forceData2[i].second;
                double t = time - _ft2[i];
                if (t < 0. || cutoff < 0.) { std::cout << "rip" << std::endl; }

                Fvals[i] = (t < cutoff) * weight * t * t;
            }
            else
            {
                double cutoff = _forceData1[i].first;
                double weight = _forceData1[i].second;
                double t = time - _ft1[i];
                if (t < 0. || cutoff < 0.) { std::cout << "rip" << std::endl; }

                Fvals[i] = (t < cutoff) * weight * t * t;
            }*/

            /*Oscillator* osc = _coupled_osc[i];
            
            // Compute forcing term.  TODO: (Kangrui) consider factoring out / vectorizing
            int forceIdx = 0;
            while (forceIdx < osc->m_forcing.size() && time >= osc->m_forcing[forceIdx].first)
            {
                if (forceIdx + 1 == osc->m_forcing.size() || time < osc->m_forcing[forceIdx + 1].first)
                {
                    Fvals[i] += osc->m_forcing[forceIdx].second->value(time - osc->m_forcing[forceIdx].first);
                    Fvals[i] /= pressures[i];
                    break;
                }
                forceIdx++;
            }*/
            //Fvals[i] /= pressures[i];
        }
        //Fvals.head(n_total) *= Kvals.head(n_total);
        
        auto coeff_end = std::chrono::steady_clock::now();
        coeff_time += coeff_end - coeff_start;
    }

    void updateData(double time1, double time2) 
    { 
        std::vector<Oscillator*> total_osc(_coupled_osc.begin(), _coupled_osc.end());
        total_osc.insert(total_osc.end(), _uncoupled_osc.begin(), _uncoupled_osc.end());
        n_coupled = _coupled_osc.size(); n_total = total_osc.size();

        _time1 = time1; _time2 = time2;
        _Data1.resize(7, n_total); _Data2.resize(7, n_total);
        for (int i = 0; i < n_total; i++)
        {
            _Data1.col(i) = total_osc[i]->interp(time1);
            _Data2.col(i) = total_osc[i]->interp(time2);
        }

        _ft1.resize(n_coupled); _ft2.resize(n_coupled);
        _forcing1.resize(n_coupled); _forcing2.resize(n_coupled);
        _forceData1.resize(n_coupled); _forceData2.resize(n_coupled);
        for (int i = 0; i < n_coupled; i++)
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
        
        curDeriv = (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
        return curState + dt * curDeriv;
    }
    
    std::chrono::duration<double> coeff_time;
    std::chrono::duration<double> mass_time;
    std::chrono::duration<double> solve_time;
    
    Eigen::ArrayXd curDeriv;
protected:
    std::vector<Oscillator*> &_coupled_osc;
    std::vector<Oscillator*> &_uncoupled_osc;
    int n_coupled = 0; int n_total = 0;

    Eigen::MatrixXd _Data1, _Data2;
    double _time1 = -1., _time2 = -1.;

    std::vector<double> _ft1, _ft2;
    std::vector<std::shared_ptr<ForcingFunction>> _forcing1, _forcing2;
    std::vector<std::pair<double, double>> _forceData1, _forceData2;
    
    Eigen::ArrayXd Kvals;
    Eigen::ArrayXd Cvals;
    Eigen::ArrayXd Fvals;
    Eigen::ArrayXd radii;
};


class Coupled_Direct : public Integrator
{
public:
    Coupled_Direct(std::vector<Oscillator*>& uncoupled_osc, std::vector<Oscillator*>& coupled_osc);

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
    Uncoupled(std::vector<Oscillator*> &uncoupled_osc, std::vector<Oscillator*> &coupled_osc)
        : Integrator(uncoupled_osc, coupled_osc)
    { }
    
    void refactor() { }
    Eigen::ArrayXd solve(const Eigen::ArrayXd &state, double time)
    {
        computeKCF(time);
        
        auto solve_start = std::chrono::steady_clock::now();
        
        curDeriv.resize(2 * n_total);
        curDeriv.segment(0, n_total) = state.segment(n_total, n_total);
        curDeriv.segment(n_total, n_total) = Fvals.head(n_total) - Cvals.head(n_total) * state.segment(n_total, n_total)
            - Kvals.head(n_total) * state.segment(0, n_total);
        
        auto solve_end = std::chrono::steady_clock::now();
        solve_time += solve_end - solve_start;
        
        return curDeriv;
    }
    
private:
};

} // namespace FluidSound

#endif // _INTEGRATOR_H