/** (c) 2024 Kangrui Xue
 *
 * \file Integrator.cpp
 */

#include "Integrators.h"

namespace FluidSound {

/**  */
template <typename T>
void Integrator<T>::step(double time)
{
    Eigen::ArrayX<T> k1 = solve(_States, time);
    Eigen::ArrayX<T> k2 = solve(_States + _dt / 2. * k1, time + _dt / 2.);
    Eigen::ArrayX<T> k3 = solve(_States + _dt / 2. * k2, time + _dt / 2.);
    Eigen::ArrayX<T> k4 = solve(_States + _dt * k3, time + _dt);

    _Derivs = (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
    _States += _dt * _Derivs;
}

/** */
template <typename T>
void Integrator<T>::updateData(const std::vector<Oscillator<T>*>& coupled_osc, const std::vector<Oscillator<T>*>& uncoupled_osc,
    double time1, double time2)
{ 
    std::vector<Oscillator<T>*> total_osc(coupled_osc.begin(), coupled_osc.end());
    total_osc.insert(total_osc.end(), uncoupled_osc.begin(), uncoupled_osc.end());
    
    //
    _N_coupled = coupled_osc.size(); _N_total = total_osc.size();
    _t1 = time1; _t2 = time2;

    //
    _solveData1.resize(6, _N_total); _solveData2.resize(6, _N_total);
    for (int i = 0; i < _N_total; i++)
    {
        _solveData1.col(i) = total_osc[i]->interp(time1);
        _solveData2.col(i) = total_osc[i]->interp(time2);
    }

    //
    _forceData1.resize(3, _N_coupled); _forceData2.resize(3, _N_coupled);
    for (int i = 0; i < _N_coupled; i++)
    {
        for (int forceIdx = 0; forceIdx < total_osc[i]->forceData.cols(); forceIdx++)
        {
            if (forceIdx == total_osc[i]->forceData.cols() - 1)
            {
                _forceData1.col(i) = total_osc[i]->forceData.col(forceIdx);
                _forceData2.col(i) = _forceData1.col(i);
            }
            else if (time1 < total_osc[i]->forceData(0, forceIdx + 1))
            {
                _forceData1.col(i) = total_osc[i]->forceData.col(forceIdx);
                _forceData2.col(i) = total_osc[i]->forceData.col(forceIdx + 1);
                break;
            }
        }
    }

    //
    _States.resize(2 * _N_total);
    for (int i = 0; i < _N_total; i++)
    {
        _States(i) = total_osc[i]->state(0);            // v - volume displacement
        _States(i + _N_total) = total_osc[i]->state(1); // v'- volume velocity
    }
    _Derivs.resize(2 * _N_total);
}

/** */
template <typename T>
void Integrator<T>::_computeKCF(double time)
{
    auto coeff_start = std::chrono::steady_clock::now();

    double alpha = (time - _t1) / (_t2 - _t1);

    //
    Eigen::ArrayX<T> w0 = (1. - alpha) * _solveData1.row(1) + alpha * _solveData2.row(1);
    _Kvals = w0 * w0;

    //
    _Cvals = (1. - alpha) * _solveData1.row(5) + alpha * _solveData2.row(5);

    //
    _Fvals = Eigen::ArrayX<T>::Zero(_N_total);
    for (int i = 0; i < _N_coupled; i++)
    {
        if (time > _forceData2(0, i))
        {
            T cutoff = _forceData2(1, i);
            T weight = _forceData2(2, i);
            T t = time - _forceData2(0, i);
            _Fvals[i] = (t < cutoff) * weight * t * t;
        }
        else
        {
            T cutoff = _forceData1(1, i);
            T weight = _forceData1(2, i);
            T t = time - _forceData1(0, i);
            _Fvals[i] = (t < cutoff) * weight * t * t;
        }
    }

    auto coeff_end = std::chrono::steady_clock::now();
    coeff_time += coeff_end - coeff_start;
}

template class Integrator<float>;
template class Integrator<double>;


/**  */
template <typename T>
void Coupled_Direct<T>::_constructMass(double time)
{
    double alpha = (time - _t1) / (_t2 - _t1);

    // Precompute bubble centers and radii
    _centers.resize(3, _N_total);
    _centers.row(0) = (1. - alpha) * _solveData1.row(2) + alpha * _solveData2.row(2);  // x
    _centers.row(1) = (1. - alpha) * _solveData1.row(3) + alpha * _solveData2.row(3);  // y
    _centers.row(2) = (1. - alpha) * _solveData1.row(4) + alpha * _solveData2.row(4);  // z

    _radii = (1. - alpha) * _solveData1.row(0) + alpha * _solveData2.row(0);

    // Dense, symmetric mass matrix M
    _M.resize(_N_coupled, _N_coupled);
    for (int i = 0; i < _N_coupled; ++i)
    {
        T r_i = _radii[i];
        for (int j = i; j < _N_coupled; ++j)
        {
            if (i == j) { _M(i, j) = 1.; } // diagonal entries
            else {
                T r_j = _radii[j];
                T distSq = (_centers.col(j) - _centers.col(i)).squaredNorm();
                _M(i, j) = 1. / std::sqrt(distSq / (r_i * r_j) + _epsSq);
                _M(j, i) = _M(i, j);
            }
        }
    }
}

/**  */
template <typename T>
void Coupled_Direct<T>::refactor()
{
    auto mass_start = std::chrono::steady_clock::now();

    _constructMass(_t1);
    _factor1.compute(_M.topLeftCorner(_N_coupled, _N_coupled));
    if (_factor1.info() == Eigen::NumericalIssue) { throw std::runtime_error("Non positive definite matrix!"); }

    _constructMass(_t2);
    _factor2.compute(_M.topLeftCorner(_N_coupled, _N_coupled));
    if (_factor2.info() == Eigen::NumericalIssue) { throw std::runtime_error("Non positive definite matrix!"); }

    auto mass_end = std::chrono::steady_clock::now();
    mass_time += mass_end - mass_start;
}

/** We re-use _Derivs to avoid having to allocate array memory */
template <typename T>
Eigen::ArrayX<T> Coupled_Direct<T>::solve(const Eigen::ArrayX<T>& States, double time)
{
    _computeKCF(time);

    auto solve_start = std::chrono::steady_clock::now();

    double alpha = (time - _t1) / (_t2 - _t1);
    _radii = (1. - alpha) * _solveData1.row(0) + alpha * _solveData2.row(0);

    // solve for y'' | My'' = (F/sqrt(r) - Cy' - Ky)
    _RHS = (_Fvals - _Cvals * States.segment(_N_total, _N_total) - _Kvals * States.segment(0, _N_total)) / _radii.sqrt();
    _RHS.head(_N_coupled) = (1. - alpha) * _factor1.solve(_RHS.head(_N_coupled)) + alpha * _factor2.solve(_RHS.head(_N_coupled));

    _Derivs.segment(0, _N_total) = States.segment(_N_total, _N_total);
    _Derivs.segment(_N_total, _N_total) = _RHS.array() * _radii.sqrt();

    auto solve_end = std::chrono::steady_clock::now();
    solve_time += solve_end - solve_start;

    return _Derivs;
}

template class Coupled_Direct<float>;
template class Coupled_Direct<double>;


/** */
template <typename T>
Eigen::ArrayX<T> Uncoupled<T>::solve(const Eigen::ArrayX<T>& States, double time)
{
    _computeKCF(time);

    auto solve_start = std::chrono::steady_clock::now();

    _Derivs.segment(0, _N_total) = States.segment(_N_total, _N_total);
    _Derivs.segment(_N_total, _N_total) = _Fvals - _Cvals * States.segment(_N_total, _N_total) - _Kvals * States.segment(0, _N_total);

    auto solve_end = std::chrono::steady_clock::now();
    solve_time += solve_end - solve_start;

    return _Derivs;
}

template class Uncoupled<float>;
template class Uncoupled<double>;

} // namespace FluidSound