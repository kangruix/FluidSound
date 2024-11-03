/** (c) 2024 Kangrui Xue
 *
 * @file Integrator.cpp
 */

#include "Integrator.h"

namespace FluidSound {


/**  */
void Coupled_Direct::_constructMass(double time)
{
    double alpha = (time - _t1) / (_t2 - _t1);

    // Precompute bubble centers and radii
    _centers.resize(3, _N_total);
    _centers.row(0) = (1. - alpha) * _solveData1.row(2) + alpha * _solveData2.row(2);  // x
    _centers.row(1) = (1. - alpha) * _solveData1.row(3) + alpha * _solveData2.row(3);  // y
    _centers.row(2) = (1. - alpha) * _solveData1.row(4) + alpha * _solveData2.row(4);  // z

    _radii = (1. - alpha) * _solveData1.row(0) + alpha * _solveData2.row(0);

    // dense, symmetric mass matrix M
    _M.resize(_N_coupled, _N_coupled);
    for (int i = 0; i < _N_coupled; ++i)
    {
        REAL r_i = _radii[i];
        for (int j = i; j < _N_coupled; ++j)
        {
            if (i == j) { _M(i, j) = 1.; } // diagonal entries
            else {
                REAL r_j = _radii[j];
                REAL distSq = (_centers.col(j) - _centers.col(i)).squaredNorm();
                _M(i, j) = 1. / std::sqrt(distSq / (r_i * r_j) + _epsSq);
                _M(j, i) = _M(i, j);
            }
        }
    }
}

/**  */
void Coupled_Direct::refactor()
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
Eigen::ArrayX<REAL> Coupled_Direct::solve(const Eigen::ArrayX<REAL>& States, double time)
{
    computeKCF(time);

    auto solve_start = std::chrono::steady_clock::now();

    double alpha = (time - _t1) / (_t2 - _t1);
    _radii = (1. - alpha) * _solveData1.row(0) + alpha * _solveData2.row(0);

    // solve for y'' | My'' = (F/sqrt(r) - Cy' - Ky)
    _RHS.head(_N_total) = (Fvals.head(_N_total) - Cvals.head(_N_total) * States.segment(_N_total, _N_total) -
        Kvals.head(_N_total) * States.segment(0, _N_total)) / _radii.head(_N_total).sqrt();

    _RHS.head(_N_coupled) = (1. - alpha) * _factor1.solve(_RHS.head(_N_coupled))
        + alpha * _factor2.solve(_RHS.head(_N_coupled));

    _Derivs.segment(0, _N_total) = States.segment(_N_total, _N_total);
    _Derivs.segment(_N_total, _N_total) = _RHS.head(_N_total).array() * _radii.head(_N_total).sqrt();

    auto solve_end = std::chrono::steady_clock::now();
    solve_time += solve_end - solve_start;

    return _Derivs;
}

} // namespace FluidSound