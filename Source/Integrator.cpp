/** (c) 2024 Kangrui Xue
 *
 * @file Integrator.cpp
 */

#include "Integrator.h"

namespace FluidSound {

//##############################################################################
Coupled_Direct::Coupled_Direct(double dt) : Integrator(dt)
{
    RHS.resize(1024);
}

//##############################################################################
void Coupled_Direct::constructMass(double time, double epsSq)
{
    double alpha = (time - _t1) / (_t2 - _t1);

    // Precompute bubble centers and radii
    centers.resize(3, _N_total);
    centers.row(0) = (1. - alpha) * _solveData1.row(2) + (alpha) * _solveData2.row(2);  // x
    centers.row(1) = (1. - alpha) * _solveData1.row(3) + (alpha) * _solveData2.row(3);  // y
    centers.row(2) = (1. - alpha) * _solveData1.row(4) + (alpha) * _solveData2.row(4);  // z
    radii.head(_N_total) = (1. - alpha) * _solveData1.row(0) + alpha * _solveData2.row(0);

    // dense, symmetric mass matrix M
    M.resize(_N_coupled, _N_coupled);
    for (int i = 0; i < _N_coupled; ++i)
    {
        double r_i = radii[i];
        for (int j = i; j < _N_coupled; ++j)
        {
            if (i == j) { M(i, j) = 1.; } // diagonal entries
            else {
                double r_j = radii[j];
                double distSq = (centers.col(j) - centers.col(i)).squaredNorm();
                M(i, j) = 1. / std::sqrt( distSq/(r_i*r_j) + epsSq );
                M(j, i) = M(i, j);
            }
        }
    }
}

//##############################################################################
void Coupled_Direct::refactor()
{
    auto mass_start = std::chrono::steady_clock::now();
    
    constructMass(_t1, epsSq);
    factor1.compute(M.topLeftCorner(_N_coupled, _N_coupled));
    if (factor1.info() == Eigen::NumericalIssue) throw std::runtime_error("Possibly non positive definite matrix!");
    
    constructMass(_t2, epsSq);
    factor2.compute(M.topLeftCorner(_N_coupled, _N_coupled));
    if (factor2.info() == Eigen::NumericalIssue) throw std::runtime_error("Possibly non positive definite matrix!");
    
    auto mass_end = std::chrono::steady_clock::now();
    mass_time += mass_end - mass_start;
}


Eigen::ArrayXd Coupled_Direct::solve(const Eigen::ArrayXd& state, double time)
{
    computeKCF(time);
    
    // solve for y'' | My'' = (F/sqrt(r) - Cy' - Ky)
    auto solve_start = std::chrono::steady_clock::now();
    
    RHS.head(_N_total) = (Fvals.head(_N_total) - Cvals.head(_N_total) * state.segment(_N_total, _N_total) -
        Kvals.head(_N_total) * state.segment(0, _N_total)) / radii.head(_N_total).sqrt();

    double alpha = (time - _t1) / (_t2 - _t1);
    RHS.head(_N_coupled) = (1. - alpha) * factor1.solve(RHS.head(_N_coupled)) + alpha * factor2.solve(RHS.head(_N_coupled));

    _Derivs.resize(2 * _N_total);
    _Derivs.segment(0, _N_total) = state.segment(_N_total, _N_total);
    _Derivs.segment(_N_total, _N_total) = RHS.head(_N_total).array() * radii.head(_N_total).sqrt();
    
    auto solve_end = std::chrono::steady_clock::now();
    solve_time += solve_end - solve_start;
    
    return _Derivs;
}

} // namespace FluidSound