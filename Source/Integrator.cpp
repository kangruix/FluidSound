#include "Integrator.h"

namespace FluidSound {

//##############################################################################
Coupled_Direct::Coupled_Direct(std::vector<Oscillator*>& uncoupled_osc, std::vector<Oscillator*>& coupled_osc)
    : Integrator(uncoupled_osc, coupled_osc)
{
    RHS.resize(1024);
}

//##############################################################################
void Coupled_Direct::constructMass(double time, double epsSq)
{
    double alpha = (time - _time1) / (_time2 - _time1);

    // Precompute bubble centers and radii
    centers.resize(3, n_total);
    centers.row(0) = (1. - alpha) * _Data1.row(2) + (alpha) * _Data2.row(2);  // x
    centers.row(1) = (1. - alpha) * _Data1.row(3) + (alpha) * _Data2.row(3);  // y
    centers.row(2) = (1. - alpha) * _Data1.row(4) + (alpha) * _Data2.row(4);  // z
    radii.head(n_total) = (1. - alpha) * _Data1.row(0) + (alpha)*_Data2.row(0);

    // dense, symmetric mass matrix M
    M.resize(n_coupled, n_coupled);
    for (int i = 0; i < n_coupled; ++i)
    {
        double r_i = radii[i];
        for (int j = i; j < n_coupled; ++j)
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
    
    constructMass(_time1, epsSq);
    factor1.compute(M.topLeftCorner(n_coupled, n_coupled));
    if (factor1.info() == Eigen::NumericalIssue) throw std::runtime_error("Possibly non positive definite matrix!");
    
    constructMass(_time2, epsSq);
    factor2.compute(M.topLeftCorner(n_coupled, n_coupled));
    if (factor2.info() == Eigen::NumericalIssue) throw std::runtime_error("Possibly non positive definite matrix!");
    
    auto mass_end = std::chrono::steady_clock::now();
    mass_time += mass_end - mass_start;
}


//##############################################################################
Eigen::ArrayXd Coupled_Direct::solve(const Eigen::ArrayXd& state, double time)
{
    computeKCF(time);
    
    // solve for y'' | My'' = (F/sqrt(r) - Cy' - Ky)
    auto solve_start = std::chrono::steady_clock::now();
    
    RHS.head(n_total) = (Fvals.head(n_total) - Cvals.head(n_total) * state.segment(n_total, n_total) -
        Kvals.head(n_total) * state.segment(0, n_total)) / radii.head(n_total).sqrt();

    double alpha = (time - _time1) / (_time2 - _time1);
    RHS.head(n_coupled) = (1. - alpha) * factor1.solve(RHS.head(n_coupled)) + (alpha)*factor2.solve(RHS.head(n_coupled));

    curDeriv.resize(2 * n_total);
    curDeriv.segment(0, n_total) = state.segment(n_total, n_total);
    curDeriv.segment(n_total, n_total) = RHS.head(n_total).array() * radii.head(n_total).sqrt();
    
    auto solve_end = std::chrono::steady_clock::now();
    solve_time += solve_end - solve_start;
    
    return curDeriv;
}

} // namespace FluidSound