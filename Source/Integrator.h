/** (c) 2024 Kangrui Xue
 *
 * \file Integrator.h
 * \brief Defines classes for numerically integrating coupled (or uncoupled) Oscillator systems
 */

#ifndef _FS_INTEGRATOR_H
#define _FS_INTEGRATOR_H


#include <chrono>

#include "Oscillator.h"


namespace FluidSound {

/**
 * \class Integrator
 * \brief Base RK4 integrator for the Oscillator system Mv''(t) + Cv'(t) + Kv(t) = F(t)
 */
template <typename T>
class Integrator
{
public:
    Integrator(double dt) : _dt(dt) { }
     
    /** \brief Takes an RK4 integration step */
    void step(double time);

    /**
     * \brief Copies all Oscillator data needed for the integration batch; must be called at start of batch.
     * \param[in]  coupled_osc    Oscillators to treat as coupled
     * \param[in]  uncoupled_osc  Oscillators to treat as uncoupled
     * \param[in]  time1          integration batch start time
     * \param[in]  time2          integration batch end time
     */
    void updateData(const std::vector<Oscillator<T>*>& coupled_osc, const std::vector<Oscillator<T>*>& uncoupled_osc,
        double time1, double time2);

    /** \brief Computes and factorizes mass matrix at batch endpoints _t1 and _t2 */
    virtual void refactor() = 0;

    /**
     * \brief Solves for v''(t) = M^-1 ( F(t) - Cv'(t) - Kv(t) )
     * \param[in]  State  packed state vectors [v v']
     * \param[in]  time   solve time t (must satisfy _t1 <= t <= _t2)
     */
    virtual Eigen::ArrayX<T> solve(const Eigen::ArrayX<T>& State, double time) = 0;
    
    const Eigen::ArrayX<T>& States() { return _States; }
    const Eigen::ArrayX<T>& Derivs() { return _Derivs; }

protected:
    double _dt = 0.;        //!< timestep size
    int _N_coupled = 0;     //!< number of Oscillators to treat as coupled for the batch
    int _N_total = 0;       //!< total number of Oscillators for the batch

    // 
    Eigen::ArrayXd _Kvals, _Cvals, _Fvals;
    void computeKCF(double time);

    Eigen::ArrayX<T> _States;    //!< packed state vectors [v v'] (over all active Oscillators)
    Eigen::ArrayX<T> _Derivs;    //!< packed derivatives [v' v''] (over all active Oscillators)

    // Batch endpoint times
    double _t1 = -1., _t2 = -1.;

    // Solve data (over all active Oscillators) at times '_t1' and '_t2'
    Eigen::Array<T, 6, Eigen::Dynamic> _solveData1, _solveData2;

    // Force data
    Eigen::Array<T, 3, Eigen::Dynamic> _forceData1, _forceData2;

public:
    std::chrono::duration<double> coeff_time = std::chrono::duration<double>::zero();
    std::chrono::duration<double> mass_time = std::chrono::duration<double>::zero();
    std::chrono::duration<double> solve_time = std::chrono::duration<double>::zero();
};

/**
 * \class Coupled_Direct
 * \brief
 */
template <typename T>
class Coupled_Direct : public Integrator<T>
{
public:
    Coupled_Direct(double dt) : Integrator<T>(dt) { }

    void refactor();
    Eigen::ArrayX<T> solve(const Eigen::ArrayX<T>& States, double time);
    
private:
    const T _epsSq = 4.;  //!< regularization term

    /** \private */
    void _constructMass(double time);

#ifndef USE_CUDA
    // --- host (CPU) ---
    Eigen::Matrix<T, 3, Eigen::Dynamic, Eigen::RowMajor> _centers;
    Eigen::ArrayX<T> _radii;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _M;
    Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> _factor1, _factor2;
    Eigen::Vector<T, Eigen::Dynamic> _RHS;
#else
    // --- device (GPU) ---
    Eigen::MatrixXd L1, L2;
    T *d_cx, *d_cy, *d_cz, *d_r;

    T *d_M;
    T *d_RHS;
#endif
};

/**
 * \class Uncoupled
 * \brief
 */
template <typename T>
class Uncoupled : public Integrator<T>
{
public:
    Uncoupled(double dt) : Integrator<T>(dt) { }
    
    void refactor() { }     // dummy function call
    Eigen::ArrayX<T> solve(const Eigen::ArrayX<T>& State, double time);
};

} // namespace FluidSound

#endif // _FS_INTEGRATOR_H