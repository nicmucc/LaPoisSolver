#pragma once

#include "solvers/solver_base.h"

namespace poisson {

/// Successive Over-Relaxation (SOR) solver.
///
/// Accelerates Gauss-Seidel by extrapolating each update with a
/// relaxation factor ω ∈ (1, 2). The optimal ω for an N×N Laplace
/// problem on a square domain is:
///
///   ω_opt = 2 / (1 + sin(π / N))
///
/// With optimal ω, SOR converges in O(N) iterations vs O(N²) for
/// Gauss-Seidel — a dramatic speedup for large grids.
class SORSolver : public SolverBase {
public:
    /// Construct with explicit relaxation factor
    explicit SORSolver(double omega = 0.0);

    /// Compute the theoretically optimal ω for the given grid
    static double optimal_omega(const Grid& grid);

    double omega() const { return omega_; }
    std::string name() const override;

protected:
    double iterate(Grid& grid) override;

private:
    double omega_;       // Relaxation factor (0 = auto-detect)
    bool   auto_omega_;  // Whether to compute ω from grid dimensions
};

} // namespace poisson
