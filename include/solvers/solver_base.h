#pragma once

#include "geometry/grid.h"
#include <string>
#include <vector>

namespace poisson {

/// Convergence result from a solver run
struct SolverResult {
    size_t iterations;           // Number of iterations performed
    double final_residual;       // Final L2 residual norm
    bool   converged;            // Whether tolerance was reached
    std::vector<double> residual_history;  // Residual at each iteration
};

/// Abstract base for iterative Poisson/Laplace solvers.
///
/// All solvers operate on a Grid object, iterating until the L2 norm
/// of the residual drops below a given tolerance or max iterations is hit.
class SolverBase {
public:
    virtual ~SolverBase() = default;

    /// Solve the system. Modifies grid.phi in place.
    SolverResult solve(Grid& grid, double tol = 1e-6, size_t max_iter = 100000);

    /// Compute the L2 residual norm: ||∇²φ + ρ/ε₀||₂
    double compute_residual(const Grid& grid) const;

    /// Human-readable solver name
    virtual std::string name() const = 0;

protected:
    /// Perform one iteration. Returns the max absolute update.
    virtual double iterate(Grid& grid) = 0;
};

} // namespace poisson
