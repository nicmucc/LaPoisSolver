#pragma once

#include "solvers/solver_base.h"

namespace poisson {

/// Jacobi iterative solver.
///
/// The simplest iterative method: updates each interior node using
/// only values from the *previous* iteration. Highly parallelizable
/// but converges slowly.
///
/// Update rule (uniform grid, dx = dy = h):
///   φ_new(i,j) = 0.25 * [φ(i+1,j) + φ(i-1,j) + φ(i,j+1) + φ(i,j-1) + h²·ρ(i,j)]
///
/// For non-uniform grids (dx ≠ dy), the weights adjust accordingly.
class JacobiSolver : public SolverBase {
public:
    std::string name() const override { return "Jacobi"; }

protected:
    double iterate(Grid& grid) override;
};

} // namespace poisson
