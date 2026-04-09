#pragma once

#include "solvers/solver_base.h"

namespace poisson {

/// Gauss-Seidel iterative solver.
///
/// Improves on Jacobi by using *already-updated* values from the
/// current iteration as soon as they're available. This typically
/// converges ~2× faster than Jacobi for Poisson-type problems.
///
/// The update is identical to Jacobi except φ(i-1,j) and φ(i,j-1)
/// come from the current sweep rather than the previous iteration.
class GaussSeidelSolver : public SolverBase {
public:
    std::string name() const override { return "Gauss-Seidel"; }

protected:
    double iterate(Grid& grid) override;
};

} // namespace poisson
