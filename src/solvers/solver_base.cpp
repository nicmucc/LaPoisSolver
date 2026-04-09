#include "solvers/solver_base.h"
#include <cmath>
#include <iostream>

namespace poisson {

SolverResult SolverBase::solve(Grid& grid, double tol, size_t max_iter) {
    SolverResult result;
    result.converged = false;
    result.residual_history.reserve(max_iter / 10);

    for (size_t iter = 0; iter < max_iter; ++iter) {
        iterate(grid);

        // Check residual every 10 iterations to reduce overhead
        if (iter % 10 == 0 || iter == max_iter - 1) {
            double res = compute_residual(grid);
            result.residual_history.push_back(res);

            if (res < tol) {
                result.iterations = iter + 1;
                result.final_residual = res;
                result.converged = true;
                return result;
            }
        }
    }

    result.iterations = max_iter;
    result.final_residual = compute_residual(grid);
    return result;
}

double SolverBase::compute_residual(const Grid& grid) const {
    double sum_sq = 0.0;
    size_t count = 0;
    double dx2 = grid.dx() * grid.dx();
    double dy2 = grid.dy() * grid.dy();

    for (size_t i = 1; i < grid.nx() - 1; ++i) {
        for (size_t j = 1; j < grid.ny() - 1; ++j) {
            if (grid.is_boundary(i, j)) continue;

            // Residual: ∇²φ + ρ/ε₀ (we absorb ε₀ into ρ for simplicity)
            double laplacian =
                (grid.phi(i + 1, j) - 2.0 * grid.phi(i, j) + grid.phi(i - 1, j)) / dx2 +
                (grid.phi(i, j + 1) - 2.0 * grid.phi(i, j) + grid.phi(i, j - 1)) / dy2;

            double r = laplacian + grid.rho(i, j);
            sum_sq += r * r;
            ++count;
        }
    }

    return (count > 0) ? std::sqrt(sum_sq / count) : 0.0;
}

} // namespace poisson
