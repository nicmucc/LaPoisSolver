#include "solvers/gauss_seidel.h"
#include <cmath>

namespace poisson {

double GaussSeidelSolver::iterate(Grid& grid) {
    const size_t nx = grid.nx();
    const size_t ny = grid.ny();
    const double dx2 = grid.dx() * grid.dx();
    const double dy2 = grid.dy() * grid.dy();
    const double denom = 2.0 * (1.0 / dx2 + 1.0 / dy2);

    double max_diff = 0.0;

    for (size_t i = 1; i < nx - 1; ++i) {
        for (size_t j = 1; j < ny - 1; ++j) {
            if (grid.is_boundary(i, j)) continue;

            // Uses already-updated values for (i-1,j) and (i,j-1)
            double phi_new =
                (grid.phi(i + 1, j) + grid.phi(i - 1, j)) / dx2 +
                (grid.phi(i, j + 1) + grid.phi(i, j - 1)) / dy2;

            phi_new = (phi_new + grid.rho(i, j)) / denom;

            double diff = std::abs(phi_new - grid.phi(i, j));
            if (diff > max_diff) max_diff = diff;

            grid.phi(i, j) = phi_new;
        }
    }

    return max_diff;
}

} // namespace poisson
