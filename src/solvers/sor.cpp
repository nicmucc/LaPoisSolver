#include "solvers/sor.h"
#include <cmath>
#include <sstream>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace poisson {

SORSolver::SORSolver(double omega)
    : omega_(omega)
    , auto_omega_(omega <= 0.0)
{}

double SORSolver::optimal_omega(const Grid& grid) {
    // Spectral radius of the Jacobi iteration matrix for a rectangular grid
    double rho_j = 0.5 * (
        std::cos(M_PI / grid.nx()) +
        std::cos(M_PI / grid.ny())
    );
    return 2.0 / (1.0 + std::sqrt(1.0 - rho_j * rho_j));
}

std::string SORSolver::name() const {
    std::ostringstream ss;
    ss << "SOR (ω=" << omega_ << ")";
    return ss.str();
}

double SORSolver::iterate(Grid& grid) {
    // Auto-detect optimal omega on first call
    if (auto_omega_) {
        omega_ = optimal_omega(grid);
        auto_omega_ = false;
    }

    const size_t nx = grid.nx();
    const size_t ny = grid.ny();
    const double dx2 = grid.dx() * grid.dx();
    const double dy2 = grid.dy() * grid.dy();
    const double denom = 2.0 * (1.0 / dx2 + 1.0 / dy2);

    double max_diff = 0.0;

    for (size_t i = 1; i < nx - 1; ++i) {
        for (size_t j = 1; j < ny - 1; ++j) {
            if (grid.is_boundary(i, j)) continue;

            // Gauss-Seidel update
            double phi_gs =
                (grid.phi(i + 1, j) + grid.phi(i - 1, j)) / dx2 +
                (grid.phi(i, j + 1) + grid.phi(i, j - 1)) / dy2;

            phi_gs = (phi_gs + grid.rho(i, j)) / denom;

            // SOR: blend old and GS values with relaxation factor
            double phi_new = (1.0 - omega_) * grid.phi(i, j) + omega_ * phi_gs;

            double diff = std::abs(phi_new - grid.phi(i, j));
            if (diff > max_diff) max_diff = diff;

            grid.phi(i, j) = phi_new;
        }
    }

    return max_diff;
}

} // namespace poisson
