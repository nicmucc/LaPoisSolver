#include "geometry/grid.h"
#include <cmath>
#include <stdexcept>

namespace poisson {

Grid::Grid(size_t nx, size_t ny,
           double x0, double x1,
           double y0, double y1)
    : nx_(nx), ny_(ny)
    , x0_(x0), x1_(x1)
    , y0_(y0), y1_(y1)
    , dx_((x1 - x0) / (nx - 1))
    , dy_((y1 - y0) / (ny - 1))
    , phi_(nx * ny, 0.0)
    , rho_(nx * ny, 0.0)
    , is_boundary_(nx * ny, false)
    , bc_(nx * ny)
{
    if (nx < 3 || ny < 3)
        throw std::invalid_argument("Grid must be at least 3×3");
    if (x1 <= x0 || y1 <= y0)
        throw std::invalid_argument("Invalid domain bounds");
}

void Grid::set_dirichlet(size_t i, size_t j, double value) {
    size_t idx = i * ny_ + j;
    is_boundary_[idx] = true;
    bc_[idx] = {BCType::Dirichlet, value};
    phi_[idx] = value;
}

void Grid::set_boundary_rectangle(double value) {
    for (size_t i = 0; i < nx_; ++i) {
        set_dirichlet(i, 0, value);
        set_dirichlet(i, ny_ - 1, value);
    }
    for (size_t j = 0; j < ny_; ++j) {
        set_dirichlet(0, j, value);
        set_dirichlet(nx_ - 1, j, value);
    }
}

void Grid::set_charge_density(std::function<double(double, double)> rho_func) {
    for (size_t i = 0; i < nx_; ++i) {
        for (size_t j = 0; j < ny_; ++j) {
            rho_[i * ny_ + j] = rho_func(x(i), y(j));
        }
    }
}

void Grid::set_dirichlet_region(
    std::function<bool(double, double)> mask,
    std::function<double(double, double)> value_func)
{
    for (size_t i = 0; i < nx_; ++i) {
        for (size_t j = 0; j < ny_; ++j) {
            double xi = x(i), yj = y(j);
            if (mask(xi, yj)) {
                set_dirichlet(i, j, value_func(xi, yj));
            }
        }
    }
}

void Grid::reset_potential() {
    for (size_t i = 0; i < nx_; ++i) {
        for (size_t j = 0; j < ny_; ++j) {
            size_t idx = i * ny_ + j;
            if (!is_boundary_[idx]) {
                phi_[idx] = 0.0;
            }
        }
    }
}

std::pair<std::vector<double>, std::vector<double>>
Grid::compute_efield() const {
    std::vector<double> Ex(nx_ * ny_, 0.0);
    std::vector<double> Ey(nx_ * ny_, 0.0);

    for (size_t i = 0; i < nx_; ++i) {
        for (size_t j = 0; j < ny_; ++j) {
            size_t idx = i * ny_ + j;

            // E = -∇φ, using central differences (forward/backward at edges)
            if (i > 0 && i < nx_ - 1) {
                Ex[idx] = -(phi((i + 1), j) - phi((i - 1), j)) / (2.0 * dx_);
            } else if (i == 0) {
                Ex[idx] = -(phi(1, j) - phi(0, j)) / dx_;
            } else {
                Ex[idx] = -(phi(nx_ - 1, j) - phi(nx_ - 2, j)) / dx_;
            }

            if (j > 0 && j < ny_ - 1) {
                Ey[idx] = -(phi(i, j + 1) - phi(i, j - 1)) / (2.0 * dy_);
            } else if (j == 0) {
                Ey[idx] = -(phi(i, 1) - phi(i, 0)) / dy_;
            } else {
                Ey[idx] = -(phi(i, ny_ - 1) - phi(i, ny_ - 2)) / dy_;
            }
        }
    }
    return {Ex, Ey};
}

std::vector<double> Grid::compute_efield_magnitude() const {
    auto [Ex, Ey] = compute_efield();
    std::vector<double> mag(nx_ * ny_);
    for (size_t k = 0; k < nx_ * ny_; ++k) {
        mag[k] = std::sqrt(Ex[k] * Ex[k] + Ey[k] * Ey[k]);
    }
    return mag;
}

} // namespace poisson
