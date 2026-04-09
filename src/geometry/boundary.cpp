#include "geometry/boundary.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace poisson {
namespace boundary {

void parallel_plates(Grid& grid,
                     double y_lo, double y_hi,
                     double V_lo, double V_hi,
                     double plate_x_min, double plate_x_max)
{
    // Ground the outer boundary
    grid.set_boundary_rectangle(0.0);

    // Set the plate electrodes
    grid.set_dirichlet_region(
        [&](double x, double y) {
            return x >= plate_x_min && x <= plate_x_max &&
                   std::abs(y - y_lo) < grid.dy() * 0.6;
        },
        [V_lo](double, double) { return V_lo; }
    );

    grid.set_dirichlet_region(
        [&](double x, double y) {
            return x >= plate_x_min && x <= plate_x_max &&
                   std::abs(y - y_hi) < grid.dy() * 0.6;
        },
        [V_hi](double, double) { return V_hi; }
    );
}

void coaxial(Grid& grid,
             double r_inner, double V_inner,
             double r_outer, double V_outer)
{
    // Outer boundary: circle (or just the grid edges as approximation)
    grid.set_dirichlet_region(
        [r_outer](double x, double y) {
            return std::sqrt(x * x + y * y) >= r_outer;
        },
        [V_outer](double, double) { return V_outer; }
    );

    // Inner conductor
    grid.set_dirichlet_region(
        [r_inner](double x, double y) {
            return std::sqrt(x * x + y * y) <= r_inner;
        },
        [V_inner](double, double) { return V_inner; }
    );
}

void quadrupole(Grid& grid,
                double electrode_radius,
                double aperture_radius,
                double voltage)
{
    // Ground the outer boundary
    grid.set_boundary_rectangle(0.0);

    // Four electrodes at 45°, 135°, 225°, 315°
    // Alternating +V, -V, +V, -V
    double angles[] = {M_PI / 4.0, 3.0 * M_PI / 4.0,
                       5.0 * M_PI / 4.0, 7.0 * M_PI / 4.0};
    double signs[] = {1.0, -1.0, 1.0, -1.0};

    for (int k = 0; k < 4; ++k) {
        double cx = aperture_radius * std::cos(angles[k]);
        double cy = aperture_radius * std::sin(angles[k]);
        double V = signs[k] * voltage;

        grid.set_dirichlet_region(
            [cx, cy, electrode_radius](double x, double y) {
                double dx = x - cx;
                double dy = y - cy;
                return std::sqrt(dx * dx + dy * dy) <= electrode_radius;
            },
            [V](double, double) { return V; }
        );
    }
}

void gaussian_charge(Grid& grid,
                     double cx, double cy,
                     double sigma, double Q)
{
    double norm = Q / (2.0 * M_PI * sigma * sigma);
    grid.set_charge_density(
        [cx, cy, sigma, norm](double x, double y) {
            double dx = x - cx;
            double dy = y - cy;
            return -norm * std::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
        }
    );
}

} // namespace boundary
} // namespace poisson
