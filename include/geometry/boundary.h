#pragma once

#include "geometry/grid.h"

namespace poisson {

/// Convenience functions for setting up common boundary configurations
namespace boundary {

/// Parallel plate capacitor: two horizontal plates at y = y_lo and y = y_hi
/// with potentials V_lo and V_hi. Outer boundary grounded.
void parallel_plates(Grid& grid,
                     double y_lo, double y_hi,
                     double V_lo, double V_hi,
                     double plate_x_min, double plate_x_max);

/// Coaxial geometry: inner circle of radius r_inner at V_inner,
/// outer circle of radius r_outer at V_outer.
/// Grid should be centered at (0,0).
void coaxial(Grid& grid,
             double r_inner, double V_inner,
             double r_outer, double V_outer);

/// Electrostatic quadrupole: four electrodes at 45° intervals
/// with alternating ±V, inscribed in a circle of given radius.
/// Produces a field proportional to xy in the center.
void quadrupole(Grid& grid,
                double electrode_radius,
                double aperture_radius,
                double voltage);

/// Point-like charge blob at (cx, cy) with Gaussian profile
/// σ controls the width, Q is the total charge.
void gaussian_charge(Grid& grid,
                     double cx, double cy,
                     double sigma, double Q);

} // namespace boundary
} // namespace poisson
