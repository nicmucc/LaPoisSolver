#pragma once

#include "geometry/grid.h"

namespace poisson {

/// Pre-built accelerator physics scenarios, ready to solve.
namespace scenarios {

/// Parallel plate capacitor (Laplace, analytic solution available)
/// Returns a grid with plates at y = ±gap/2, voltages ±V/2
Grid parallel_plates(size_t n = 128, double voltage = 1000.0, double gap = 0.05);

/// Coaxial cable (Laplace, analytic: φ = V·ln(r/r_out)/ln(r_in/r_out))
Grid coaxial(size_t n = 128, double V_inner = 1000.0,
             double r_inner = 0.01, double r_outer = 0.05);

/// LHC-style electrostatic quadrupole (Laplace)
Grid quadrupole(size_t n = 128, double voltage = 5000.0);

/// Beam pipe with off-center Gaussian charge distribution (Poisson)
Grid beam_with_charge(size_t n = 128, double beam_offset = 0.005,
                      double beam_sigma = 0.002, double beam_charge = 1e-9);

} // namespace scenarios
} // namespace poisson
