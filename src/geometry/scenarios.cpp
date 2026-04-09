#include "geometry/scenarios.h"
#include "geometry/boundary.h"

namespace poisson {
namespace scenarios {

Grid parallel_plates(size_t n, double voltage, double gap) {
    double half_gap = gap / 2.0;
    double domain = gap * 2.0;  // Some space around the plates
    Grid grid(n, n, -domain / 2, domain / 2, -domain / 2, domain / 2);

    boundary::parallel_plates(grid,
        -half_gap, half_gap,       // plate y-positions
        -voltage / 2, voltage / 2, // plate voltages
        -domain / 4, domain / 4);  // plate x-extent

    return grid;
}

Grid coaxial(size_t n, double V_inner, double r_inner, double r_outer) {
    double domain = r_outer * 1.1;
    Grid grid(n, n, -domain, domain, -domain, domain);

    boundary::coaxial(grid, r_inner, V_inner, r_outer, 0.0);

    return grid;
}

Grid quadrupole(size_t n, double voltage) {
    double domain = 0.06;  // 60mm half-aperture
    Grid grid(n, n, -domain, domain, -domain, domain);

    boundary::quadrupole(grid,
        0.008,   // electrode radius: 8mm
        0.03,    // aperture radius: 30mm
        voltage);

    return grid;
}

Grid beam_with_charge(size_t n, double beam_offset,
                      double beam_sigma, double beam_charge)
{
    double pipe_radius = 0.04;  // 40mm beam pipe
    double domain = pipe_radius * 1.2;
    Grid grid(n, n, -domain, domain, -domain, domain);

    // Grounded beam pipe
    boundary::coaxial(grid, 0.0, 0.0, pipe_radius, 0.0);

    // Off-center beam
    boundary::gaussian_charge(grid, beam_offset, 0.0, beam_sigma, beam_charge);

    return grid;
}

} // namespace scenarios
} // namespace poisson
