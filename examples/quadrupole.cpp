#include "geometry/scenarios.h"
#include "solvers/sor.h"
#include "utils/field_export.h"
#include <iostream>
#include <cmath>

using namespace poisson;

int main() {
    std::cout << "=== Electrostatic Quadrupole ===\n\n";
    std::cout << "Geometry: 4 electrodes at 45° intervals, alternating ±5 kV\n";
    std::cout << "Aperture radius: 30 mm, Electrode radius: 8 mm\n\n";

    auto grid = scenarios::quadrupole(128, 5000.0);

    SORSolver solver;
    auto result = solver.solve(grid, 1e-4, 200000);

    std::cout << "Solver: " << solver.name() << "\n";
    std::cout << "Converged: " << (result.converged ? "yes" : "no") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final residual: " << result.final_residual << "\n\n";

    field_export::to_json(grid, "quadrupole.json");
    field_export::convergence_to_json(result.residual_history,
                                      "quadrupole_convergence.json");

    // In an ideal quadrupole, φ ∝ xy near the center
    // Check that the potential has the expected symmetry
    std::cout << "--- Quadrupole Symmetry Check ---\n";
    size_t cx = grid.nx() / 2;
    size_t cy = grid.ny() / 2;
    int offset = 5;

    double phi_pp = grid.phi(cx + offset, cy + offset); // (+x, +y)
    double phi_pn = grid.phi(cx + offset, cy - offset); // (+x, -y)
    double phi_np = grid.phi(cx - offset, cy + offset); // (-x, +y)
    double phi_nn = grid.phi(cx - offset, cy - offset); // (-x, -y)

    std::cout << "  φ(+x,+y) = " << phi_pp << " V\n";
    std::cout << "  φ(+x,-y) = " << phi_pn << " V\n";
    std::cout << "  φ(-x,+y) = " << phi_np << " V\n";
    std::cout << "  φ(-x,-y) = " << phi_nn << " V\n";
    std::cout << "  Expected: φ(+,+) ≈ φ(-,-) and φ(+,-) ≈ φ(-,+) with opposite sign\n";
    std::cout << "  Diagonal symmetry error: "
              << std::abs(phi_pp - phi_nn) << " V\n";
    std::cout << "  Anti-diagonal symmetry error: "
              << std::abs(phi_pn - phi_np) << " V\n";

    // Compute gradient at center to extract quadrupole strength
    auto [Ex, Ey] = grid.compute_efield();
    size_t center_idx = cx * grid.ny() + cy;
    std::cout << "\n  E-field at center: Ex=" << Ex[center_idx]
              << ", Ey=" << Ey[center_idx] << " V/m\n";
    std::cout << "  (Should be ≈ 0 by symmetry)\n";

    std::cout << "\nExported: quadrupole.json\n";
    return 0;
}
