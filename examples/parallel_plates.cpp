#include "geometry/scenarios.h"
#include "solvers/sor.h"
#include "utils/field_export.h"
#include "utils/output_path.h"
#include <iostream>
#include <cmath>

using namespace poisson;

int main(int argc, char* argv[]) {
    std::string out = output_dir(argc, argv);
    std::cout << "=== Parallel Plate Capacitor ===\n\n";

    // Create scenario: 128×128 grid, 1 kV across 50mm gap
    auto grid = scenarios::parallel_plates(128, 1000.0, 0.05);

    // Solve with SOR
    SORSolver solver;
    auto result = solver.solve(grid, 1e-6, 500000);

    std::cout << "Solver: " << solver.name() << "\n";
    std::cout << "Converged: " << (result.converged ? "yes" : "no") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final residual: " << result.final_residual << "\n\n";

    // Export for visualization
    field_export::to_json(grid, out + "parallel_plates.json");
    field_export::convergence_to_json(result.residual_history,
                                      out + "parallel_plates_convergence.json");

    std::cout << "Exported: " << out << "parallel_plates.json\n";
    std::cout << "Exported: " << out << "parallel_plates_convergence.json\n";

    // Validate: between the plates, E should be roughly uniform
    // and φ should be linear in y
    std::cout << "\n--- Analytic Validation ---\n";
    double V = 1000.0, gap = 0.05;
    double E_theory = V / gap;
    size_t mid_i = grid.nx() / 2;

    std::cout << "Expected E_y ≈ " << E_theory << " V/m (between plates)\n";

    auto [Ex, Ey] = grid.compute_efield();
    size_t mid_idx = mid_i * grid.ny() + grid.ny() / 2;
    double Ey_center = std::abs(Ey[mid_idx]);
    std::cout << "Computed |E_y| at center: " << Ey_center << " V/m\n";
    std::cout << "Relative error: "
              << std::abs(Ey_center - E_theory) / E_theory * 100.0
              << "%\n";

    return 0;
}
