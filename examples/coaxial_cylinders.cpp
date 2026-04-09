#include "geometry/scenarios.h"
#include "solvers/sor.h"
#include "utils/field_export.h"
#include <iostream>
#include <cmath>

using namespace poisson;

int main() {
    std::cout << "=== Coaxial Cylinders ===\n\n";

    double r_inner = 0.01;   // 10mm inner conductor
    double r_outer = 0.05;   // 50mm outer conductor
    double V_inner = 1000.0;

    auto grid = scenarios::coaxial(192, V_inner, r_inner, r_outer);

    SORSolver solver;
    auto result = solver.solve(grid, 1e-6, 500000);

    std::cout << "Solver: " << solver.name() << "\n";
    std::cout << "Converged: " << (result.converged ? "yes" : "no") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final residual: " << result.final_residual << "\n\n";

    field_export::to_json(grid, "coaxial.json");
    field_export::convergence_to_json(result.residual_history,
                                      "coaxial_convergence.json");

    // Analytic solution: φ(r) = V_inner * ln(r/r_outer) / ln(r_inner/r_outer)
    std::cout << "--- Analytic Validation (along x-axis, y=0) ---\n";
    double log_ratio = std::log(r_inner / r_outer);

    size_t j_mid = grid.ny() / 2;
    double max_error = 0.0;
    int samples = 0;

    for (size_t i = 0; i < grid.nx(); ++i) {
        double x = grid.x(i);
        double y = grid.y(j_mid);
        double r = std::sqrt(x * x + y * y);

        if (r > r_inner * 1.5 && r < r_outer * 0.9) {
            double phi_analytic = V_inner * std::log(r / r_outer) / log_ratio;
            double phi_numeric = grid.phi(i, j_mid);
            double error = std::abs(phi_numeric - phi_analytic);
            if (error > max_error) max_error = error;
            ++samples;

            if (samples <= 5) {
                std::cout << "  r=" << r * 1000 << "mm: "
                          << "analytic=" << phi_analytic
                          << "  numeric=" << phi_numeric
                          << "  err=" << error << " V\n";
            }
        }
    }
    std::cout << "Max absolute error over " << samples
              << " samples: " << max_error << " V\n";

    std::cout << "\nExported: coaxial.json\n";
    return 0;
}
