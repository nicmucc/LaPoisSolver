#include "geometry/scenarios.h"
#include "solvers/jacobi.h"
#include "solvers/gauss_seidel.h"
#include "solvers/sor.h"
#include "utils/field_export.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <vector>
#include <string>

using namespace poisson;

struct BenchmarkEntry {
    std::string solver_name;
    size_t grid_size;
    size_t iterations;
    double residual;
    double elapsed_ms;
    bool converged;
    std::vector<double> residual_history;
};

template<typename SolverT>
BenchmarkEntry run_benchmark(const std::string& label, size_t n, double tol) {
    auto grid = scenarios::quadrupole(n, 5000.0);
    SolverT solver;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto result = solver.solve(grid, tol, 200000);
    auto t1 = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    return {label, n, result.iterations, result.final_residual,
            ms, result.converged, result.residual_history};
}

int main() {
    std::cout << "=== Solver Benchmark: Quadrupole Problem ===\n\n";

    double tol = 1e-4;
    std::vector<size_t> grid_sizes = {16, 32, 64};
    std::vector<BenchmarkEntry> entries;

    std::cout << std::left << std::setw(16) << "Solver"
              << std::setw(8)  << "N"
              << std::setw(12) << "Iterations"
              << std::setw(14) << "Residual"
              << std::setw(14) << "Time (ms)"
              << std::setw(10) << "Status"
              << "\n";
    std::cout << std::string(74, '-') << "\n";

    for (size_t n : grid_sizes) {
        auto e1 = run_benchmark<JacobiSolver>("Jacobi", n, tol);
        entries.push_back(e1);

        auto e2 = run_benchmark<GaussSeidelSolver>("Gauss-Seidel", n, tol);
        entries.push_back(e2);

        auto e3 = run_benchmark<SORSolver>("SOR", n, tol);
        entries.push_back(e3);

        for (const auto& e : {e1, e2, e3}) {
            std::cout << std::left << std::setw(16) << e.solver_name
                      << std::setw(8)  << e.grid_size
                      << std::setw(12) << e.iterations
                      << std::scientific << std::setprecision(2)
                      << std::setw(14) << e.residual
                      << std::fixed << std::setprecision(1)
                      << std::setw(14) << e.elapsed_ms
                      << std::setw(10) << (e.converged ? "OK" : "FAIL")
                      << "\n";
        }
        std::cout << "\n";
    }

    // Export benchmark results as JSON for Plotly
    std::ofstream f("benchmark.json");
    f << "{\n  \"entries\": [\n";
    for (size_t i = 0; i < entries.size(); ++i) {
        const auto& e = entries[i];
        f << "    {\"solver\": \"" << e.solver_name
          << "\", \"n\": " << e.grid_size
          << ", \"iterations\": " << e.iterations
          << ", \"elapsed_ms\": " << std::fixed << std::setprecision(2) << e.elapsed_ms
          << ", \"converged\": " << (e.converged ? "true" : "false")
          << ", \"residuals\": [";
        for (size_t j = 0; j < e.residual_history.size(); ++j) {
            if (j > 0) f << ",";
            f << std::scientific << std::setprecision(6) << e.residual_history[j];
        }
        f << "]}";
        if (i + 1 < entries.size()) f << ",";
        f << "\n";
    }
    f << "  ]\n}\n";

    std::cout << "Exported: benchmark.json\n";
    return 0;
}
