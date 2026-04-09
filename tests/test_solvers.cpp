#include "geometry/grid.h"
#include "geometry/boundary.h"
#include "geometry/scenarios.h"
#include "solvers/jacobi.h"
#include "solvers/gauss_seidel.h"
#include "solvers/sor.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <string>

using namespace poisson;

static int tests_run = 0;
static int tests_passed = 0;

void check(bool condition, const std::string& name) {
    ++tests_run;
    if (condition) {
        ++tests_passed;
        std::cout << "  [PASS] " << name << "\n";
    } else {
        std::cout << "  [FAIL] " << name << "\n";
    }
}

// ── Grid tests ───────────────────────────────────────────────

void test_grid_creation() {
    std::cout << "\n--- Grid Creation ---\n";
    Grid g(10, 20, 0.0, 1.0, 0.0, 2.0);
    check(g.nx() == 10, "nx = 10");
    check(g.ny() == 20, "ny = 20");
    check(std::abs(g.dx() - 1.0 / 9.0) < 1e-12, "dx correct");
    check(std::abs(g.dy() - 2.0 / 19.0) < 1e-12, "dy correct");
}

void test_boundary_conditions() {
    std::cout << "\n--- Boundary Conditions ---\n";
    Grid g(10, 10, 0.0, 1.0, 0.0, 1.0);
    g.set_boundary_rectangle(5.0);

    check(g.is_boundary(0, 0), "Corner (0,0) is boundary");
    check(g.is_boundary(0, 5), "Edge (0,5) is boundary");
    check(!g.is_boundary(5, 5), "Interior (5,5) is not boundary");
    check(std::abs(g.phi(0, 0) - 5.0) < 1e-12, "Boundary value correct");
}

// ── Solver convergence tests ─────────────────────────────────

template<typename SolverT>
void test_solver_converges(const std::string& name) {
    std::cout << "\n--- " << name << " Convergence ---\n";

    // Simple Laplace problem: top=1V, other sides=0V
    Grid g(32, 32, 0.0, 1.0, 0.0, 1.0);
    g.set_boundary_rectangle(0.0);
    for (size_t i = 0; i < g.nx(); ++i) {
        g.set_dirichlet(i, g.ny() - 1, 1.0);
    }

    SolverT solver;
    auto result = solver.solve(g, 1e-6, 100000);

    check(result.converged, name + " converged");
    check(result.final_residual < 1e-6, name + " residual below tol");
    check(result.iterations > 0, name + " took some iterations");

    // Interior values should be between 0 and 1
    bool in_range = true;
    for (size_t i = 1; i < g.nx() - 1; ++i) {
        for (size_t j = 1; j < g.ny() - 1; ++j) {
            if (g.phi(i, j) < -0.01 || g.phi(i, j) > 1.01) {
                in_range = false;
            }
        }
    }
    check(in_range, name + " solution in [0, 1] range");

    // Symmetry: φ(nx/2, j) should be ≈ φ(nx/2 + 1, j) for central columns
    // (symmetric BCs in x)
    size_t ci = g.nx() / 2;
    double sym_err = 0.0;
    for (size_t j = 0; j < g.ny(); ++j) {
        sym_err += std::abs(g.phi(ci, j) - g.phi(ci - 1, j));
    }
    sym_err /= g.ny();
    check(sym_err < 0.01, name + " approximate x-symmetry");

    std::cout << "  Iterations: " << result.iterations
              << ", Residual: " << result.final_residual << "\n";
}

// ── Analytic validation ──────────────────────────────────────

void test_coaxial_analytic() {
    std::cout << "\n--- Coaxial Analytic Validation ---\n";

    double r_in = 0.01, r_out = 0.05, V_in = 100.0;
    auto grid = scenarios::coaxial(128, V_in, r_in, r_out);

    SORSolver solver;
    solver.solve(grid, 1e-8);

    double log_ratio = std::log(r_in / r_out);
    double max_err = 0.0;
    size_t j_mid = grid.ny() / 2;

    for (size_t i = 0; i < grid.nx(); ++i) {
        double x = grid.x(i);
        double r = std::abs(x);
        if (r > r_in * 1.5 && r < r_out * 0.85) {
            double phi_a = V_in * std::log(r / r_out) / log_ratio;
            double err = std::abs(grid.phi(i, j_mid) - phi_a);
            if (err > max_err) max_err = err;
        }
    }

    check(max_err < 5.0, "Coaxial max error < 5V (got " +
          std::to_string(max_err) + "V)");
}

// ── SOR optimal omega ────────────────────────────────────────

void test_sor_faster_than_gs() {
    std::cout << "\n--- SOR vs Gauss-Seidel Speed ---\n";

    auto grid_gs = scenarios::quadrupole(64, 1000.0);
    auto grid_sor = scenarios::quadrupole(64, 1000.0);

    GaussSeidelSolver gs;
    SORSolver sor;

    auto res_gs = gs.solve(grid_gs, 1e-6);
    auto res_sor = sor.solve(grid_sor, 1e-6);

    check(res_sor.iterations < res_gs.iterations,
          "SOR fewer iterations than GS (" +
          std::to_string(res_sor.iterations) + " vs " +
          std::to_string(res_gs.iterations) + ")");
}

// ── E-field computation ──────────────────────────────────────

void test_efield_parallel_plates() {
    std::cout << "\n--- E-field: Parallel Plates ---\n";

    auto grid = scenarios::parallel_plates(64, 1000.0, 0.05);
    SORSolver solver;
    solver.solve(grid, 1e-8);

    auto [Ex, Ey] = grid.compute_efield();
    size_t ci = grid.nx() / 2;
    size_t cj = grid.ny() / 2;
    size_t idx = ci * grid.ny() + cj;

    double E_expected = 1000.0 / 0.05; // V / gap
    double E_computed = std::abs(Ey[idx]);

    double rel_err = std::abs(E_computed - E_expected) / E_expected;
    check(rel_err < 0.15, "E-field within 15% of analytic (" +
          std::to_string(E_computed) + " vs " +
          std::to_string(E_expected) + " V/m)");
}

// ── Main ─────────────────────────────────────────────────────

int main() {
    std::cout << "╔══════════════════════════════════════╗\n";
    std::cout << "║   Poisson-Laplace Solver Tests       ║\n";
    std::cout << "╚══════════════════════════════════════╝\n";

    test_grid_creation();
    test_boundary_conditions();
    test_solver_converges<JacobiSolver>("Jacobi");
    test_solver_converges<GaussSeidelSolver>("Gauss-Seidel");
    test_solver_converges<SORSolver>("SOR");
    test_coaxial_analytic();
    test_sor_faster_than_gs();
    test_efield_parallel_plates();

    std::cout << "\n══════════════════════════════════════\n";
    std::cout << "Results: " << tests_passed << "/" << tests_run << " passed\n";

    return (tests_passed == tests_run) ? 0 : 1;
}
