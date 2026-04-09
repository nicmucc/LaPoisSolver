#pragma once

#include <vector>
#include <cstddef>
#include <functional>
#include <string>

namespace poisson {

/// Boundary condition types
enum class BCType {
    Dirichlet,  // Fixed potential value
    Neumann     // Fixed normal derivative (∂φ/∂n)
};

/// Boundary condition specification for a single boundary cell
struct BoundaryCondition {
    BCType type = BCType::Dirichlet;
    double value = 0.0;
};

/// 2D computational grid for finite-difference field solving.
///
/// Stores the scalar potential φ, charge density ρ, and boundary info.
/// Uses row-major storage: index(i, j) = i * ny + j
/// where i is the x-index and j is the y-index.
class Grid {
public:
    /// Construct a grid over [x0, x1] x [y0, y1] with nx × ny nodes
    Grid(size_t nx, size_t ny,
         double x0, double x1,
         double y0, double y1);

    // ── Accessors ────────────────────────────────────────────

    size_t nx() const { return nx_; }
    size_t ny() const { return ny_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double x(size_t i) const { return x0_ + i * dx_; }
    double y(size_t j) const { return y0_ + j * dy_; }

    /// Potential at node (i, j)
    double& phi(size_t i, size_t j) { return phi_[i * ny_ + j]; }
    double  phi(size_t i, size_t j) const { return phi_[i * ny_ + j]; }

    /// Charge density at node (i, j)
    double& rho(size_t i, size_t j) { return rho_[i * ny_ + j]; }
    double  rho(size_t i, size_t j) const { return rho_[i * ny_ + j]; }

    /// Whether node (i, j) is a boundary node
    bool is_boundary(size_t i, size_t j) const { return is_boundary_[i * ny_ + j]; }

    /// Get boundary condition at node (i, j)
    const BoundaryCondition& bc(size_t i, size_t j) const { return bc_[i * ny_ + j]; }

    // ── Setup methods ────────────────────────────────────────

    /// Set Dirichlet BC at a specific node
    void set_dirichlet(size_t i, size_t j, double value);

    /// Set all edge nodes to Dirichlet with given value
    void set_boundary_rectangle(double value = 0.0);

    /// Set charge density using a function ρ(x, y)
    void set_charge_density(std::function<double(double, double)> rho_func);

    /// Set Dirichlet BC using a mask function: if mask(x,y) is true,
    /// set φ = value_func(x,y) at that node
    void set_dirichlet_region(
        std::function<bool(double, double)> mask,
        std::function<double(double, double)> value_func);

    /// Reset potential to zero (keeping BCs intact)
    void reset_potential();

    // ── Field computation ────────────────────────────────────

    /// Compute electric field components Ex, Ey via central differences.
    /// Returns {Ex, Ey} each of size nx*ny.
    std::pair<std::vector<double>, std::vector<double>> compute_efield() const;

    /// Compute field magnitude |E| at each node
    std::vector<double> compute_efield_magnitude() const;

    // ── I/O ──────────────────────────────────────────────────

    /// Raw potential data (row-major)
    const std::vector<double>& phi_data() const { return phi_; }

    /// Domain bounds
    double x0() const { return x0_; }
    double x1() const { return x1_; }
    double y0() const { return y0_; }
    double y1() const { return y1_; }

private:
    size_t nx_, ny_;
    double x0_, x1_, y0_, y1_;
    double dx_, dy_;

    std::vector<double> phi_;             // Scalar potential
    std::vector<double> rho_;             // Charge density
    std::vector<bool>   is_boundary_;     // Boundary flag
    std::vector<BoundaryCondition> bc_;   // Boundary conditions
};

} // namespace poisson
