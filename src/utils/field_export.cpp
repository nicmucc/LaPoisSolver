#include "utils/field_export.h"
#include <fstream>
#include <iomanip>
#include <cmath>

namespace poisson {
namespace field_export {

void to_json(const Grid& grid, const std::string& filename) {
    std::ofstream f(filename);
    f << std::setprecision(8);

    auto [Ex, Ey] = grid.compute_efield();

    // x array
    f << "{\n  \"nx\": " << grid.nx() << ",\n";
    f << "  \"ny\": " << grid.ny() << ",\n";
    f << "  \"x0\": " << grid.x0() << ",\n";
    f << "  \"x1\": " << grid.x1() << ",\n";
    f << "  \"y0\": " << grid.y0() << ",\n";
    f << "  \"y1\": " << grid.y1() << ",\n";

    // Potential as flattened 2D array (row-major: [ix][iy])
    f << "  \"phi\": [";
    for (size_t i = 0; i < grid.nx(); ++i) {
        f << (i > 0 ? ",\n    " : "\n    ") << "[";
        for (size_t j = 0; j < grid.ny(); ++j) {
            if (j > 0) f << ",";
            f << grid.phi(i, j);
        }
        f << "]";
    }
    f << "\n  ],\n";

    // Ex
    f << "  \"Ex\": [";
    for (size_t i = 0; i < grid.nx(); ++i) {
        f << (i > 0 ? ",\n    " : "\n    ") << "[";
        for (size_t j = 0; j < grid.ny(); ++j) {
            if (j > 0) f << ",";
            f << Ex[i * grid.ny() + j];
        }
        f << "]";
    }
    f << "\n  ],\n";

    // Ey
    f << "  \"Ey\": [";
    for (size_t i = 0; i < grid.nx(); ++i) {
        f << (i > 0 ? ",\n    " : "\n    ") << "[";
        for (size_t j = 0; j < grid.ny(); ++j) {
            if (j > 0) f << ",";
            f << Ey[i * grid.ny() + j];
        }
        f << "]";
    }
    f << "\n  ]\n";

    f << "}\n";
}

void to_csv(const Grid& grid, const std::string& filename) {
    std::ofstream f(filename);
    f << std::setprecision(8);
    f << "x,y,phi,Ex,Ey,E_mag\n";

    auto [Ex, Ey] = grid.compute_efield();

    for (size_t i = 0; i < grid.nx(); ++i) {
        for (size_t j = 0; j < grid.ny(); ++j) {
            size_t idx = i * grid.ny() + j;
            double emag = std::sqrt(Ex[idx] * Ex[idx] + Ey[idx] * Ey[idx]);
            f << grid.x(i) << "," << grid.y(j) << ","
              << grid.phi(i, j) << ","
              << Ex[idx] << "," << Ey[idx] << ","
              << emag << "\n";
        }
    }
}

void convergence_to_json(const std::vector<double>& residual_history,
                         const std::string& filename) {
    std::ofstream f(filename);
    f << std::setprecision(12);
    f << "{\n  \"residuals\": [";
    for (size_t i = 0; i < residual_history.size(); ++i) {
        if (i > 0) f << ", ";
        f << residual_history[i];
    }
    f << "]\n}\n";
}

} // namespace field_export
} // namespace poisson
