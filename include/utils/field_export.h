#pragma once

#include "geometry/grid.h"
#include <string>

namespace poisson {
namespace field_export {

/// Export potential and field data as JSON (for Plotly visualization)
/// Format: { "x": [...], "y": [...], "phi": [[...]], "Ex": [[...]], "Ey": [[...]] }
void to_json(const Grid& grid, const std::string& filename);

/// Export as CSV: x, y, phi, Ex, Ey, |E|
void to_csv(const Grid& grid, const std::string& filename);

/// Export solver convergence history as JSON
/// Format: { "iterations": [...], "residuals": [...] }
void convergence_to_json(const std::vector<double>& residual_history,
                         const std::string& filename);

} // namespace field_export
} // namespace poisson
