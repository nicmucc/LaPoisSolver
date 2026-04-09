# Poisson-Laplace Field Solver

A 2D finite-difference solver for electrostatic and magnetostatic field problems, with a focus on **accelerator physics** geometries inspired by CERN's LHC infrastructure.

Built in **C++17** with interactive **Plotly** visualizations via Python.

---

## Physics Background

This solver computes the scalar potential **φ** for two fundamental PDEs of classical electromagnetism:

| Equation | Form | Use case |
|---|---|---|
| **Laplace** | ∇²φ = 0 | Regions free of charge (electrode problems) |
| **Poisson** | ∇²φ = −ρ/ε₀ | Regions with space charge (beam dynamics) |

Once **φ** is known, the electric field follows from **E** = −∇φ.

The solver discretizes the domain on a uniform 2D grid and applies iterative relaxation methods until the residual drops below a configurable tolerance.

## Solver Methods

Three iterative solvers are implemented, each building on the last:

### Jacobi
The simplest approach: update every interior node from its four neighbours using **only previous-iteration values**. Highly parallelizable but slow to converge.

### Gauss-Seidel
Uses **already-updated values** from the current sweep as soon as they are available. Roughly 2× faster convergence than Jacobi with no extra memory cost.

### Successive Over-Relaxation (SOR)
Accelerates Gauss-Seidel by extrapolating each update with a relaxation factor **ω ∈ (1, 2)**. With the optimal ω (computed from the grid's spectral radius), convergence improves from O(N²) to O(N) iterations — a dramatic speedup for large grids.

```
ω_opt = 2 / (1 + sin(π/N))
```

## Accelerator Physics Scenarios

Pre-built scenarios with physically meaningful parameters:

| Scenario | Type | Description |
|---|---|---|
| **Parallel Plates** | Laplace | Classic capacitor with analytic validation |
| **Coaxial Cylinders** | Laplace | Cylindrical geometry, validates against φ = V·ln(r/r₂)/ln(r₁/r₂) |
| **Electrostatic Quadrupole** | Laplace | Four electrodes at ±5 kV, LHC-scale aperture (30 mm) |
| **Beam with Space Charge** | Poisson | Off-center Gaussian beam in a grounded pipe |

## Project Structure

```
poisson-laplace-solver/
├── CMakeLists.txt              # Build system
├── include/
│   ├── solvers/
│   │   ├── solver_base.h       # Abstract solver interface
│   │   ├── jacobi.h            # Jacobi iteration
│   │   ├── gauss_seidel.h      # Gauss-Seidel iteration
│   │   └── sor.h               # SOR with optimal ω
│   ├── geometry/
│   │   ├── grid.h              # 2D grid + potential storage
│   │   ├── boundary.h          # Boundary condition helpers
│   │   └── scenarios.h         # Pre-built accelerator geometries
│   └── utils/
│       └── field_export.h      # JSON/CSV export for visualization
├── src/                        # Implementation files (.cpp)
├── examples/
│   ├── parallel_plates.cpp     # Capacitor with analytic validation
│   ├── coaxial_cylinders.cpp   # Coaxial cable with analytic check
│   ├── quadrupole.cpp          # Electrostatic quadrupole
│   └── solver_benchmark.cpp    # Head-to-head solver comparison
├── tests/
│   └── test_solvers.cpp        # Unit + validation tests
├── scripts/
│   └── visualize.py            # Plotly interactive visualization
└── docs/                       # Theory notes
```

## Building

### Requirements
- C++17 compiler (GCC ≥ 7, Clang ≥ 5, MSVC ≥ 2017)
- CMake ≥ 3.16
- Python 3.8+ with `plotly` and `numpy` (for visualization)

### Compile

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Run Tests

```bash
cd build
ctest --output-on-failure
# or directly:
./test_solvers
```

### Run Examples

```bash
cd build

# Solve and export field data as JSON
./parallel_plates
./coaxial_cylinders
./quadrupole
./solver_benchmark

# Visualize (from project root)
cd ..
pip install plotly numpy
python scripts/visualize.py build/quadrupole.json -o plots/ -p quadrupole
python scripts/visualize.py build/coaxial.json -o plots/ -p coaxial
python scripts/visualize.py --benchmark build/benchmark.json -o plots/
```

## Visualization Examples

The Python script generates interactive HTML plots:

- **Potential heatmap** with equipotential contour overlay
- **Electric field magnitude** (|E|) with vector arrows
- **Line cuts** through the potential along x and y axes
- **Convergence curves** comparing solver methods (log-scale residual vs iteration)
- **Benchmark bar charts** showing iterations and wall-clock time per solver

## Validation

Each scenario includes quantitative checks against known analytic solutions:

- **Parallel plates**: Linear potential between plates, uniform E-field
- **Coaxial cylinders**: Logarithmic radial potential φ(r) = V·ln(r/r₂)/ln(r₁/r₂)
- **Quadrupole**: Diagonal symmetry (φ(+x,+y) ≈ φ(−x,−y)) and vanishing center field

## Extending

Adding a new geometry is straightforward:

1. Write a setup function in `src/geometry/boundary.cpp`
2. Add a scenario preset in `src/geometry/scenarios.cpp`
3. Create an example in `examples/`

The solver interface is polymorphic — implement `SolverBase::iterate()` to add new methods (e.g., multigrid, conjugate gradient).

## References

- J.D. Jackson, *Classical Electrodynamics*, 3rd ed. (Wiley, 1999)
- CERN Accelerator School proceedings on beam optics and magnet design
- W.H. Press et al., *Numerical Recipes*, Ch. 20: Partial Differential Equations
- S. Russenschuck, *Field Computation for Accelerator Magnets* (Wiley, 2010)

## License

MIT
