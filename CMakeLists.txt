cmake_minimum_required(VERSION 3.16)
project(PoissonLaplaceSolver
    VERSION 1.0.0
    DESCRIPTION "2D Poisson & Laplace field solver for accelerator physics"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_EXAMPLES "Build example executables" ON)

# ── Core library ──────────────────────────────────────────────
add_library(poisson_solver
    src/solvers/jacobi.cpp
    src/solvers/gauss_seidel.cpp
    src/solvers/sor.cpp
    src/solvers/solver_base.cpp
    src/geometry/grid.cpp
    src/geometry/boundary.cpp
    src/geometry/scenarios.cpp
    src/utils/field_export.cpp
)

target_include_directories(poisson_solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(poisson_solver PRIVATE
    -Wall -Wextra -O2
)

# ── Examples ──────────────────────────────────────────────────
if(BUILD_EXAMPLES)
    add_executable(parallel_plates examples/parallel_plates.cpp)
    target_link_libraries(parallel_plates PRIVATE poisson_solver)

    add_executable(coaxial_cylinders examples/coaxial_cylinders.cpp)
    target_link_libraries(coaxial_cylinders PRIVATE poisson_solver)

    add_executable(quadrupole examples/quadrupole.cpp)
    target_link_libraries(quadrupole PRIVATE poisson_solver)

    add_executable(solver_benchmark examples/solver_benchmark.cpp)
    target_link_libraries(solver_benchmark PRIVATE poisson_solver)
endif()

# ── Tests ─────────────────────────────────────────────────────
if(BUILD_TESTS)
    enable_testing()
    add_executable(test_solvers tests/test_solvers.cpp)
    target_link_libraries(test_solvers PRIVATE poisson_solver)
    add_test(NAME SolverTests COMMAND test_solvers)
endif()
