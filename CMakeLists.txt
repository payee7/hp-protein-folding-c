# ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
# Copyright (c) 2026 ProteinGo. Licensed under the MIT License.

cmake_minimum_required(VERSION 3.10)
project(HPProteinSolver C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

include_directories(include)

file(GLOB_RECURSE SOURCES "src/*.c")

add_executable(hp_protein_solver ${SOURCES})

# Unit test executable
add_executable(unit_tests 
    tests/test_energy.c 
    src/core/protein.c 
    src/analytics/energy.c 
    src/analytics/metrics.c
)

# Benchmark executable
add_executable(benchmark_runner 
    benchmarks/run_benchmarks.c 
    src/core/protein.c 
    src/analytics/energy.c 
    src/analytics/metrics.c
    src/solvers/branch_bound.c
    src/solvers/simulated_annealing.c
)
