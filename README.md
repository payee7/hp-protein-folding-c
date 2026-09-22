# HP Lattice Protein Folding Prototyping Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/Build-Passing-brightgreen.svg)](CMakeLists.txt)
[![Language](https://img.shields.io/badge/Language-C99%20%2F%20C%2B%2B11-blue.svg)](src/main.c)

A modular, zero-dependency computational biophysics engine written in C/C++. The toolkit provides exact **Branch-and-Bound** global minimum energy searches alongside stochastic **Metropolis Simulated Annealing** (Pivot & Rigid Subchain Monte Carlo) for 2D/3D Hydrophobic-Polar (HP) polymer lattice models.

Developed as an internal software prototyping toolkit for algorithm testing, code analysis, and high-performance polymer simulations.

---

## Directory Structure

```text
hp-protein-folding-c/
├── include/                   # Header definitions
│   ├── core/                  # Polymer sequence & grid coordinate types
│   ├── solvers/               # Exact B&B & Simulated Annealing solvers
│   ├── analytics/             # Energy scoring & biophysical metric functions
│   └── io/                    # ASCII terminal renderer & PDB exporters
├── src/                       # Source code modules
│   ├── core/                  # Sequence initialization & self-avoiding walk logic
│   ├── solvers/               # Branch-and-Bound & Pivot Monte Carlo algorithms
│   ├── analytics/             # Contact matrix energy & Radius of Gyration
│   ├── io/                    # Visualizers & structure file generators
│   └── main.c                 # CLI application entry point
├── tests/                     # Unit test suite
├── benchmarks/                # Dill benchmark suite harness
├── docs/                      # Architectural & biophysical theory documentation
├── CMakeLists.txt             # CMake build configuration
├── Makefile                   # GNU Makefile
├── LICENSE                    # MIT License
└── README.md                  # Project documentation
```

---

## Core Features

1. **Decoupled Architecture**: Clean header/source separation for sequence modeling, energy analytics, solvers, and exporters.
2. **Dual Optimization Solvers (2D & 3D)**:
   - **Exact Branch-and-Bound**: Finds global minimum energy states and ground-state degeneracy for 2D square ($z=4$) and 3D cubic ($z=6$) polymer chains. Includes 2D/3D symmetry breaking and dynamic energy lower-bounding ($\hat{E}(s)$ pruning).
   - **Metropolis Simulated Annealing**: Stochastic Pivot Monte Carlo search supporting single-bond updates and rigid-body subchain rotations (2D planar & 3D $SO(3)$ rotations about X, Y, Z axes).
3. **Biophysical Structural Analytics**:
   - **Radius of Gyration ($R_g$)**: Computes 2D/3D spatial polymer compactness and folding density.
   - **Hydrophobic Core Index**: Calculates topological non-sequential H-H contact energy $E = -N_{\text{contacts}}$ across 4 planar or 6 spatial nearest neighbors.
4. **Visualizers & Exporters**:
   - **ASCII 2D & 3D Terminal Render**: Visualizes amino acid spatial layouts, hydrophobic cores (`H`), polar shells (`P`), intra-plane links, and 3D layer-by-layer $Z$-plane projections.
   - **3D PDB Exporter (`.pdb`)**: Generates Standard Protein Data Bank files containing pseudo C-$\alpha$ atom coordinates compatible with PyMOL, VMD, and ChimeraX.

---

## Quick Start & Installation

### Cloning Repository
```bash
git clone https://github.com/payee7/hp-protein-folding-c.git
cd hp-protein-folding-c
```

### Build with GCC / G++
```bash
g++ -O3 -Wall -Iinclude src/core/*.c src/analytics/*.c src/solvers/*.c src/io/*.c src/main.c -o hp_protein_solver
```

### Build with Makefile
```bash
make            # Build main solver executable
make test       # Run comprehensive 2D & 3D unit tests
make benchmark  # Run 2D vs 3D thermodynamic benchmark suite
```

### Build with CMake
```bash
mkdir build && cd build
cmake ..
make
```

### Run Solver
```bash
./hp_protein_solver "HPHPPHHPHPPHPHHPPHPH"         # Run both 2D and 3D solvers
./hp_protein_solver "HPHPPHHPHPPHPHHPPHPH" --2d    # Run 2D square lattice only
./hp_protein_solver "HPHPPHHPHPPHPHHPPHPH" --3d    # Run 3D cubic lattice only
```

---

## Verification & Benchmark Results (Dill's 20-mer)

- **Global Minimum Energy ($E$)**: `-9` (9 H-H Topological Contacts)
- **Search Space Explored**: 68,536 states (41,234 pruned via energy lower-bounding)
- **B&B Execution Time**: **0.0020 s** (2 ms)
- **Radius of Gyration ($R_g$)**: 1.844 grid units

```text
==== ASCII FOLDED CONFORMATION (2D HP GRID) ====
  P-H-P P-P
  |   | | |
  P-H H-H H
    |     |
    P-H H-P
      |
    P-H H
    |   |
    P-H-P
```

---

## License

This project is licensed under the [MIT License](LICENSE). Copyright (c) 2026 payee7.
