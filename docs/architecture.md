<!--
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
-->

# Modular System Architecture

The **HP Lattice Protein Folding Solver** is structured as a decoupled, multi-module computational engine:

```text
hp-protein-folding-c/
├── include/               # Public API headers
│   ├── core/              # Type definitions & sequence structures
│   ├── solvers/           # B&B & Simulated Annealing solvers
│   ├── analytics/         # Energy evaluation & biophysical metrics
│   └── io/                # ASCII terminal renderer & PDB exporters
├── src/                   # Implementation modules
│   ├── core/
│   ├── solvers/
│   ├── analytics/
│   └── io/
├── tests/                 # Modular unit test suites
├── benchmarks/            # Dill benchmark suite harness
└── docs/                  # Architectural & theoretical documentation
```

## Module Descriptions

1. **`core`**: Contains lattice coordinate representation, self-avoiding walk (SAW) spatial grid lookups, and sequence parsing.
2. **`analytics`**: Implements topological contact energy calculation $E = -N_{\text{contacts}}$ and Radius of Gyration ($R_g$).
3. **`solvers`**:
   - `branch_bound.c`: Exact solver with dynamic lower-bound pruning and 8-fold lattice symmetry breaking.
   - `simulated_annealing.c`: Stochastic Pivot & Rigid Subchain Monte Carlo engine.
4. **`io`**: ASCII 2D visualization and Standard PDB format export.
