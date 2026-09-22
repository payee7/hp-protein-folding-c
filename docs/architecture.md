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

1. **`core`**: Contains lattice coordinate representation (`Point3D`), self-avoiding walk (SAW) spatial grid lookups (2D $512\times 512$ and 3D $128\times 128\times 128$), and sequence parsing.
2. **`analytics`**: Implements 2D/3D topological contact energy calculation $E = -N_{\text{contacts}}$ across 4 planar or 6 spatial neighbors, and 3D Radius of Gyration ($R_g$).
3. **`solvers`**:
   - `branch_bound.c`: Exact solvers (`solve_exact_bb` and `solve_exact_bb_3d`) with dynamic lower-bound pruning and 2D/3D lattice symmetry breaking.
   - `simulated_annealing.c`: Stochastic Pivot & Rigid Subchain Monte Carlo engine supporting 2D planar pivots and 3D $SO(3)$ cubic subchain rotations.
4. **`io`**: ASCII 2D and 3D layer-by-layer terminal visualization, and Standard PDB format export with pseudo $C_\alpha$ coordinates for PyMOL, ChimeraX, and VMD.
