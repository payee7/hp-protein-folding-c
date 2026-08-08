<!--
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
-->

# Hydrophobic-Polar (HP) Lattice Model Theory

The **Hydrophobic-Polar (HP) model** (Dill et al., 1985) is a fundamental computational abstraction of protein folding thermodynamics.

## Thermodynamic Formulation

Primary amino acid sequences are simplified into a binary alphabet:
- **H (Hydrophobic)**: Valine, Leucine, Isoleucine, Phenylalanine, etc.
- **P (Polar / Hydrophilic)**: Alanine, Glycine, Serine, etc.

The protein backbone is represented as a self-avoiding walk (SAW) on a discrete 2D/3D square/cubic lattice.

### Energy Function

The free energy of a conformation is determined strictly by non-bonded topological contacts between hydrophobic residues:

$$E = - \sum_{i < j - 1} \delta(r_i, r_j) \cdot \mathbb{I}(\text{Amino}_i = H \land \text{Amino}_j = H)$$

Where $\delta(r_i, r_j) = 1$ if residues $i$ and $j$ occupy adjacent lattice coordinates ($||r_i - r_j||_1 = 1$), and $0$ otherwise.

## Global Minimum Core Formation

Thermodynamically, proteins minimize free energy by sequestering hydrophobic ($H$) residues into an interior compact core, shielded from water by polar ($P$) residues on the surface.
