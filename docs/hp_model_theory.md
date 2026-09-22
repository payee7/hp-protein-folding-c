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

The protein backbone is represented as a self-avoiding walk (SAW) on a discrete 2D square lattice ($z = 4$) or 3D simple cubic lattice ($z = 6$).
- In 2D, an internal residue has coordination number 4, with 2 directions occupied by backbone bonds, leaving up to 2 directions for non-covalent contacts.
- In 3D, coordination number is 6 ($z = 6$). Each internal residue has up to 4 directions available for topological contacts, enabling substantially deeper hydrophobic core packing ($E_{3D} < E_{2D}$) and lower Radius of Gyration ($R_{g, 3D} < R_{g, 2D}$).

### Energy Function

The free energy of a conformation is determined strictly by non-bonded topological contacts between hydrophobic residues:

$$E = - \sum_{i < j - 1} \delta(r_i, r_j) \cdot \mathbb{I}(\text{Amino}_i = H \land \text{Amino}_j = H)$$

Where $\delta(r_i, r_j) = 1$ if residues $i$ and $j$ occupy adjacent lattice coordinates ($||r_i - r_j||_1 = 1$), and $0$ otherwise.

## Global Minimum Core Formation

Thermodynamically, proteins minimize free energy by sequestering hydrophobic ($H$) residues into an interior compact core, shielded from water by polar ($P$) residues on the surface. In 3D, ground-state conformations achieve higher ground-state degeneracy and tighter sphere packing than is geometrically possible in 2D.
