/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#ifndef BRANCH_BOUND_H
#define BRANCH_BOUND_H

#include "core/types.h"

Conformation solve_exact_bb(const ProteinSequence *prot, SolverStats *stats);
Conformation solve_exact_bb_3d(const ProteinSequence *prot, SolverStats *stats);

#endif // BRANCH_BOUND_H
