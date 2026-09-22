/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include "core/types.h"

Conformation solve_simulated_annealing(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats);
Conformation solve_simulated_annealing_3d(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats);

#endif // SIMULATED_ANNEALING_H
