/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#ifndef ENERGY_H
#define ENERGY_H

#include "core/types.h"

int calculate_energy(const ProteinSequence *prot, const Conformation *conf);
int calculate_energy_3d(const ProteinSequence *prot, const Conformation *conf);

#endif // ENERGY_H
