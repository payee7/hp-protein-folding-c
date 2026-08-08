/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#ifndef PROTEIN_H
#define PROTEIN_H

#include "types.h"

void init_protein(ProteinSequence *prot, const char *seq_str);
bool build_coords_from_directions(Conformation *conf);

#endif // PROTEIN_H
