/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#ifndef PDB_EXPORTER_H
#define PDB_EXPORTER_H

#include "core/types.h"

void export_pdb(const char *filename, const ProteinSequence *prot, const Conformation *conf);

#endif // PDB_EXPORTER_H
