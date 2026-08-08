/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "io/pdb_exporter.h"
#include <stdio.h>

void export_pdb(const char *filename, const ProteinSequence *prot, const Conformation *conf) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Error: Unable to create output PDB file %s\n", filename);
        return;
    }
    
    fprintf(fp, "REMARK   1 HP LATTICE PROTEIN FOLDING SOLVER OUTPUT\n");
    fprintf(fp, "REMARK   2 SEQUENCE: %s\n", prot->sequence);
    fprintf(fp, "REMARK   3 ENERGY: %d, RADIUS OF GYRATION: %.3f\n", conf->energy, conf->radius_of_gyration);
    
    const double SPACING = 3.8;
    
    for (int i = 0; i < conf->length; i++) {
        double x = conf->coords[i].x * SPACING;
        double y = conf->coords[i].y * SPACING;
        double z = 0.0;
        
        const char *res_name = (prot->amino[i] == AMINO_H) ? "VAL" : "ALA";
        
        fprintf(fp, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n",
                i + 1, res_name, i + 1, x, y, z);
    }
    
    for (int i = 1; i < conf->length; i++) {
        fprintf(fp, "CONECT%5d%5d\n", i, i + 1);
    }
    
    fprintf(fp, "END\n");
    fclose(fp);
    printf("Saved PDB structure: %s\n", filename);
}
