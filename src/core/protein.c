/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "core/protein.h"
#include <string.h>

static int grid_occupancy[GRID_SIZE][GRID_SIZE];

void init_protein(ProteinSequence *prot, const char *seq_str) {
    prot->length = (int)strlen(seq_str);
    if (prot->length >= MAX_LEN) {
        prot->length = MAX_LEN - 1;
    }
    strncpy(prot->sequence, seq_str, (size_t)prot->length);
    prot->sequence[prot->length] = '\0';
    
    prot->h_count = 0;
    for (int i = 0; i < prot->length; i++) {
        if (prot->sequence[i] == 'H' || prot->sequence[i] == 'h') {
            prot->amino[i] = AMINO_H;
            prot->h_indices[prot->h_count] = i;
            prot->h_count++;
        } else {
            prot->amino[i] = AMINO_P;
        }
    }
    
    int count = 0;
    for (int i = prot->length - 1; i >= 0; i--) {
        if (prot->amino[i] == AMINO_H) count++;
        prot->remaining_h_from[i] = count;
    }
}

bool build_coords_from_directions(Conformation *conf) {
    conf->coords[0].x = GRID_OFFSET;
    conf->coords[0].y = GRID_OFFSET;
    
    memset(grid_occupancy, 0, sizeof(grid_occupancy));
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 1;
    
    for (int i = 0; i < conf->length - 1; i++) {
        int dir = conf->directions[i];
        int nx = conf->coords[i].x + DX[dir];
        int ny = conf->coords[i].y + DY[dir];
        
        if (grid_occupancy[nx][ny] != 0) {
            memset(grid_occupancy, 0, sizeof(grid_occupancy));
            return false;
        }
        
        conf->coords[i + 1].x = nx;
        conf->coords[i + 1].y = ny;
        grid_occupancy[nx][ny] = i + 2;
    }
    
    memset(grid_occupancy, 0, sizeof(grid_occupancy));
    return true;
}
