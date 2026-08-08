/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "analytics/energy.h"
#include <string.h>

static int grid_occupancy[GRID_SIZE][GRID_SIZE];

int calculate_energy(const ProteinSequence *prot, const Conformation *conf) {
    int contacts = 0;
    
    for (int i = 0; i < conf->length; i++) {
        grid_occupancy[conf->coords[i].x][conf->coords[i].y] = i + 1;
    }
    
    for (int i = 0; i < conf->length; i++) {
        if (prot->amino[i] != AMINO_H) continue;
        
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        
        for (int d = 0; d < 4; d++) {
            int nx = x + DX[d];
            int ny = y + DY[d];
            
            int neighbor_idx = grid_occupancy[nx][ny] - 1;
            
            if (neighbor_idx >= 0 && neighbor_idx != i - 1 && neighbor_idx != i + 1) {
                if (neighbor_idx > i && prot->amino[neighbor_idx] == AMINO_H) {
                    contacts++;
                }
            }
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        grid_occupancy[conf->coords[i].x][conf->coords[i].y] = 0;
    }
    
    return -contacts;
}
