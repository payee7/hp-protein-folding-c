/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "analytics/energy.h"
#include <string.h>

static int grid_occupancy[GRID_SIZE][GRID_SIZE];
static int grid_occupancy_3d[GRID_SIZE_3D][GRID_SIZE_3D][GRID_SIZE_3D];

static int calculate_energy_2d(const ProteinSequence *prot, const Conformation *conf) {
    int contacts = 0;
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        if (x >= 0 && x < GRID_SIZE && y >= 0 && y < GRID_SIZE) {
            grid_occupancy[x][y] = i + 1;
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        if (prot->amino[i] != AMINO_H) continue;
        
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        
        for (int d = 0; d < 4; d++) {
            int nx = x + DX[d];
            int ny = y + DY[d];
            
            if (nx < 0 || nx >= GRID_SIZE || ny < 0 || ny >= GRID_SIZE) continue;
            
            int neighbor_idx = grid_occupancy[nx][ny] - 1;
            
            if (neighbor_idx >= 0 && neighbor_idx != i - 1 && neighbor_idx != i + 1) {
                if (neighbor_idx > i && prot->amino[neighbor_idx] == AMINO_H) {
                    contacts++;
                }
            }
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        if (x >= 0 && x < GRID_SIZE && y >= 0 && y < GRID_SIZE) {
            grid_occupancy[x][y] = 0;
        }
    }
    
    return -contacts;
}

int calculate_energy_3d(const ProteinSequence *prot, const Conformation *conf) {
    int contacts = 0;
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        int z = conf->coords[i].z;
        if (x >= 0 && x < GRID_SIZE_3D && y >= 0 && y < GRID_SIZE_3D && z >= 0 && z < GRID_SIZE_3D) {
            grid_occupancy_3d[x][y][z] = i + 1;
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        if (prot->amino[i] != AMINO_H) continue;
        
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        int z = conf->coords[i].z;
        
        for (int d = 0; d < 6; d++) {
            int nx = x + DX3D[d];
            int ny = y + DY3D[d];
            int nz = z + DZ3D[d];
            
            if (nx < 0 || nx >= GRID_SIZE_3D || ny < 0 || ny >= GRID_SIZE_3D || nz < 0 || nz >= GRID_SIZE_3D) {
                continue;
            }
            
            int neighbor_idx = grid_occupancy_3d[nx][ny][nz] - 1;
            
            if (neighbor_idx >= 0 && neighbor_idx != i - 1 && neighbor_idx != i + 1) {
                if (neighbor_idx > i && prot->amino[neighbor_idx] == AMINO_H) {
                    contacts++;
                }
            }
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        int z = conf->coords[i].z;
        if (x >= 0 && x < GRID_SIZE_3D && y >= 0 && y < GRID_SIZE_3D && z >= 0 && z < GRID_SIZE_3D) {
            grid_occupancy_3d[x][y][z] = 0;
        }
    }
    
    return -contacts;
}

int calculate_energy(const ProteinSequence *prot, const Conformation *conf) {
    if (conf->dimension == 3) {
        return calculate_energy_3d(prot, conf);
    }
    return calculate_energy_2d(prot, conf);
}
