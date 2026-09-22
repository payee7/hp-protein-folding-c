/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "core/protein.h"
#include <string.h>

static int grid_occupancy[GRID_SIZE][GRID_SIZE];
static int grid_occupancy_3d[GRID_SIZE_3D][GRID_SIZE_3D][GRID_SIZE_3D];

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

static bool build_coords_from_directions_2d(Conformation *conf) {
    conf->dimension = 2;
    conf->coords[0].x = GRID_OFFSET;
    conf->coords[0].y = GRID_OFFSET;
    conf->coords[0].z = 0;
    
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 1;
    
    for (int i = 0; i < conf->length - 1; i++) {
        int dir = conf->directions[i];
        if (dir < 0 || dir >= 4) {
            dir = 0;
        }
        int nx = conf->coords[i].x + DX[dir];
        int ny = conf->coords[i].y + DY[dir];
        
        if (nx < 0 || nx >= GRID_SIZE || ny < 0 || ny >= GRID_SIZE || grid_occupancy[nx][ny] != 0) {
            for (int k = 0; k <= i; k++) {
                grid_occupancy[conf->coords[k].x][conf->coords[k].y] = 0;
            }
            return false;
        }
        
        conf->coords[i + 1].x = nx;
        conf->coords[i + 1].y = ny;
        conf->coords[i + 1].z = 0;
        grid_occupancy[nx][ny] = i + 2;
    }
    
    for (int k = 0; k < conf->length; k++) {
        grid_occupancy[conf->coords[k].x][conf->coords[k].y] = 0;
    }
    return true;
}

bool build_coords_from_directions_3d(Conformation *conf) {
    conf->dimension = 3;
    conf->coords[0].x = GRID_OFFSET_3D;
    conf->coords[0].y = GRID_OFFSET_3D;
    conf->coords[0].z = GRID_OFFSET_3D;
    
    grid_occupancy_3d[GRID_OFFSET_3D][GRID_OFFSET_3D][GRID_OFFSET_3D] = 1;
    
    for (int i = 0; i < conf->length - 1; i++) {
        int dir = conf->directions[i];
        if (dir < 0 || dir >= 6) {
            dir = 0;
        }
        int nx = conf->coords[i].x + DX3D[dir];
        int ny = conf->coords[i].y + DY3D[dir];
        int nz = conf->coords[i].z + DZ3D[dir];
        
        if (nx < 0 || nx >= GRID_SIZE_3D || ny < 0 || ny >= GRID_SIZE_3D || nz < 0 || nz >= GRID_SIZE_3D ||
            grid_occupancy_3d[nx][ny][nz] != 0) {
            for (int k = 0; k <= i; k++) {
                grid_occupancy_3d[conf->coords[k].x][conf->coords[k].y][conf->coords[k].z] = 0;
            }
            return false;
        }
        
        conf->coords[i + 1].x = nx;
        conf->coords[i + 1].y = ny;
        conf->coords[i + 1].z = nz;
        grid_occupancy_3d[nx][ny][nz] = i + 2;
    }
    
    for (int k = 0; k < conf->length; k++) {
        grid_occupancy_3d[conf->coords[k].x][conf->coords[k].y][conf->coords[k].z] = 0;
    }
    return true;
}

bool build_coords_from_directions(Conformation *conf) {
    if (conf->dimension == 3) {
        return build_coords_from_directions_3d(conf);
    }
    return build_coords_from_directions_2d(conf);
}
