/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "solvers/branch_bound.h"
#include "analytics/metrics.h"
#include <string.h>
#include <time.h>

static int grid_occupancy[GRID_SIZE][GRID_SIZE];
static int grid_occupancy_3d[GRID_SIZE_3D][GRID_SIZE_3D][GRID_SIZE_3D];

static void bb_search_2d(const ProteinSequence *prot, Conformation *current, int idx, 
                         int current_contacts, SolverStats *stats, Conformation *best_conf) {
    stats->states_explored++;
    
    int len = prot->length;
    if (idx == len - 1) {
        int energy = -current_contacts;
        current->energy = energy;
        current->radius_of_gyration = calculate_radius_of_gyration(current);
        
        if (energy < stats->best_energy) {
            stats->best_energy = energy;
            stats->ground_state_count = 1;
            *best_conf = *current;
        } else if (energy == stats->best_energy) {
            stats->ground_state_count++;
        }
        return;
    }
    
    int remaining_h = prot->remaining_h_from[idx + 1];
    int max_possible_contacts = current_contacts + 2 * remaining_h;
    
    // Prune branch when theoretical upper bound cannot reach best_energy
    if (-max_possible_contacts > stats->best_energy) {
        stats->states_pruned_bounds++;
        return;
    }
    
    for (int dir = 0; dir < 4; dir++) {
        if (idx == 0 && dir != 1) continue;
        if (idx == 1 && (dir == 2 || dir == 3)) continue;
        
        int nx = current->coords[idx].x + DX[dir];
        int ny = current->coords[idx].y + DY[dir];
        
        if (nx < 0 || nx >= GRID_SIZE || ny < 0 || ny >= GRID_SIZE) continue;
        
        if (grid_occupancy[nx][ny] != 0) {
            stats->states_pruned_collision++;
            continue;
        }
        
        int new_contacts = 0;
        if (prot->amino[idx + 1] == AMINO_H) {
            for (int d = 0; d < 4; d++) {
                int adj_x = nx + DX[d];
                int adj_y = ny + DY[d];
                if (adj_x < 0 || adj_x >= GRID_SIZE || adj_y < 0 || adj_y >= GRID_SIZE) continue;
                int neighbor_pos = grid_occupancy[adj_x][adj_y];
                if (neighbor_pos > 0) {
                    int neighbor_idx = neighbor_pos - 1;
                    if (neighbor_idx != idx && prot->amino[neighbor_idx] == AMINO_H) {
                        new_contacts++;
                    }
                }
            }
        }
        
        current->coords[idx + 1].x = nx;
        current->coords[idx + 1].y = ny;
        current->coords[idx + 1].z = 0;
        current->directions[idx] = dir;
        grid_occupancy[nx][ny] = idx + 2;
        
        bb_search_2d(prot, current, idx + 1, current_contacts + new_contacts, stats, best_conf);
        
        grid_occupancy[nx][ny] = 0;
    }
}

Conformation solve_exact_bb(const ProteinSequence *prot, SolverStats *stats) {
    Conformation current;
    Conformation best_conf;
    
    current.length = prot->length;
    current.dimension = 2;
    best_conf.length = prot->length;
    best_conf.dimension = 2;
    stats->states_explored = 0;
    stats->states_pruned_bounds = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = 1;
    stats->ground_state_count = 0;
    
    current.coords[0].x = GRID_OFFSET;
    current.coords[0].y = GRID_OFFSET;
    current.coords[0].z = 0;
    
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 1;
    
    clock_t start_t = clock();
    bb_search_2d(prot, &current, 0, 0, stats, &best_conf);
    clock_t end_t = clock();
    
    stats->elapsed_time_sec = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 0;
    return best_conf;
}

static void bb_search_3d(const ProteinSequence *prot, Conformation *current, int idx, 
                         int current_contacts, SolverStats *stats, Conformation *best_conf) {
    stats->states_explored++;
    
    int len = prot->length;
    if (idx == len - 1) {
        int energy = -current_contacts;
        current->energy = energy;
        current->radius_of_gyration = calculate_radius_of_gyration(current);
        
        if (energy < stats->best_energy) {
            stats->best_energy = energy;
            stats->ground_state_count = 1;
            *best_conf = *current;
        } else if (energy == stats->best_energy) {
            stats->ground_state_count++;
        }
        return;
    }
    
    int remaining_h = prot->remaining_h_from[idx + 1];
    // In 3D (coordination number z = 6), an internal residue has at most 4 non-sequential contacts
    int max_possible_contacts = current_contacts + 4 * remaining_h;
    
    if (-max_possible_contacts > stats->best_energy) {
        stats->states_pruned_bounds++;
        return;
    }
    
    for (int dir = 0; dir < 6; dir++) {
        // 3D Cubic Symmetry Breaking:
        // 1. Move 0: Fixed to +X (dir 1)
        if (idx == 0 && dir != 1) continue;
        // 2. Move 1: Rotational symmetry around X axis restricts turn to +X (1) or +Y (0)
        if (idx == 1 && (dir != 1 && dir != 0)) continue;
        // 3. Move 2: If previous turn was in XY plane, reflection symmetry across XY plane restricts Z to +Z (dir 4)
        if (idx == 2 && current->directions[1] == 0 && dir == 5) continue;
        
        int nx = current->coords[idx].x + DX3D[dir];
        int ny = current->coords[idx].y + DY3D[dir];
        int nz = current->coords[idx].z + DZ3D[dir];
        
        if (nx < 0 || nx >= GRID_SIZE_3D || ny < 0 || ny >= GRID_SIZE_3D || nz < 0 || nz >= GRID_SIZE_3D) {
            continue;
        }
        
        if (grid_occupancy_3d[nx][ny][nz] != 0) {
            stats->states_pruned_collision++;
            continue;
        }
        
        int new_contacts = 0;
        if (prot->amino[idx + 1] == AMINO_H) {
            for (int d = 0; d < 6; d++) {
                int adj_x = nx + DX3D[d];
                int adj_y = ny + DY3D[d];
                int adj_z = nz + DZ3D[d];
                if (adj_x < 0 || adj_x >= GRID_SIZE_3D || adj_y < 0 || adj_y >= GRID_SIZE_3D || adj_z < 0 || adj_z >= GRID_SIZE_3D) {
                    continue;
                }
                int neighbor_pos = grid_occupancy_3d[adj_x][adj_y][adj_z];
                if (neighbor_pos > 0) {
                    int neighbor_idx = neighbor_pos - 1;
                    if (neighbor_idx != idx && prot->amino[neighbor_idx] == AMINO_H) {
                        new_contacts++;
                    }
                }
            }
        }
        
        current->coords[idx + 1].x = nx;
        current->coords[idx + 1].y = ny;
        current->coords[idx + 1].z = nz;
        current->directions[idx] = dir;
        grid_occupancy_3d[nx][ny][nz] = idx + 2;
        
        bb_search_3d(prot, current, idx + 1, current_contacts + new_contacts, stats, best_conf);
        
        grid_occupancy_3d[nx][ny][nz] = 0;
    }
}

Conformation solve_exact_bb_3d(const ProteinSequence *prot, SolverStats *stats) {
    Conformation current;
    Conformation best_conf;
    
    current.length = prot->length;
    current.dimension = 3;
    best_conf.length = prot->length;
    best_conf.dimension = 3;
    stats->states_explored = 0;
    stats->states_pruned_bounds = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = 1;
    stats->ground_state_count = 0;
    
    current.coords[0].x = GRID_OFFSET_3D;
    current.coords[0].y = GRID_OFFSET_3D;
    current.coords[0].z = GRID_OFFSET_3D;
    
    grid_occupancy_3d[GRID_OFFSET_3D][GRID_OFFSET_3D][GRID_OFFSET_3D] = 1;
    
    clock_t start_t = clock();
    bb_search_3d(prot, &current, 0, 0, stats, &best_conf);
    clock_t end_t = clock();
    
    stats->elapsed_time_sec = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    grid_occupancy_3d[GRID_OFFSET_3D][GRID_OFFSET_3D][GRID_OFFSET_3D] = 0;
    return best_conf;
}
