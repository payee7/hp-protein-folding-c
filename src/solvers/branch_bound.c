/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "solvers/branch_bound.h"
#include "analytics/metrics.h"
#include <string.h>
#include <time.h>

static int grid_occupancy[GRID_SIZE][GRID_SIZE];

static void bb_search(const ProteinSequence *prot, Conformation *current, int idx, 
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
    
    if (-max_possible_contacts >= stats->best_energy) {
        stats->states_pruned_bounds++;
        return;
    }
    
    for (int dir = 0; dir < 4; dir++) {
        if (idx == 0 && dir != 1) continue;
        if (idx == 1 && (dir == 2 || dir == 3)) continue;
        
        int nx = current->coords[idx].x + DX[dir];
        int ny = current->coords[idx].y + DY[dir];
        
        if (grid_occupancy[nx][ny] != 0) {
            stats->states_pruned_collision++;
            continue;
        }
        
        int new_contacts = 0;
        if (prot->amino[idx + 1] == AMINO_H) {
            for (int d = 0; d < 4; d++) {
                int adj_x = nx + DX[d];
                int adj_y = ny + DY[d];
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
        current->directions[idx] = dir;
        grid_occupancy[nx][ny] = idx + 2;
        
        bb_search(prot, current, idx + 1, current_contacts + new_contacts, stats, best_conf);
        
        grid_occupancy[nx][ny] = 0;
    }
}

Conformation solve_exact_bb(const ProteinSequence *prot, SolverStats *stats) {
    Conformation current;
    Conformation best_conf;
    
    current.length = prot->length;
    best_conf.length = prot->length;
    stats->states_explored = 0;
    stats->states_pruned_bounds = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = 1;
    stats->ground_state_count = 0;
    
    current.coords[0].x = GRID_OFFSET;
    current.coords[0].y = GRID_OFFSET;
    
    memset(grid_occupancy, 0, sizeof(grid_occupancy));
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 1;
    
    clock_t start_t = clock();
    bb_search(prot, &current, 0, 0, stats, &best_conf);
    clock_t end_t = clock();
    
    stats->elapsed_time_sec = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    memset(grid_occupancy, 0, sizeof(grid_occupancy));
    return best_conf;
}
