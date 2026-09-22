/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "solvers/simulated_annealing.h"
#include "core/protein.h"
#include "analytics/energy.h"
#include "analytics/metrics.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

// 3D Rotation lookup tables for 90-degree rotations about X, Y, Z axes
static const int ROT_X_POS[6] = { 4, 1, 5, 3, 2, 0 };
static const int ROT_X_NEG[6] = { 5, 1, 4, 3, 0, 2 };
static const int ROT_Y_POS[6] = { 0, 5, 2, 4, 1, 3 };
static const int ROT_Y_NEG[6] = { 0, 4, 2, 5, 3, 1 };
static const int ROT_Z_POS[6] = { 3, 0, 1, 2, 4, 5 };
static const int ROT_Z_NEG[6] = { 1, 2, 3, 0, 4, 5 };

Conformation solve_simulated_annealing(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats) {
    Conformation current;
    current.length = prot->length;
    current.dimension = 2;
    
    for (int i = 0; i < prot->length - 1; i++) {
        current.directions[i] = 1; // All RIGHT
    }
    build_coords_from_directions(&current);
    current.energy = calculate_energy(prot, &current);
    current.radius_of_gyration = calculate_radius_of_gyration(&current);
    
    Conformation best_conf = current;
    
    stats->states_explored = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = current.energy;
    
    if (prot->length < 3) {
        stats->elapsed_time_sec = 0.0;
        return best_conf;
    }
    
    double temp = initial_temp;
    clock_t start_t = clock();
    
    for (int step = 0; step < max_steps; step++) {
        stats->states_explored++;
        
        int pivot_idx = rand() % (prot->length - 1);
        
        Conformation candidate = current;
        int move_type = rand() % 3;
        
        if (move_type == 0) {
            candidate.directions[pivot_idx] = rand() % 4;
        } else if (move_type == 1) {
            for (int k = pivot_idx; k < prot->length - 1; k++) {
                candidate.directions[k] = (candidate.directions[k] + 1) % 4;
            }
        } else {
            for (int k = pivot_idx; k < prot->length - 1; k++) {
                candidate.directions[k] = (candidate.directions[k] + 3) % 4;
            }
        }
        
        if (!build_coords_from_directions(&candidate)) {
            stats->states_pruned_collision++;
            continue;
        }
        
        candidate.energy = calculate_energy(prot, &candidate);
        int delta_e = candidate.energy - current.energy;
        
        if (delta_e <= 0 || ((double)rand() / RAND_MAX) < exp(-delta_e / temp)) {
            current = candidate;
            if (current.energy < best_conf.energy) {
                best_conf = current;
                best_conf.radius_of_gyration = calculate_radius_of_gyration(&best_conf);
                stats->best_energy = best_conf.energy;
            }
        }
        
        temp *= cooling_rate;
    }
    
    clock_t end_t = clock();
    stats->elapsed_time_sec = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    
    return best_conf;
}

Conformation solve_simulated_annealing_3d(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats) {
    Conformation current;
    current.length = prot->length;
    current.dimension = 3;
    
    for (int i = 0; i < prot->length - 1; i++) {
        current.directions[i] = 1; // All +X
    }
    build_coords_from_directions_3d(&current);
    current.energy = calculate_energy_3d(prot, &current);
    current.radius_of_gyration = calculate_radius_of_gyration(&current);
    
    Conformation best_conf = current;
    
    stats->states_explored = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = current.energy;
    
    if (prot->length < 3) {
        stats->elapsed_time_sec = 0.0;
        return best_conf;
    }
    
    double temp = initial_temp;
    clock_t start_t = clock();
    
    for (int step = 0; step < max_steps; step++) {
        stats->states_explored++;
        
        int pivot_idx = rand() % (prot->length - 1);
        
        Conformation candidate = current;
        int move_type = rand() % 7;
        
        if (move_type == 0) {
            // Single bond 3D direction alteration
            candidate.directions[pivot_idx] = rand() % 6;
        } else {
            // Rigid subchain 3D rotation around X, Y, or Z axis
            const int *rot_map = NULL;
            switch (move_type) {
                case 1: rot_map = ROT_X_POS; break;
                case 2: rot_map = ROT_X_NEG; break;
                case 3: rot_map = ROT_Y_POS; break;
                case 4: rot_map = ROT_Y_NEG; break;
                case 5: rot_map = ROT_Z_POS; break;
                default: rot_map = ROT_Z_NEG; break;
            }
            for (int k = pivot_idx; k < prot->length - 1; k++) {
                int d = candidate.directions[k];
                if (d >= 0 && d < 6) {
                    candidate.directions[k] = rot_map[d];
                }
            }
        }
        
        if (!build_coords_from_directions_3d(&candidate)) {
            stats->states_pruned_collision++;
            continue;
        }
        
        candidate.energy = calculate_energy_3d(prot, &candidate);
        int delta_e = candidate.energy - current.energy;
        
        if (delta_e <= 0 || ((double)rand() / RAND_MAX) < exp(-delta_e / temp)) {
            current = candidate;
            if (current.energy < best_conf.energy) {
                best_conf = current;
                best_conf.radius_of_gyration = calculate_radius_of_gyration(&best_conf);
                stats->best_energy = best_conf.energy;
            }
        }
        
        temp *= cooling_rate;
    }
    
    clock_t end_t = clock();
    stats->elapsed_time_sec = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    
    return best_conf;
}
