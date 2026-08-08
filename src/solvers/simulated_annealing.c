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

Conformation solve_simulated_annealing(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats) {
    Conformation current;
    current.length = prot->length;
    
    for (int i = 0; i < prot->length - 1; i++) {
        current.directions[i] = 1;
    }
    build_coords_from_directions(&current);
    current.energy = calculate_energy(prot, &current);
    
    Conformation best_conf = current;
    
    double temp = initial_temp;
    stats->states_explored = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = current.energy;
    
    clock_t start_t = clock();
    
    for (int step = 0; step < max_steps; step++) {
        stats->states_explored++;
        
        int pivot_idx = 1 + rand() % (prot->length - 2);
        
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
