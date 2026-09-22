/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include <stdio.h>
#include <stdlib.h>
#include "core/protein.h"
#include "solvers/branch_bound.h"
#include "solvers/simulated_annealing.h"

int main() {
    printf("=================================================================\n");
    printf("        HP LATTICE SOLVER BENCHMARK SUITE (2D vs 3D)             \n");
    printf("=================================================================\n\n");
    
    const char *suite[] = {
        "HPHPPHHPHPPHPHHPPHPH",                          // 20-mer
        "HHHPPHPHPHPPHPHPHPPH",                          // 20-mer (2)
        "HHPPHPPHPPHPPHPPHPPHPPH"                        // 23-mer
    };
    int n_bench = 3;
    
    for (int i = 0; i < n_bench; i++) {
        ProteinSequence prot;
        init_protein(&prot, suite[i]);
        
        printf("Benchmark %d: Sequence [%d-mer] %s\n", i + 1, prot.length, prot.sequence);
        
        // 2D Exact B&B
        SolverStats stats_2d;
        solve_exact_bb(&prot, &stats_2d);
        printf("  [2D B&B] Min Energy = %2d (Contacts = %2d), States = %8llu, Time = %.4fs\n",
               stats_2d.best_energy, -stats_2d.best_energy, stats_2d.states_explored, stats_2d.elapsed_time_sec);
        
        // 3D Simulated Annealing
        SolverStats stats_3d;
        Conformation conf_3d = solve_simulated_annealing_3d(&prot, 1500000, 10.0, 0.999996, &stats_3d);
        printf("  [3D SA ] Min Energy = %2d (Contacts = %2d), States = %8llu, Time = %.4fs, Rg = %.3f\n\n",
               conf_3d.energy, -conf_3d.energy, stats_3d.states_explored, stats_3d.elapsed_time_sec, conf_3d.radius_of_gyration);
    }
    
    return 0;
}
