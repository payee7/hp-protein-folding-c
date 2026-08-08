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
    printf("             HP LATTICE SOLVER BENCHMARK SUITE                   \n");
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
        
        SolverStats stats;
        solve_exact_bb(&prot, &stats);
        
        printf("Benchmark %d [%d-mer]: Min Energy = %d, States Explored = %llu, Time = %.4fs\n",
               i + 1, prot.length, stats.best_energy, stats.states_explored, stats.elapsed_time_sec);
    }
    
    return 0;
}
