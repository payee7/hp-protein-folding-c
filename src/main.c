/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "core/protein.h"
#include "solvers/branch_bound.h"
#include "solvers/simulated_annealing.h"
#include "io/ascii_render.h"
#include "io/pdb_exporter.h"

int main(int argc, char *argv[]) {
    srand((unsigned int)time(NULL));
    
    printf("=================================================================\n");
    printf("        C/C++ HP LATTICE PROTEIN FOLDING SOLVER ENGINE           \n");
    printf("=================================================================\n\n");
    
    const char *benchmarks[] = {
        "HPHPPHHPHPPHPHHPPHPH",                          // 20-mer (Ground State E = -9)
        "HHHPPHPHPHPPHPHPHPPH",                          // 20-mer 2
        "PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPH",  // 44-mer
        "HHPPHPPHPPHPPHPPHPPHPPH"                        // 23-mer
    };
    
    const char *seq_input = benchmarks[0];
    if (argc > 1) {
        seq_input = argv[1];
    }
    
    ProteinSequence prot;
    init_protein(&prot, seq_input);
    
    printf("Target Sequence : %s\n", prot.sequence);
    printf("Length          : %d residues\n", prot.length);
    printf("Hydrophobic (H) : %d (%.1f%%)\n", prot.h_count, 100.0 * prot.h_count / prot.length);
    printf("Polar (P)       : %d\n", prot.length - prot.h_count);
    printf("-----------------------------------------------------------------\n\n");
    
    SolverStats stats_bb;
    Conformation best_bb;
    
    if (prot.length <= 25) {
        printf("[Solver 1] Exact Branch-and-Bound (Bounding + Symmetry Pruning)...\n");
        best_bb = solve_exact_bb(&prot, &stats_bb);
        
        printf("  > States Explored    : %llu\n", stats_bb.states_explored);
        printf("  > Energy Bound Pruned: %llu\n", stats_bb.states_pruned_bounds);
        printf("  > Collision Pruned   : %llu\n", stats_bb.states_pruned_collision);
        printf("  > Global Min Energy  : %d (H-H Contacts = %d)\n", stats_bb.best_energy, -stats_bb.best_energy);
        printf("  > Ground State Count : %d degenerate conformations\n", stats_bb.ground_state_count);
        printf("  > Execution Time     : %.4f seconds\n", stats_bb.elapsed_time_sec);
        printf("  > Radius of Gyration : %.3f grid units\n", best_bb.radius_of_gyration);
        
        render_ascii_folded_protein(&prot, &best_bb);
        export_pdb("protein_folded_exact.pdb", &prot, &best_bb);
    } else {
        printf("[Solver 1] Length %d > 25 (Exceeds exact B&B fast threshold). Skipping.\n", prot.length);
    }
    
    printf("\n[Solver 2] Stochastic Metropolis Simulated Annealing (Pivot MC)...\n");
    SolverStats stats_sa;
    int sa_steps = (prot.length > 30) ? 2000000 : 800000;
    Conformation best_sa = solve_simulated_annealing(&prot, sa_steps, 8.0, 0.999995, &stats_sa);
    
    printf("  > SA Steps Evaluated : %llu\n", stats_sa.states_explored);
    printf("  > Collision Pruned   : %llu\n", stats_sa.states_pruned_collision);
    printf("  > Best Energy Found  : %d (H-H Contacts = %d)\n", stats_sa.best_energy, -stats_sa.best_energy);
    printf("  > Execution Time     : %.4f seconds\n", stats_sa.elapsed_time_sec);
    printf("  > Radius of Gyration : %.3f grid units\n", best_sa.radius_of_gyration);
    
    render_ascii_folded_protein(&prot, &best_sa);
    export_pdb("protein_folded_sa.pdb", &prot, &best_sa);
    
    printf("\n=================================================================\n");
    printf("                     SOLVER RUN COMPLETED                        \n");
    printf("=================================================================\n");
    
    return 0;
}
