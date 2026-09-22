#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "core/protein.h"
#include "solvers/branch_bound.h"
#include "solvers/simulated_annealing.h"
#include "io/ascii_render.h"
#include "io/pdb_exporter.h"

int main(int argc, char *argv[]) {
    srand((unsigned int)time(NULL));
    
    printf("=================================================================\n");
    printf("       C/C++ HP LATTICE PROTEIN FOLDING SOLVER ENGINE (2D/3D)    \n");
    printf("=================================================================\n\n");
    
    const char *benchmarks[] = {
        "HPHPPHHPHPPHPHHPPHPH",                          // 20-mer (Ground State 2D E = -9)
        "HHHPPHPHPHPPHPHPHPPH",                          // 20-mer 2
        "PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPH",  // 44-mer
        "HHPPHPPHPPHPPHPPHPPHPPH"                        // 23-mer
    };
    
    const char *seq_input = benchmarks[0];
    int run_dim = 0; // 0: both, 2: 2D only, 3: 3D only
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--2d") == 0) {
            run_dim = 2;
        } else if (strcmp(argv[i], "--3d") == 0) {
            run_dim = 3;
        } else if (argv[i][0] != '-') {
            seq_input = argv[i];
        }
    }
    
    ProteinSequence prot;
    init_protein(&prot, seq_input);
    
    printf("Target Sequence : %s\n", prot.sequence);
    printf("Length          : %d residues\n", prot.length);
    printf("Hydrophobic (H) : %d (%.1f%%)\n", prot.h_count, 100.0 * prot.h_count / prot.length);
    printf("Polar (P)       : %d\n", prot.length - prot.h_count);
    printf("-----------------------------------------------------------------\n\n");
    
    // =========================================================================
    // 2D LATTICE SOLVER
    // =========================================================================
    if (run_dim == 0 || run_dim == 2) {
        printf(">>>>>>>>>>>>>>>>>>>> [ 2D SQUARE LATTICE (z = 4) ] >>>>>>>>>>>>>>>>>>>>\n");
        
        SolverStats stats_bb2d;
        Conformation best_bb2d;
        if (prot.length <= 25) {
            printf("[2D Solver 1] Exact Branch-and-Bound...\n");
            best_bb2d = solve_exact_bb(&prot, &stats_bb2d);
            
            printf("  > States Explored    : %llu\n", stats_bb2d.states_explored);
            printf("  > Energy Bound Pruned: %llu\n", stats_bb2d.states_pruned_bounds);
            printf("  > Collision Pruned   : %llu\n", stats_bb2d.states_pruned_collision);
            printf("  > Global Min Energy  : %d (H-H Contacts = %d)\n", stats_bb2d.best_energy, -stats_bb2d.best_energy);
            printf("  > Ground State Count : %d degenerate conformations\n", stats_bb2d.ground_state_count);
            printf("  > Execution Time     : %.4f seconds\n", stats_bb2d.elapsed_time_sec);
            printf("  > Radius of Gyration : %.3f grid units\n", best_bb2d.radius_of_gyration);
            
            render_ascii_folded_protein(&prot, &best_bb2d);
            export_pdb("protein_folded_exact_2d.pdb", &prot, &best_bb2d);
        } else {
            printf("[2D Solver 1] Length %d > 25 (Exceeds exact B&B fast threshold). Skipping.\n", prot.length);
        }
        
        printf("\n[2D Solver 2] Stochastic Metropolis Simulated Annealing (Pivot MC)...\n");
        SolverStats stats_sa2d;
        int sa_steps = (prot.length > 30) ? 2000000 : 800000;
        Conformation best_sa2d = solve_simulated_annealing(&prot, sa_steps, 8.0, 0.999995, &stats_sa2d);
        
        printf("  > SA Steps Evaluated : %llu\n", stats_sa2d.states_explored);
        printf("  > Collision Pruned   : %llu\n", stats_sa2d.states_pruned_collision);
        printf("  > Best Energy Found  : %d (H-H Contacts = %d)\n", stats_sa2d.best_energy, -stats_sa2d.best_energy);
        printf("  > Execution Time     : %.4f seconds\n", stats_sa2d.elapsed_time_sec);
        printf("  > Radius of Gyration : %.3f grid units\n", best_sa2d.radius_of_gyration);
        
        export_pdb("protein_folded_sa_2d.pdb", &prot, &best_sa2d);
        printf("\n");
    }
    
    // =========================================================================
    // 3D CUBIC LATTICE SOLVER
    // =========================================================================
    if (run_dim == 0 || run_dim == 3) {
        printf(">>>>>>>>>>>>>>>>>>>> [ 3D CUBIC LATTICE (z = 6) ] >>>>>>>>>>>>>>>>>>>>\n");
        
        SolverStats stats_bb3d;
        Conformation best_bb3d;
        if (prot.length <= 18) {
            printf("[3D Solver 1] Exact Branch-and-Bound (3D Symmetry Pruning)...\n");
            best_bb3d = solve_exact_bb_3d(&prot, &stats_bb3d);
            
            printf("  > States Explored    : %llu\n", stats_bb3d.states_explored);
            printf("  > Energy Bound Pruned: %llu\n", stats_bb3d.states_pruned_bounds);
            printf("  > Collision Pruned   : %llu\n", stats_bb3d.states_pruned_collision);
            printf("  > Global Min Energy  : %d (H-H Contacts = %d)\n", stats_bb3d.best_energy, -stats_bb3d.best_energy);
            printf("  > Ground State Count : %d degenerate conformations\n", stats_bb3d.ground_state_count);
            printf("  > Execution Time     : %.4f seconds\n", stats_bb3d.elapsed_time_sec);
            printf("  > Radius of Gyration : %.3f grid units\n", best_bb3d.radius_of_gyration);
            
            render_ascii_folded_protein(&prot, &best_bb3d);
            export_pdb("protein_folded_exact_3d.pdb", &prot, &best_bb3d);
        } else {
            printf("[3D Solver 1] Length %d > 18 (Exceeds exact 3D B&B fast threshold). Skipping.\n", prot.length);
        }
        
        printf("\n[3D Solver 2] 3D Stochastic Simulated Annealing (SO(3) Rigid Subchain Rotations)...\n");
        SolverStats stats_sa3d;
        int sa_steps_3d = (prot.length > 30) ? 2500000 : 1000000;
        Conformation best_sa3d = solve_simulated_annealing_3d(&prot, sa_steps_3d, 10.0, 0.999996, &stats_sa3d);
        
        printf("  > 3D Steps Evaluated : %llu\n", stats_sa3d.states_explored);
        printf("  > Collision Pruned   : %llu\n", stats_sa3d.states_pruned_collision);
        printf("  > Best Energy Found  : %d (H-H Contacts = %d)\n", stats_sa3d.best_energy, -stats_sa3d.best_energy);
        printf("  > Execution Time     : %.4f seconds\n", stats_sa3d.elapsed_time_sec);
        printf("  > Radius of Gyration : %.3f grid units\n", best_sa3d.radius_of_gyration);
        
        render_ascii_folded_protein(&prot, &best_sa3d);
        export_pdb("protein_folded_sa_3d.pdb", &prot, &best_sa3d);
        printf("\n");
    }
    
    printf("=================================================================\n");
    printf("                     SOLVER RUN COMPLETED                        \n");
    printf("=================================================================\n");
    
    return 0;
}
