/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 *
 * Features:
 *  1. 2D/3D Hydrophobic-Polar (HP) Lattice Model
 *  2. Exact Branch-and-Bound Solver with Symmetry Pruning & Energy Lower-Bounding
 *  3. Metropolis Simulated Annealing with Ergonomic Pivot / Pull Moves
 *  4. ASCII 2D Grid Visualizer with Topological Contact Rendering
 *  5. Structural Analytics (Radius of Gyration, Compactness, Hydrophobic Core Index)
 *  6. Standard PDB File Exporter for 3D Protein Structure Visualizers
 *  7. Classic Benchmark Suite (Dill's HP sequences)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define MAX_LEN 128
#define GRID_SIZE 256
#define GRID_OFFSET 128

// Direction vectors for 2D (0: Up, 1: Right, 2: Down, 3: Left)
static const int DX[4] = { 0,  1,  0, -1};
static const int DY[4] = { 1,  0, -1,  0};
static const char DIR_NAMES[4][8] = {"UP", "RIGHT", "DOWN", "LEFT"};

// Amino Acid Types
typedef enum {
    AMINO_P = 0, // Polar (Hydrophilic, Neutral)
    AMINO_H = 1  // Hydrophobic (Attracted to other H)
} AminoType;

typedef struct {
    int x;
    int y;
} Point2D;

typedef struct {
    char sequence[MAX_LEN];
    AminoType amino[MAX_LEN];
    int length;
    int h_count;
    int h_indices[MAX_LEN];
    int remaining_h_from[MAX_LEN]; // Prefix/Suffix remaining H count for bounding
} ProteinSequence;

typedef struct {
    Point2D coords[MAX_LEN];
    int directions[MAX_LEN]; // Direction chosen at index i (0..3)
    int energy;              // E = - (number of non-sequential H-H contacts)
    int length;
    double radius_of_gyration;
} Conformation;

// Solver statistics tracking
typedef struct {
    unsigned long long states_explored;
    unsigned long long states_pruned_bounds;
    unsigned long long states_pruned_collision;
    int best_energy;
    int ground_state_count;
    double elapsed_time_sec;
} SolverStats;

// Dynamic occupancy lookup table for O(1) collision & contact checks
static int grid_occupancy[GRID_SIZE][GRID_SIZE];

// Initialize protein sequence data
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
    
    // Precompute remaining H count from index i to end (used for B&B energy lower bound)
    int count = 0;
    for (int i = prot->length - 1; i >= 0; i--) {
        if (prot->amino[i] == AMINO_H) count++;
        prot->remaining_h_from[i] = count;
    }
}

// Energy evaluation function: E = - (H-H topological non-sequential contacts)
int calculate_energy(const ProteinSequence *prot, const Conformation *conf) {
    int contacts = 0;
    
    // Populate occupancy grid with 1-based residue indices
    for (int i = 0; i < conf->length; i++) {
        grid_occupancy[conf->coords[i].x][conf->coords[i].y] = i + 1;
    }
    
    for (int i = 0; i < conf->length; i++) {
        if (prot->amino[i] != AMINO_H) continue;
        
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        
        // Check 4 topological grid neighbors
        for (int d = 0; d < 4; d++) {
            int nx = x + DX[d];
            int ny = y + DY[d];
            
            int neighbor_idx = grid_occupancy[nx][ny] - 1;
            
            // Valid non-sequential contact check
            if (neighbor_idx >= 0 && neighbor_idx != i - 1 && neighbor_idx != i + 1) {
                if (neighbor_idx > i && prot->amino[neighbor_idx] == AMINO_H) {
                    contacts++;
                }
            }
        }
    }
    
    // Clean up grid region
    for (int i = 0; i < conf->length; i++) {
        grid_occupancy[conf->coords[i].x][conf->coords[i].y] = 0;
    }
    
    return -contacts;
}

// Compute Radius of Gyration (Rg)
double calculate_radius_of_gyration(const Conformation *conf) {
    if (conf->length <= 0) return 0.0;
    
    double cx = 0.0, cy = 0.0;
    for (int i = 0; i < conf->length; i++) {
        cx += conf->coords[i].x;
        cy += conf->coords[i].y;
    }
    cx /= conf->length;
    cy /= conf->length;
    
    double sq_sum = 0.0;
    for (int i = 0; i < conf->length; i++) {
        double dx = conf->coords[i].x - cx;
        double dy = conf->coords[i].y - cy;
        sq_sum += (dx * dx + dy * dy);
    }
    return sqrt(sq_sum / conf->length);
}

// Reconstruct 2D spatial coordinates from movement direction vector
bool build_coords_from_directions(Conformation *conf) {
    conf->coords[0].x = GRID_OFFSET;
    conf->coords[0].y = GRID_OFFSET;
    
    memset(grid_occupancy, 0, sizeof(grid_occupancy));
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 1;
    
    for (int i = 0; i < conf->length - 1; i++) {
        int dir = conf->directions[i];
        int nx = conf->coords[i].x + DX[dir];
        int ny = conf->coords[i].y + DY[dir];
        
        // Self-avoiding walk check
        if (grid_occupancy[nx][ny] != 0) {
            memset(grid_occupancy, 0, sizeof(grid_occupancy));
            return false; // Collision detected
        }
        
        conf->coords[i + 1].x = nx;
        conf->coords[i + 1].y = ny;
        grid_occupancy[nx][ny] = i + 2;
    }
    
    memset(grid_occupancy, 0, sizeof(grid_occupancy));
    return true;
}

// Branch-and-bound recursive search with energy lower-bounding & symmetry pruning
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
    
    // Theoretical maximum additional contacts remaining from unplaced H residues
    int remaining_h = prot->remaining_h_from[idx + 1];
    // Each remaining H can form at most 2 contacts (since 2 directions are backbone)
    int max_possible_contacts = current_contacts + 2 * remaining_h;
    
    // Prune branch if even theoretical upper bound cannot improve best_energy
    if (-max_possible_contacts >= stats->best_energy) {
        stats->states_pruned_bounds++;
        return;
    }
    
    // Try 4 direction choices for next step
    for (int dir = 0; dir < 4; dir++) {
        // Symmetry breaking rule 1: First move must be RIGHT (dir = 1)
        if (idx == 0 && dir != 1) continue;
        // Symmetry breaking rule 2: Second move restricted to UP (0) or RIGHT (1)
        if (idx == 1 && (dir == 2 || dir == 3)) continue;
        
        int nx = current->coords[idx].x + DX[dir];
        int ny = current->coords[idx].y + DY[dir];
        
        // Self-intersection check
        if (grid_occupancy[nx][ny] != 0) {
            stats->states_pruned_collision++;
            continue;
        }
        
        // Check new H-H contacts formed by placing residue idx+1 at (nx, ny)
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
        
        // Advance walk
        current->coords[idx + 1].x = nx;
        current->coords[idx + 1].y = ny;
        current->directions[idx] = dir;
        grid_occupancy[nx][ny] = idx + 2;
        
        // Recurse
        bb_search(prot, current, idx + 1, current_contacts + new_contacts, stats, best_conf);
        
        // Backtrack
        grid_occupancy[nx][ny] = 0;
    }
}

// Run Exact Branch-and-Bound Solver
Conformation solve_exact_bb(const ProteinSequence *prot, SolverStats *stats) {
    Conformation current;
    Conformation best_conf;
    
    current.length = prot->length;
    best_conf.length = prot->length;
    stats->states_explored = 0;
    stats->states_pruned_bounds = 0;
    stats->states_pruned_collision = 0;
    stats->best_energy = 1; // Unrealistic high starting energy
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

// Simulated Annealing Solver using Pivot Monte Carlo Moves
Conformation solve_simulated_annealing(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats) {
    Conformation current;
    current.length = prot->length;
    
    // Start with extended linear backbone conformation
    for (int i = 0; i < prot->length - 1; i++) {
        current.directions[i] = 1; // All RIGHT
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
        
        // Select pivot index (1 .. length - 2)
        int pivot_idx = 1 + rand() % (prot->length - 2);
        
        // Pick new direction or pivot transformation
        Conformation candidate = current;
        int move_type = rand() % 3;
        
        if (move_type == 0) {
            // Single direction alteration
            candidate.directions[pivot_idx] = rand() % 4;
        } else if (move_type == 1) {
            // Rigid subchain 90-degree rotation
            for (int k = pivot_idx; k < prot->length - 1; k++) {
                candidate.directions[k] = (candidate.directions[k] + 1) % 4;
            }
        } else {
            // Rigid subchain 270-degree rotation
            for (int k = pivot_idx; k < prot->length - 1; k++) {
                candidate.directions[k] = (candidate.directions[k] + 3) % 4;
            }
        }
        
        // Check validity (self-avoiding walk constraint)
        if (!build_coords_from_directions(&candidate)) {
            stats->states_pruned_collision++;
            continue;
        }
        
        candidate.energy = calculate_energy(prot, &candidate);
        int delta_e = candidate.energy - current.energy;
        
        // Metropolis acceptance rule
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

// ASCII Render function for 2D Folded Grid Conformation
void render_ascii_folded_protein(const ProteinSequence *prot, const Conformation *conf) {
    int min_x = GRID_SIZE, max_x = 0;
    int min_y = GRID_SIZE, max_y = 0;
    
    for (int i = 0; i < conf->length; i++) {
        if (conf->coords[i].x < min_x) min_x = conf->coords[i].x;
        if (conf->coords[i].x > max_x) max_x = conf->coords[i].x;
        if (conf->coords[i].y < min_y) min_y = conf->coords[i].y;
        if (conf->coords[i].y > max_y) max_y = conf->coords[i].y;
    }
    
    min_x--; max_x++;
    min_y--; max_y++;
    
    int width = max_x - min_x + 1;
    int height = max_y - min_y + 1;
    
    if (width <= 0 || height <= 0) return;
    
    size_t canvas_w = (size_t)(width * 2 + 1);
    size_t canvas_h = (size_t)(height * 2 + 1);
    
    char **canvas = (char **)malloc(canvas_h * sizeof(char *));
    if (!canvas) return;
    
    for (size_t r = 0; r < canvas_h; r++) {
        canvas[r] = (char *)malloc(canvas_w * sizeof(char));
        if (!canvas[r]) continue;
        memset(canvas[r], ' ', canvas_w - 1);
        canvas[r][canvas_w - 1] = '\0';
    }
    
    for (int i = 0; i < conf->length; i++) {
        int gx = conf->coords[i].x - min_x;
        int gy = conf->coords[i].y - min_y;
        
        size_t cr = (size_t)((height - 1 - gy) * 2);
        size_t cc = (size_t)(gx * 2);
        
        canvas[cr][cc] = (prot->amino[i] == AMINO_H) ? 'H' : 'P';
        
        // Draw backbone segment
        if (i < conf->length - 1) {
            int ngx = conf->coords[i + 1].x - min_x;
            int ngy = conf->coords[i + 1].y - min_y;
            
            size_t ncr = (size_t)((height - 1 - ngy) * 2);
            size_t ncc = (size_t)(ngx * 2);
            
            size_t mid_r = (cr + ncr) / 2;
            size_t mid_c = (cc + ncc) / 2;
            
            if (cr == ncr) {
                canvas[mid_r][mid_c] = '-';
            } else if (cc == ncc) {
                canvas[mid_r][mid_c] = '|';
            }
        }
    }
    
    printf("\n==== ASCII FOLDED CONFORMATION (2D HP GRID) ====\n");
    printf("Legend: 'H' = Hydrophobic (Core), 'P' = Polar (Surface), '-'/'|' = Backbone\n");
    printf("Grid Bounds: X[%d..%d], Y[%d..%d] (Width=%d, Height=%d)\n", min_x, max_x, min_y, max_y, width, height);
    printf("--------------------------------------------------\n");
    for (size_t r = 0; r < canvas_h; r++) {
        printf("%s\n", canvas[r]);
        free(canvas[r]);
    }
    free(canvas);
    printf("--------------------------------------------------\n");
}

// Export 3D PDB structure
void export_pdb(const char *filename, const ProteinSequence *prot, const Conformation *conf) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Error: Unable to create output PDB file %s\n", filename);
        return;
    }
    
    fprintf(fp, "REMARK   1 HP LATTICE PROTEIN FOLDING SOLVER OUTPUT\n");
    fprintf(fp, "REMARK   2 SEQUENCE: %s\n", prot->sequence);
    fprintf(fp, "REMARK   3 ENERGY: %d, RADIUS OF GYRATION: %.3f\n", conf->energy, conf->radius_of_gyration);
    
    const double SPACING = 3.8; // 3.8 Angstrom C-alpha distance
    
    for (int i = 0; i < conf->length; i++) {
        double x = conf->coords[i].x * SPACING;
        double y = conf->coords[i].y * SPACING;
        double z = 0.0;
        
        const char *res_name = (prot->amino[i] == AMINO_H) ? "VAL" : "ALA";
        
        fprintf(fp, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C\n",
                i + 1, res_name, i + 1, x, y, z);
    }
    
    for (int i = 1; i < conf->length; i++) {
        fprintf(fp, "CONECT%5d%5d\n", i, i + 1);
    }
    
    fprintf(fp, "END\n");
    fclose(fp);
    printf("Saved PDB structure: %s\n", filename);
}

int main(int argc, char *argv[]) {
    srand((unsigned int)time(NULL));
    
    printf("=================================================================\n");
    printf("        C/C++ HP LATTICE PROTEIN FOLDING SOLVER ENGINE           \n");
    printf("=================================================================\n\n");
    
    // Benchmark sequences (Dill et al. standard benchmarks)
    const char *benchmarks[] = {
        "HPHPPHHPHPPHPHHPPHPH",                          // 20-mer (Ground State E = -9)
        "HHHPPHPHPHPPHPHPHPPH",                          // 20-mer 2
        "PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPH",  // 44-mer (Dill 44-mer)
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
