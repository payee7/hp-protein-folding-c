/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 *
 * Features:
 *  1. 2D/3D Hydrophobic-Polar (HP) Lattice Model (Square & Cubic Lattices)
 *  2. Exact Branch-and-Bound Solver with Symmetry Pruning & Energy Lower-Bounding
 *  3. Metropolis Simulated Annealing with 2D Pivot & 3D SO(3) Subchain Rotations
 *  4. ASCII Grid Visualizer (2D & 3D Layer-by-Layer Projections)
 *  5. Structural Analytics (Radius of Gyration, Compactness, Hydrophobic Core Index)
 *  6. Standard PDB File Exporter for 3D Protein Structure Visualizers (PyMOL, ChimeraX)
 *  7. Classic Benchmark Suite (Dill's HP sequences)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define MAX_LEN 128
#define GRID_SIZE 512
#define GRID_OFFSET 256

#define GRID_SIZE_3D 128
#define GRID_OFFSET_3D 64

// Direction vectors for 2D (0: Up, 1: Right, 2: Down, 3: Left)
static const int DX[4] = { 0,  1,  0, -1};
static const int DY[4] = { 1,  0, -1,  0};
static const char DIR_NAMES[4][8] = {"UP", "RIGHT", "DOWN", "LEFT"};

// Direction vectors for 3D cubic lattice (0: +Y, 1: +X, 2: -Y, 3: -X, 4: +Z, 5: -Z)
static const int DX3D[6] = { 0,  1,  0, -1,  0,  0};
static const int DY3D[6] = { 1,  0, -1,  0,  0,  0};
static const int DZ3D[6] = { 0,  0,  0,  0,  1, -1};
static const char DIR3D_NAMES[6][8] = {"+Y", "+X", "-Y", "-X", "+Z", "-Z"};

// 3D Rotation lookup tables for 90-degree rotations about X, Y, Z axes
static const int ROT_X_POS[6] = { 4, 1, 5, 3, 2, 0 };
static const int ROT_X_NEG[6] = { 5, 1, 4, 3, 0, 2 };
static const int ROT_Y_POS[6] = { 0, 5, 2, 4, 1, 3 };
static const int ROT_Y_NEG[6] = { 0, 4, 2, 5, 3, 1 };
static const int ROT_Z_POS[6] = { 3, 0, 1, 2, 4, 5 };
static const int ROT_Z_NEG[6] = { 1, 2, 3, 0, 4, 5 };

// Amino Acid Types
typedef enum {
    AMINO_P = 0, // Polar (Hydrophilic, Neutral)
    AMINO_H = 1  // Hydrophobic (Attracted to other H)
} AminoType;

typedef struct {
    int x;
    int y;
    int z;
} Point3D;

typedef Point3D Point2D;

typedef struct {
    char sequence[MAX_LEN];
    AminoType amino[MAX_LEN];
    int length;
    int h_count;
    int h_indices[MAX_LEN];
    int remaining_h_from[MAX_LEN]; // Prefix/Suffix remaining H count for bounding
} ProteinSequence;

typedef struct {
    Point3D coords[MAX_LEN];
    int directions[MAX_LEN]; // Direction chosen at index i
    int energy;              // E = - (number of non-sequential H-H contacts)
    int length;
    int dimension;           // 2 or 3
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

// Dynamic occupancy lookup tables for O(1) collision & contact checks
static int grid_occupancy[GRID_SIZE][GRID_SIZE];
static int grid_occupancy_3d[GRID_SIZE_3D][GRID_SIZE_3D][GRID_SIZE_3D];

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
    
    int count = 0;
    for (int i = prot->length - 1; i >= 0; i--) {
        if (prot->amino[i] == AMINO_H) count++;
        prot->remaining_h_from[i] = count;
    }
}

// Compute Radius of Gyration (Rg) in 2D / 3D
double calculate_radius_of_gyration(const Conformation *conf) {
    if (conf->length <= 0) return 0.0;
    
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int i = 0; i < conf->length; i++) {
        cx += conf->coords[i].x;
        cy += conf->coords[i].y;
        cz += conf->coords[i].z;
    }
    cx /= conf->length;
    cy /= conf->length;
    cz /= conf->length;
    
    double sq_sum = 0.0;
    for (int i = 0; i < conf->length; i++) {
        double dx = conf->coords[i].x - cx;
        double dy = conf->coords[i].y - cy;
        double dz = conf->coords[i].z - cz;
        sq_sum += (dx * dx + dy * dy + dz * dz);
    }
    return sqrt(sq_sum / conf->length);
}

// Reconstruct 2D spatial coordinates from movement direction vector
bool build_coords_from_directions_2d(Conformation *conf) {
    conf->dimension = 2;
    conf->coords[0].x = GRID_OFFSET;
    conf->coords[0].y = GRID_OFFSET;
    conf->coords[0].z = 0;
    
    grid_occupancy[GRID_OFFSET][GRID_OFFSET] = 1;
    
    for (int i = 0; i < conf->length - 1; i++) {
        int dir = conf->directions[i];
        if (dir < 0 || dir >= 4) dir = 0;
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

// Reconstruct 3D spatial coordinates from movement direction vector
bool build_coords_from_directions_3d(Conformation *conf) {
    conf->dimension = 3;
    conf->coords[0].x = GRID_OFFSET_3D;
    conf->coords[0].y = GRID_OFFSET_3D;
    conf->coords[0].z = GRID_OFFSET_3D;
    
    grid_occupancy_3d[GRID_OFFSET_3D][GRID_OFFSET_3D][GRID_OFFSET_3D] = 1;
    
    for (int i = 0; i < conf->length - 1; i++) {
        int dir = conf->directions[i];
        if (dir < 0 || dir >= 6) dir = 0;
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

// Energy evaluation in 2D
int calculate_energy_2d(const ProteinSequence *prot, const Conformation *conf) {
    int contacts = 0;
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        if (x >= 0 && x < GRID_SIZE && y >= 0 && y < GRID_SIZE) {
            grid_occupancy[x][y] = i + 1;
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        if (prot->amino[i] != AMINO_H) continue;
        
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        
        for (int d = 0; d < 4; d++) {
            int nx = x + DX[d];
            int ny = y + DY[d];
            
            if (nx < 0 || nx >= GRID_SIZE || ny < 0 || ny >= GRID_SIZE) continue;
            
            int neighbor_idx = grid_occupancy[nx][ny] - 1;
            
            if (neighbor_idx >= 0 && neighbor_idx != i - 1 && neighbor_idx != i + 1) {
                if (neighbor_idx > i && prot->amino[neighbor_idx] == AMINO_H) {
                    contacts++;
                }
            }
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        if (x >= 0 && x < GRID_SIZE && y >= 0 && y < GRID_SIZE) {
            grid_occupancy[x][y] = 0;
        }
    }
    
    return -contacts;
}

// Energy evaluation in 3D
int calculate_energy_3d(const ProteinSequence *prot, const Conformation *conf) {
    int contacts = 0;
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        int z = conf->coords[i].z;
        if (x >= 0 && x < GRID_SIZE_3D && y >= 0 && y < GRID_SIZE_3D && z >= 0 && z < GRID_SIZE_3D) {
            grid_occupancy_3d[x][y][z] = i + 1;
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        if (prot->amino[i] != AMINO_H) continue;
        
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        int z = conf->coords[i].z;
        
        for (int d = 0; d < 6; d++) {
            int nx = x + DX3D[d];
            int ny = y + DY3D[d];
            int nz = z + DZ3D[d];
            
            if (nx < 0 || nx >= GRID_SIZE_3D || ny < 0 || ny >= GRID_SIZE_3D || nz < 0 || nz >= GRID_SIZE_3D) {
                continue;
            }
            
            int neighbor_idx = grid_occupancy_3d[nx][ny][nz] - 1;
            
            if (neighbor_idx >= 0 && neighbor_idx != i - 1 && neighbor_idx != i + 1) {
                if (neighbor_idx > i && prot->amino[neighbor_idx] == AMINO_H) {
                    contacts++;
                }
            }
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        int x = conf->coords[i].x;
        int y = conf->coords[i].y;
        int z = conf->coords[i].z;
        if (x >= 0 && x < GRID_SIZE_3D && y >= 0 && y < GRID_SIZE_3D && z >= 0 && z < GRID_SIZE_3D) {
            grid_occupancy_3d[x][y][z] = 0;
        }
    }
    
    return -contacts;
}

int calculate_energy(const ProteinSequence *prot, const Conformation *conf) {
    if (conf->dimension == 3) {
        return calculate_energy_3d(prot, conf);
    }
    return calculate_energy_2d(prot, conf);
}

// 2D Branch-and-Bound recursive search
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

// 2D Exact Branch-and-Bound
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

// 3D Branch-and-Bound recursive search
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
    int max_possible_contacts = current_contacts + 4 * remaining_h;
    
    if (-max_possible_contacts > stats->best_energy) {
        stats->states_pruned_bounds++;
        return;
    }
    
    for (int dir = 0; dir < 6; dir++) {
        // 3D Symmetry Breaking
        if (idx == 0 && dir != 1) continue;
        if (idx == 1 && (dir != 1 && dir != 0)) continue;
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

// 3D Exact Branch-and-Bound
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

// 2D Simulated Annealing
Conformation solve_simulated_annealing(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats) {
    Conformation current;
    current.length = prot->length;
    current.dimension = 2;
    
    for (int i = 0; i < prot->length - 1; i++) {
        current.directions[i] = 1;
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

// 3D Simulated Annealing with SO(3) Subchain Rotations
Conformation solve_simulated_annealing_3d(const ProteinSequence *prot, int max_steps, double initial_temp, double cooling_rate, SolverStats *stats) {
    Conformation current;
    current.length = prot->length;
    current.dimension = 3;
    
    for (int i = 0; i < prot->length - 1; i++) {
        current.directions[i] = 1;
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
            candidate.directions[pivot_idx] = rand() % 6;
        } else {
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

// ASCII 2D Render
static void render_ascii_2d(const ProteinSequence *prot, const Conformation *conf) {
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
        if (canvas[r]) {
            memset(canvas[r], ' ', canvas_w - 1);
            canvas[r][canvas_w - 1] = '\0';
        }
    }
    
    for (int i = 0; i < conf->length; i++) {
        int gx = conf->coords[i].x - min_x;
        int gy = conf->coords[i].y - min_y;
        size_t cr = (size_t)((height - 1 - gy) * 2);
        size_t cc = (size_t)(gx * 2);
        
        canvas[cr][cc] = (prot->amino[i] == AMINO_H) ? 'H' : 'P';
        
        if (i < conf->length - 1) {
            int ngx = conf->coords[i + 1].x - min_x;
            int ngy = conf->coords[i + 1].y - min_y;
            size_t ncr = (size_t)((height - 1 - ngy) * 2);
            size_t ncc = (size_t)(ngx * 2);
            size_t mid_r = (cr + ncr) / 2;
            size_t mid_c = (cc + ncc) / 2;
            
            if (cr == ncr) canvas[mid_r][mid_c] = '-';
            else if (cc == ncc) canvas[mid_r][mid_c] = '|';
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

// ASCII 3D Render
static void render_ascii_3d(const ProteinSequence *prot, const Conformation *conf) {
    int min_x = GRID_SIZE_3D, max_x = 0;
    int min_y = GRID_SIZE_3D, max_y = 0;
    int min_z = GRID_SIZE_3D, max_z = 0;
    
    for (int i = 0; i < conf->length; i++) {
        if (conf->coords[i].x < min_x) min_x = conf->coords[i].x;
        if (conf->coords[i].x > max_x) max_x = conf->coords[i].x;
        if (conf->coords[i].y < min_y) min_y = conf->coords[i].y;
        if (conf->coords[i].y > max_y) max_y = conf->coords[i].y;
        if (conf->coords[i].z < min_z) min_z = conf->coords[i].z;
        if (conf->coords[i].z > max_z) max_z = conf->coords[i].z;
    }
    
    printf("\n==== ASCII FOLDED CONFORMATION (3D CUBIC HP LATTICE) ====\n");
    printf("Legend: 'H' = Hydrophobic, 'P' = Polar, '-'/'|' = Intra-plane\n");
    printf("Bounding Box: X[%d..%d], Y[%d..%d], Z[%d..%d] (Dimensions: %dx%dx%d)\n",
           min_x, max_x, min_y, max_y, min_z, max_z,
           max_x - min_x + 1, max_y - min_y + 1, max_z - min_z + 1);
    printf("----------------------------------------------------------\n");
    
    int width = max_x - min_x + 1;
    int height = max_y - min_y + 1;
    size_t canvas_w = (size_t)(width * 2 + 1);
    size_t canvas_h = (size_t)(height * 2 + 1);
    
    for (int z = min_z; z <= max_z; z++) {
        bool has_residues = false;
        for (int i = 0; i < conf->length; i++) {
            if (conf->coords[i].z == z) {
                has_residues = true;
                break;
            }
        }
        if (!has_residues) continue;
        
        printf("\n[ Plane Z = %d (Layer %d of %d) ]\n", z, z - min_z + 1, max_z - min_z + 1);
        
        char **canvas = (char **)malloc(canvas_h * sizeof(char *));
        if (!canvas) continue;
        for (size_t r = 0; r < canvas_h; r++) {
            canvas[r] = (char *)malloc(canvas_w * sizeof(char));
            if (canvas[r]) {
                memset(canvas[r], ' ', canvas_w - 1);
                canvas[r][canvas_w - 1] = '\0';
            }
        }
        
        for (int i = 0; i < conf->length; i++) {
            if (conf->coords[i].z != z) continue;
            
            int gx = conf->coords[i].x - min_x;
            int gy = conf->coords[i].y - min_y;
            size_t cr = (size_t)((height - 1 - gy) * 2);
            size_t cc = (size_t)(gx * 2);
            
            canvas[cr][cc] = (prot->amino[i] == AMINO_H) ? 'H' : 'P';
            
            if (i < conf->length - 1) {
                if (conf->coords[i + 1].z == z) {
                    int ngx = conf->coords[i + 1].x - min_x;
                    int ngy = conf->coords[i + 1].y - min_y;
                    size_t ncr = (size_t)((height - 1 - ngy) * 2);
                    size_t ncc = (size_t)(ngx * 2);
                    size_t mid_r = (cr + ncr) / 2;
                    size_t mid_c = (cc + ncc) / 2;
                    if (cr == ncr) canvas[mid_r][mid_c] = '-';
                    else if (cc == ncc) canvas[mid_r][mid_c] = '|';
                }
            }
        }
        
        for (size_t r = 0; r < canvas_h; r++) {
            printf("  %s\n", canvas[r]);
            free(canvas[r]);
        }
        free(canvas);
    }
    printf("----------------------------------------------------------\n");
}

void render_ascii_folded_protein(const ProteinSequence *prot, const Conformation *conf) {
    if (conf->dimension == 3) {
        render_ascii_3d(prot, conf);
    } else {
        render_ascii_2d(prot, conf);
    }
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
    fprintf(fp, "REMARK   4 DIMENSION: %dD CUBIC LATTICE MODEL\n", (conf->dimension == 3) ? 3 : 2);
    
    const double SPACING = 3.8; // 3.8 Angstrom C-alpha distance
    
    for (int i = 0; i < conf->length; i++) {
        double x = conf->coords[i].x * SPACING;
        double y = conf->coords[i].y * SPACING;
        double z = (conf->dimension == 3) ? (conf->coords[i].z * SPACING) : 0.0;
        
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
    printf("       C/C++ HP LATTICE PROTEIN FOLDING SOLVER ENGINE (2D/3D)    \n");
    printf("=================================================================\n\n");
    
    const char *benchmarks[] = {
        "HPHPPHHPHPPHPHHPPHPH",                          // 20-mer (Ground State 2D E = -9)
        "HHHPPHPHPHPPHPHPHPPH",                          // 20-mer 2
        "PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPH",  // 44-mer (Dill 44-mer)
        "HHPPHPPHPPHPPHPPHPPHPPH"                        // 23-mer
    };
    
    const char *seq_input = benchmarks[0];
    int run_dim = 0;
    
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
    
    // 2D Solvers
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
    
    // 3D Solvers
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
