/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>

#define MAX_LEN 128
#define GRID_SIZE 256
#define GRID_OFFSET 128

// Direction vectors for 2D (0: Up, 1: Right, 2: Down, 3: Left)
static const int DX[4] = { 0,  1,  0, -1};
static const int DY[4] = { 1,  0, -1,  0};

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
    int remaining_h_from[MAX_LEN];
} ProteinSequence;

typedef struct {
    Point2D coords[MAX_LEN];
    int directions[MAX_LEN];
    int energy;
    int length;
    double radius_of_gyration;
} Conformation;

typedef struct {
    unsigned long long states_explored;
    unsigned long long states_pruned_bounds;
    unsigned long long states_pruned_collision;
    int best_energy;
    int ground_state_count;
    double elapsed_time_sec;
} SolverStats;

#endif // TYPES_H
