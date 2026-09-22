/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "core/protein.h"
#include "analytics/energy.h"
#include "analytics/metrics.h"

void test_sequential_contact_zero() {
    ProteinSequence prot;
    init_protein(&prot, "HH");
    
    Conformation conf;
    conf.length = 2;
    conf.dimension = 2;
    conf.directions[0] = 1; // RIGHT
    assert(build_coords_from_directions(&conf));
    
    int energy = calculate_energy(&prot, &conf);
    // Sequential H-H contact does NOT count towards potential energy
    assert(energy == 0);
    printf("[TEST PASSED] Sequential H-H contact energy is zero (E = %d)\n", energy);
}

void test_2d_closed_loop_contact() {
    ProteinSequence prot;
    init_protein(&prot, "HHHH");
    
    Conformation conf;
    conf.length = 4;
    conf.dimension = 2;
    // 0: (0,0) -> +X -> (1,0) -> -Y -> (1,-1) -> -X -> (0,-1)
    // Residue 0 at (0,0) and residue 3 at (0,-1) are neighbors (|0-3| = 3 > 1)
    conf.directions[0] = 1; // +X (Right)
    conf.directions[1] = 2; // -Y (Down)
    conf.directions[2] = 3; // -X (Left)
    assert(build_coords_from_directions(&conf));
    
    int energy = calculate_energy(&prot, &conf);
    assert(energy == -1);
    printf("[TEST PASSED] 2D non-sequential H-H contact detected (E = %d)\n", energy);
}

void test_3d_z_contact() {
    ProteinSequence prot;
    init_protein(&prot, "HHHHHH");
    
    Conformation conf;
    conf.length = 6;
    conf.dimension = 3;
    // Walk in 3D:
    // 0: (0,0,0)
    // 1: (1,0,0) (+X, dir 1)
    // 2: (1,1,0) (+Y, dir 0)
    // 3: (0,1,0) (-X, dir 3)
    // 4: (0,1,1) (+Z, dir 4)
    // 5: (0,0,1) (-Y, dir 2)
    // Residue 5 (0,0,1) is adjacent along Z to residue 0 (0,0,0)!
    conf.directions[0] = 1; // +X
    conf.directions[1] = 0; // +Y
    conf.directions[2] = 3; // -X
    conf.directions[3] = 4; // +Z
    conf.directions[4] = 2; // -Y
    assert(build_coords_from_directions_3d(&conf));
    
    int energy = calculate_energy_3d(&prot, &conf);
    // Non-sequential contact between 0 and 5 along Z axis
    assert(energy == -1);
    printf("[TEST PASSED] 3D inter-plane Z-axis contact detected (E = %d)\n", energy);
}

void test_saw_collision() {
    Conformation conf;
    conf.length = 5;
    conf.dimension = 2;
    // Attempt square loop that collides back into origin:
    // (0,0) -> +X -> +Y -> -X -> -Y -> (0,0) [COLLISION]
    conf.directions[0] = 1; // +X
    conf.directions[1] = 0; // +Y
    conf.directions[2] = 3; // -X
    conf.directions[3] = 2; // -Y
    bool valid = build_coords_from_directions(&conf);
    assert(!valid);
    printf("[TEST PASSED] Self-avoiding walk collision correctly detected\n");
}

void test_radius_of_gyration() {
    Conformation conf;
    conf.length = 3;
    conf.dimension = 3;
    conf.coords[0].x = 0; conf.coords[0].y = 0; conf.coords[0].z = 0;
    conf.coords[1].x = 1; conf.coords[1].y = 0; conf.coords[1].z = 0;
    conf.coords[2].x = 2; conf.coords[2].y = 0; conf.coords[2].z = 0;
    double rg = calculate_radius_of_gyration(&conf);
    // cx = 1, cy = 0, cz = 0. dx = {-1, 0, 1}. sum = 2. rg = sqrt(2/3) ~ 0.8165
    assert(fabs(rg - sqrt(2.0 / 3.0)) < 1e-6);
    printf("[TEST PASSED] 3D Radius of Gyration calculation (Rg = %.4f)\n", rg);
}

int main() {
    printf("=================================================================\n");
    printf("               RUNNING HP SOLVER UNIT TEST SUITE                 \n");
    printf("=================================================================\n\n");
    test_sequential_contact_zero();
    test_2d_closed_loop_contact();
    test_3d_z_contact();
    test_saw_collision();
    test_radius_of_gyration();
    printf("\nAll 2D & 3D Unit Tests Passed Successfully!\n");
    return 0;
}
