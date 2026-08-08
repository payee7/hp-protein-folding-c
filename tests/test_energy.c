/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include <stdio.h>
#include <assert.h>
#include "core/protein.h"
#include "analytics/energy.h"
#include "analytics/metrics.h"

void test_energy_calculation() {
    ProteinSequence prot;
    init_protein(&prot, "HH");
    
    Conformation conf;
    conf.length = 2;
    conf.directions[0] = 1; // RIGHT
    build_coords_from_directions(&conf);
    
    int energy = calculate_energy(&prot, &conf);
    // Sequential H-H contact does NOT count towards potential energy (must be non-sequential)
    assert(energy == 0);
    printf("[TEST PASSED] Sequential H-H contact energy test (E = %d)\n", energy);
}

int main() {
    printf("Running Energy Unit Tests...\n");
    test_energy_calculation();
    printf("All Unit Tests Passed Successfully!\n");
    return 0;
}
