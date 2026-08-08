/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "analytics/metrics.h"
#include <math.h>

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
