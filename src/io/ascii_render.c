/*
 * ProteinGo: HP Lattice Protein Folding Prototyping Toolkit
 * Copyright (c) 2026 ProteinGo. Licensed under the MIT License.
 */

#include "io/ascii_render.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
