#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

int main(void) {

    //inputs reads
    Real *residues = NULL;
    su3_vector **multi_x = NULL;
    Q_path *q_paths = NULL;
    Q_path *qpaths_forward = malloc(sizeof(Q_path) * FORW_Q_PATHS); //only with forwback == 1 
    PathDisp disp_table[FORW_Q_PATHS];
    int *axis = malloc(sizeof(int) * FORW_Q_PATHS); //for input pass (pio hw friendly)
    int *steps = malloc(sizeof(int) * FORW_Q_PATHS); //for input pass (pio hw friendly)
    int *sign = malloc(sizeof(int) * FORW_Q_PATHS); //for input pass (pio hw friendly)
    su3_matrix **links = NULL;
    anti_hermitmat **mom_before = NULL;

    printf("Loading binary inputs..\n\n");

    residues = read_residues("binaries/ff_residues.bin");
    links = read_links("binaries/ff_links.bin");
    multi_x = read_multi_x("binaries/ff_multi_x.bin");
    q_paths = read_qpaths("binaries/ff_qpaths.bin");
    mom_before = read_mom("binaries/ff_mom_before.bin");

    printf("\nAll binary data loaded\n");

    //input only the paths with forwback == 1
    int idx = 0;
    for (int i = 0; i < NUM_Q_PATHS; i++) {
        if (q_paths[i].forwback == 1) {
            qpaths_forward[idx++] = q_paths[i];
        }
    }

    //compute net disp and check for invalid paths (should be zero)
    int bad_paths = 0;
    for (int i = 0; i < FORW_Q_PATHS; i++) {
        compute_net_disp(&qpaths_forward[i], &disp_table[i].axis, &disp_table[i].steps, &disp_table[i].sign);
        axis[i] = disp_table[i].axis;
        steps[i] = disp_table[i].steps;
        sign[i] = disp_table[i].sign;
        if (disp_table[i].axis == -1) {
            printf("Warning: Path %d has multi-axis displacement, skipping.\n", i);
            bad_paths++;
        }
    }
    printf("Total invalid paths: %d\n", bad_paths);

    fermion_force_fn_multi_hw_friendly(residues, multi_x, qpaths_forward, axis, steps, sign, links, mom_before);

    free(residues);
    for(int t=0; t<NTERMS; t++) {
        free(multi_x[t]);
    }
    free(multi_x);
    free(q_paths);
    free(qpaths_forward);
    for(size_t i=0; i<SITES_ON_NODE; i++) {
        free(links[i]);
        free(mom_before[i]);
    }
    free(links);
    free(mom_before);
    
    return 0;
}
