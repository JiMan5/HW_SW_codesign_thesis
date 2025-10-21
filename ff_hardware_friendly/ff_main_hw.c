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
    anti_hermitmat **mom_main = NULL;
    anti_hermitmat **mom_after_tb = NULL;

    printf("Loading binary inputs..\n\n");

    residues = read_residues("binaries/ff_residues.bin");
    links = read_links("binaries/ff_links.bin");
    multi_x = read_multi_x("binaries/ff_multi_x.bin");
    q_paths = read_qpaths("binaries/ff_qpaths.bin");
    mom_main = read_mom("binaries/ff_mom_before.bin"); //mom_before
    mom_after_tb = read_mom("binaries/ff_mom_after.bin"); //mom_after for check

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

    //call the hw_friendly function
    fermion_force_fn_multi_hw_friendly(residues, multi_x, qpaths_forward, axis, steps, sign, links, mom_main);


        // Compare mom_main (computed) vs mom_after_tb (reference)
    double max_diff = 0.0;
    double total_diff = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < SITES_ON_NODE; i++) {
        for (int dir = 0; dir < 4; dir++) {
            const anti_hermitmat *a = &mom_main[i][dir];
            const anti_hermitmat *b = &mom_after_tb[i][dir];

            // Compare all 3 complex off-diagonals
            double d01r = a->m01.real - b->m01.real;
            double d01i = a->m01.imag - b->m01.imag;
            double d02r = a->m02.real - b->m02.real;
            double d02i = a->m02.imag - b->m02.imag;
            double d12r = a->m12.real - b->m12.real;
            double d12i = a->m12.imag - b->m12.imag;

            // Compare the 3 imaginary diagonal elements
            double d00i = a->m00im - b->m00im;
            double d11i = a->m11im - b->m11im;
            double d22i = a->m22im - b->m22im;

            // Sum up squared differences
            double diff_sq =
                d01r*d01r + d01i*d01i +
                d02r*d02r + d02i*d02i +
                d12r*d12r + d12i*d12i +
                d00i*d00i + d11i*d11i + d22i*d22i;

            total_diff += diff_sq;
            if (diff_sq > max_diff) max_diff = diff_sq;
            count += 9; // 9 values per momentum matrix
        }
    }

    double rms_diff = sqrt(total_diff / (double)count);
    printf("\n=== Momentum field comparison ===\n");
    printf("Total compared elements: %zu\n", count);
    printf("RMS difference: %.6e\n", rms_diff);
    printf("Max squared difference: %.6e\n", max_diff);

    const double tol = 1e-8; // for single precision
    if (rms_diff < tol)
        printf("Result matches reference within tolerance (%.1e)\n", tol);
    else
        printf("WARNING: Difference exceeds tolerance (%.1e)\n", tol);




    free(residues);
    for(int t=0; t<NTERMS; t++) {
        free(multi_x[t]);
    }
    free(multi_x);
    free(q_paths);
    free(qpaths_forward);
    for(size_t i=0; i<SITES_ON_NODE; i++) {
        free(links[i]);
        free(mom_main[i]);
        free(mom_after_tb[i]);
    }
    free(links);
    free(mom_main);
    free(mom_after_tb);
    
    return 0;
}
