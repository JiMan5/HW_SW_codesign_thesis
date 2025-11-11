#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ff_headers.h"

int main(void) {

    //inputs reads
    int *netbackdir_table = (int *)malloc( FORW_Q_PATHS * sizeof(int) );
    Real *residues = NULL;
    su3_vector **multi_x = NULL;
    Q_path *q_paths = NULL;
    Q_path *qpaths_forward = malloc(sizeof(Q_path) * FORW_Q_PATHS); //only with forwback == 1
    Q_path *qpaths_sorted = malloc(sizeof(Q_path) * NUM_Q_PATHS);
    su3_matrix (*links)[4] = NULL;
    anti_hermitmat (*mom_main)[4] = NULL;
    anti_hermitmat (*mom_after_tb)[4] = NULL;

    printf("Loading binary inputs..\n\n");

    residues = read_residues("binaries/ff_residues.bin");
    links = read_links("binaries/ff_links.bin");
    multi_x = read_multi_x("binaries/ff_multi_x.bin");
    q_paths = read_qpaths("binaries/ff_qpaths.bin");
    mom_main = read_mom("binaries/ff_mom_before.bin"); //mom_before
    mom_after_tb = read_mom("binaries/ff_mom_after.bin"); //mom_after for check

    printf("\nAll binary data loaded\n");
    sort_quark_paths(q_paths, qpaths_sorted, NUM_Q_PATHS);

    //input only the paths with forwback == 1
    int idx = 0;
    for (int i = 0; i < NUM_Q_PATHS; i++) {
        if (qpaths_sorted[i].forwback == 1) {
            qpaths_forward[idx] = qpaths_sorted[i];
            netbackdir_table[idx] = find_backwards_gather_hw( &(qpaths_forward[idx]) );
            idx++;
        }
    }
    
    //call the hw_friendly function
    fermion_force_fn_multi_hw_friendly(netbackdir_table, residues, multi_x, qpaths_forward, links, mom_main);
    printf("Finished with the hw_friendly call!\n");
/*

        // Compare mom_main (computed) vs mom_after_tb (reference)
    double total_sq_diff = 0.0;
    double max_sq_diff = 0.0;
    size_t total_values = 0;

    for (size_t i = 0; i < SITES_ON_NODE; i++) {
        for (int dir = 0; dir < 4; dir++) {
            const anti_hermitmat *hw = &mom_main[i][dir];
            const anti_hermitmat *ref = &mom_after_tb[i][dir];

            // Compare all 3 complex off-diagonals
            double d01r = hw->m01.real - ref->m01.real;
            double d01i = hw->m01.imag - ref->m01.imag;
            double d02r = hw->m02.real - ref->m02.real;
            double d02i = hw->m02.imag - ref->m02.imag;
            double d12r = hw->m12.real - ref->m12.real;
            double d12i = hw->m12.imag - ref->m12.imag;

            // Compare the 3 imaginary diagonal elements
            double d00i = hw->m00im - ref->m00im;
            double d11i = hw->m11im - ref->m11im;
            double d22i = hw->m22im - ref->m22im;

            // Sum up squared differences
            double diff_sq =
                d01r*d01r + d01i*d01i +
                d02r*d02r + d02i*d02i +
                d12r*d12r + d12i*d12i +
                d00i*d00i + d11i*d11i + d22i*d22i;

            total_sq_diff += diff_sq;
            if (diff_sq > max_sq_diff) max_sq_diff = diff_sq;
            total_values += 9; // 9 values per momentum matrix
        }
    }

    double rms_diff = sqrt(total_sq_diff / (double)total_values);
    printf("\n=== Momentum field comparison ===\n");
    printf("Total compared elements: %zu\n", total_values);
    printf("RMS difference: %.8e\n", rms_diff);
    printf("Max squared difference: %.8e\n", max_sq_diff);

    const double tol = 1e-5; // for single precision
    if (rms_diff < tol)
        printf("Result matches reference within tolerance (%.1e)\n", tol);
    else
        printf("WARNING: Difference exceeds tolerance (%.1e)\n", tol);*/




    free(residues);
    for(int t=0; t<NTERMS; t++) {
        free(multi_x[t]);
    }
    free(multi_x);
    free(q_paths);
    free(qpaths_sorted);
    free(qpaths_forward);
    free(links);
    free(mom_main);
    free(mom_after_tb);
    
    return 0;
}
