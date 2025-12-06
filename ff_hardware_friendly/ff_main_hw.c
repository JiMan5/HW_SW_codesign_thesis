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

    for(int ipath = 0; ipath<FORW_Q_PATHS; ipath++){
        if(qpaths_forward[ipath].length == 1){
            printf("path no %d ---> (", ipath);
            for(int j = 0; j < qpaths_forward[ipath].length; j++){
                printf("%d ", qpaths_forward[ipath].dir[j]);
            }
            printf(")\n");
        }
    }
    //call the hw_friendly function
    fermion_force_fn_multi_hw_friendly(netbackdir_table, residues, multi_x, qpaths_forward, links, mom_main);
    printf("Finished with the hw_friendly call!\n");

    
    int mismatch = 0;
    const float eps = 1e-5;

    for (int dir = 0; dir < 4; dir++) {
        for (size_t i = 0; i < SITES_ON_NODE; i++) {

            anti_hermitmat *A = &mom_main[i][dir];
            anti_hermitmat *B = &mom_after_tb[i][dir];

            if (fabsf(A->m01.real - B->m01.real) > eps ||
                fabsf(A->m01.imag - B->m01.imag) > eps ||
                fabsf(A->m02.real - B->m02.real) > eps ||
                fabsf(A->m02.imag - B->m02.imag) > eps ||
                fabsf(A->m12.real - B->m12.real) > eps ||
                fabsf(A->m12.imag - B->m12.imag) > eps ||
                fabsf(A->m00im    - B->m00im)    > eps ||
                fabsf(A->m11im    - B->m11im)    > eps ||
                fabsf(A->m22im    - B->m22im)    > eps)
            {
                printf("MISMATCH at dir=%d i=%lu\n",
                       dir, (unsigned long)i);

                printf("HW : (%g,%g) (%g,%g) (%g,%g)  %g %g %g\n",
                       A->m01.real, A->m01.imag,
                       A->m02.real, A->m02.imag,
                       A->m12.real, A->m12.imag,
                       A->m00im, A->m11im, A->m22im);

                printf("MILC: (%g,%g) (%g,%g) (%g,%g)  %g %g %g\n",
                       B->m01.real, B->m01.imag,
                       B->m02.real, B->m02.imag,
                       B->m12.real, B->m12.imag,
                       B->m00im, B->m11im, B->m22im);

                mismatch = 1;
            }
        }
    }

    if (mismatch)
        printf("MISMATCHES FOUND\n");
    else
        printf("TEST PASSED\n");


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
