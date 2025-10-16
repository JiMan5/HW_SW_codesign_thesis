#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

int main(void) {

    //inputs reads
    Real *residues = NULL;
    su3_vector **multi_x = NULL;
    Q_path *q_paths = NULL;
    su3_matrix **links = NULL;
    anti_hermitmat **mom_before = NULL;

    printf("Loading binary inputs..\n\n");

    residues = read_residues("binaries/ff_residues.bin");
    links = read_links("binaries/ff_links.bin");
    multi_x = read_multi_x("binaries/ff_multi_x.bin");
    q_paths = read_qpaths("binaries/ff_qpaths.bin");

    mom_before = read_mom("binaries/ff_mom_before.bin");

    printf("\nAll binary data loaded\n");

    free(residues);
    for(int t=0; t<NTERMS_CONST; t++) {
        free(multi_x[t]);
    }
    free(multi_x);
    free(q_paths);
    for(size_t i=0; i<SITES_ON_NODE_CONST; i++) {
        free(links[i]);
        free(mom_before[i]);
    }
    free(links);
    free(mom_before);
    
    return 0;
}
