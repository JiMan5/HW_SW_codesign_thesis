#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

int main(void) {

    //inputs reads
    Real eps;
    int nterms, num_q_paths;
    size_t sites_on_node, sites_check;
    Real *residues = NULL;
    su3_vector **multi_x = NULL;
    Q_path *q_paths = NULL;
    su3_matrix **links = NULL;
    anti_hermitmat **mom_before = NULL;

    printf("Loading binary inputs..\n\n");

    eps = read_eps("binaries/ff_eps.bin");
    nterms = read_nterms("binaries/ff_nterms.bin");
    residues = read_residues("binaries/ff_residues.bin", nterms);
    links = read_links("binaries/ff_links.bin", &sites_on_node);
    multi_x = read_multi_x("binaries/ff_multi_x.bin", nterms, sites_on_node);
    q_paths = read_qpaths("binaries/ff_qpaths.bin", &num_q_paths);

    mom_before = read_mom("binaries/ff_mom_before.bin", &sites_check);

    if(sites_check != sites_on_node) printf("Mismatch!\n");
    else printf("Sites on node matches!\n");

    printf("\nAll binary data loaded\n");

    free(residues);
    for(int t=0; t<nterms; t++) {
        free(multi_x[t]);
    }
    free(multi_x);
    free(q_paths);
    for(size_t i=0; i<sites_on_node; i++) {
        free(links[i]);
        free(mom_before[i]);
    }
    free(links);
    free(mom_before);
    
    return 0;
}
