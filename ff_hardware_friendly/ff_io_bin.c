#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "ff_headers.h"

//params read
Real read_eps(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    Real eps;
    fread(&eps, sizeof(Real), 1, f);
    fclose(f);
    printf("read eps\n");
    return eps;
}

int read_nterms(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    int nterms;
    fread(&nterms, sizeof(int), 1, f);
    fclose(f);
    printf("read nterms\n");
    return nterms;
}

Real *read_residues(const char *fname, int nterms) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    Real *res = malloc(sizeof(Real)*nterms);
    fread(res, sizeof(Real), nterms, f);
    fclose(f);
    printf("read residues\n");
    return res;
}

//multi_x read
su3_vector **read_multi_x(const char *fname, int nterms, size_t sites_on_node) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    su3_vector **multi_x = malloc(sizeof(su3_vector*)*nterms);
    for(int t=0;t<nterms;t++){
        multi_x[t] = malloc(sizeof(su3_vector)*sites_on_node);
        fread(multi_x[t], sizeof(su3_vector), sites_on_node, f);
    }
    fclose(f);
    printf("read multi_x\n");
    return multi_x;
}

//qpaths and num_q_paths
Q_path *read_qpaths(const char *fname, int *num_q_paths) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    fread(num_q_paths, sizeof(int), 1, f);
    Q_path *paths = malloc(sizeof(Q_path)*(*num_q_paths));
    fread(paths, sizeof(Q_path), *num_q_paths, f);
    fclose(f);
    printf("read qpaths\n");
    return paths;
}

//links data
su3_matrix **read_links(const char *fname, size_t *sites_on_node) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    fread(sites_on_node, sizeof(size_t), 1, f);
    su3_matrix **links = malloc(sizeof(su3_matrix*)*(*sites_on_node));
    for(size_t i=0;i<*sites_on_node;i++){
        links[i] = malloc(sizeof(su3_matrix)*4); // 4 directions
        for(int dir=0; dir<4; dir++){
            fread(&links[i][dir], sizeof(su3_matrix), 1, f);
        }
    }
    fclose(f);
    printf("read links\n");
    return links;
}


//momenta read
anti_hermitmat **read_mom(const char *fname, size_t *sites_on_node) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    fread(sites_on_node, sizeof(size_t), 1, f);
    anti_hermitmat **mom = malloc(sizeof(anti_hermitmat*)*(*sites_on_node));
    for(size_t i=0;i<*sites_on_node;i++){
        mom[i] = malloc(sizeof(anti_hermitmat)*4); // 4 dirs
        for(int dir=0; dir<4; dir++){
            fread(&mom[i][dir], sizeof(anti_hermitmat), 1, f);
        }
    }
    fclose(f);
    printf("read moms\n");
    return mom;
}

/*
void *read_lattice(const char *fname, size_t *sites_on_node, size_t site_size) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    fread(sites_on_node, sizeof(size_t), 1, f);
    void *lat = malloc(site_size * (*sites_on_node));
    fread(lat, site_size, *sites_on_node, f);
    fclose(f);
    return lat;
}
*/