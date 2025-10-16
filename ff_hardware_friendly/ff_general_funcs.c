#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "ff_headers.h"


//params read
Real *read_residues(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    Real *res = malloc(sizeof(Real) * NTERMS_CONST);
    fread(res, sizeof(Real), NTERMS_CONST, f);
    fclose(f);
    printf("read residues (nterms = %d)\n", NTERMS_CONST);
    return res;
}


//multi_x read
su3_vector **read_multi_x(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    su3_vector **multi_x = malloc(sizeof(su3_vector*) * NTERMS_CONST);
    for(int t = 0; t < NTERMS_CONST; t++) {
        multi_x[t] = malloc(sizeof(su3_vector) * SITES_ON_NODE_CONST);
        fread(multi_x[t], sizeof(su3_vector), SITES_ON_NODE_CONST, f);
    }
    fclose(f);
    printf("read multi_x (nterms = %d, sites = %lu)\n",
           NTERMS_CONST, (unsigned long)SITES_ON_NODE_CONST);
    return multi_x;
}

//qpaths and num_q_paths
Q_path *read_qpaths(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    int tmp;
    fread(&tmp, sizeof(int), 1, f);
    if(tmp != NUM_Q_PATHS_CONST) {
        printf("Warning: expected %d q_paths, file has %d\n",
               NUM_Q_PATHS_CONST, tmp);
    }
    Q_path *paths = malloc(sizeof(Q_path) * NUM_Q_PATHS_CONST);
    fread(paths, sizeof(Q_path), NUM_Q_PATHS_CONST, f);
    fclose(f);
    printf("read qpaths (num_q_paths = %d)\n", NUM_Q_PATHS_CONST);
    return paths;
}

//links data
su3_matrix **read_links(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    size_t tmp_sites;
    fread(&tmp_sites, sizeof(size_t), 1, f);
    if(tmp_sites != SITES_ON_NODE_CONST) {
        printf("Warning: expected %lu sites, file has %lu\n",
               (unsigned long)SITES_ON_NODE_CONST, (unsigned long)tmp_sites);
    }
    su3_matrix **links = malloc(sizeof(su3_matrix*) * SITES_ON_NODE_CONST);
    for(size_t i = 0; i < SITES_ON_NODE_CONST; i++) {
        links[i] = malloc(sizeof(su3_matrix) * 4); // 4 directions
        fread(links[i], sizeof(su3_matrix), 4, f);
    }
    fclose(f);
    printf("read links (sites = %lu)\n", (unsigned long)SITES_ON_NODE_CONST);
    return links;
}


//momenta read
anti_hermitmat **read_mom(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    size_t tmp_sites;
    fread(&tmp_sites, sizeof(size_t), 1, f);
    if(tmp_sites != SITES_ON_NODE_CONST) {
        printf("Warning: expected %lu sites, file has %lu\n",
               (unsigned long)SITES_ON_NODE_CONST, (unsigned long)tmp_sites);
    }
    anti_hermitmat **mom = malloc(sizeof(anti_hermitmat*) * SITES_ON_NODE_CONST);
    for(size_t i = 0; i < SITES_ON_NODE_CONST; i++) {
        mom[i] = malloc(sizeof(anti_hermitmat) * 4); // 4 directions
        fread(mom[i], sizeof(anti_hermitmat), 4, f);
    }
    fclose(f);
    printf("read mom (sites = %lu)\n", (unsigned long)SITES_ON_NODE_CONST);
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

/////////////////////////////////////////
// small funcs used for calculations

//clear su3 matrix
void clear_su3mat( su3_matrix *dest ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    dest->e[i][j].real = dest->e[i][j].imag = 0.0;
    }
}

//projector
void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c ){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    CMUL_J( a->c[i], b->c[j], c->e[i][j] );
    }
}

//scalar multiply add su3 matrix
void scalar_mult_add_su3_matrix(su3_matrix *a,su3_matrix *b,Real s,
	su3_matrix *c){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
        c->e[i][j].real = a->e[i][j].real + s*b->e[i][j].real;
        c->e[i][j].imag = a->e[i][j].imag + s*b->e[i][j].imag;
    }
}

//multiply two su3 matrices normal*normal
void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j,k;
    fcomplex x,y;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    x.real=x.imag=0.0;
	    for(k=0;k<3;k++){
	        CMUL( a->e[i][k] , b->e[k][j] , y );
	        CSUM( x , y );
	    }
	    c->e[i][j] = x;
    }
}

//multiply two su3 matrices normal*adj
void mult_su3_na(  su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j,k;
    fcomplex x,y;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    x.real=x.imag=0.0;
	    for(k=0;k<3;k++){
	        CMUL_J( a->e[i][k] , b->e[j][k] , y );
	        CSUM( x , y );
	    }   
	    c->e[i][j] = x;
    }
}

//multiply two su3 matrices adj*normal
void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j,k;
    fcomplex x,y;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    x.real=x.imag=0.0;
	    for(k=0;k<3;k++){
	        CMULJ_( a->e[k][i] , b->e[k][j], y );
	        CSUM( x , y );
	    }
	    c->e[i][j] = x;
    }
}

//adjoint of a matrix
void su3_adjoint( su3_matrix *a, su3_matrix *b ){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    CONJG( a->e[j][i], b->e[i][j] );
    }
}

//uncompress anti hermitian su3 matrix
void uncompress_anti_hermitian( const anti_hermitmat * const mat_antihermit,
	su3_matrix *mat_su3 ) {
        Real temp1;
	mat_su3->e[0][0].imag=mat_antihermit->m00im;
	mat_su3->e[0][0].real=0.;
	mat_su3->e[1][1].imag=mat_antihermit->m11im;
	mat_su3->e[1][1].real=0.;
	mat_su3->e[2][2].imag=mat_antihermit->m22im;
	mat_su3->e[2][2].real=0.;
	mat_su3->e[0][1].imag=mat_antihermit->m01.imag;
	temp1=mat_antihermit->m01.real;
	mat_su3->e[0][1].real=temp1;
	mat_su3->e[1][0].real= -temp1;
	mat_su3->e[1][0].imag=mat_antihermit->m01.imag;
	mat_su3->e[0][2].imag=mat_antihermit->m02.imag;
	temp1=mat_antihermit->m02.real;
	mat_su3->e[0][2].real=temp1;
	mat_su3->e[2][0].real= -temp1;
	mat_su3->e[2][0].imag=mat_antihermit->m02.imag;
	mat_su3->e[1][2].imag=mat_antihermit->m12.imag;
	temp1=mat_antihermit->m12.real;
	mat_su3->e[1][2].real=temp1;
	mat_su3->e[2][1].real= -temp1;
	mat_su3->e[2][1].imag=mat_antihermit->m12.imag;
}

//adds two su3 matrices
void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c ) {
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	    CADD( a->e[i][j], b->e[i][j], c->e[i][j] );
    }
}

//make antihermitian matrix
void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3 ) {
    Real temp;
	
	temp = (m3->e[0][0].imag + m3->e[1][1].imag + m3->e[2][2].imag)*0.33333333333333333;
	ah3->m00im = m3->e[0][0].imag - temp;
	ah3->m11im = m3->e[1][1].imag - temp;
	ah3->m22im = m3->e[2][2].imag - temp;
	ah3->m01.real = (m3->e[0][1].real - m3->e[1][0].real)*0.5;
	ah3->m02.real = (m3->e[0][2].real - m3->e[2][0].real)*0.5;
	ah3->m12.real = (m3->e[1][2].real - m3->e[2][1].real)*0.5;
	ah3->m01.imag = (m3->e[0][1].imag + m3->e[1][0].imag)*0.5;
	ah3->m02.imag = (m3->e[0][2].imag + m3->e[2][0].imag)*0.5;
	ah3->m12.imag = (m3->e[1][2].imag + m3->e[2][1].imag)*0.5;

}/* make_anti_hermitian */