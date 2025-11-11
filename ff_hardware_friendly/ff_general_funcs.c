#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "ff_headers.h"


//params read
Real *read_residues(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    Real *res = malloc(sizeof(Real) * NTERMS);
    fread(res, sizeof(Real), NTERMS, f);
    fclose(f);
    printf("read residues (nterms = %d)\n", NTERMS);
    return res;
}


//multi_x read
su3_vector **read_multi_x(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    su3_vector **multi_x = malloc(sizeof(su3_vector*) * NTERMS);
    for(int t = 0; t < NTERMS; t++) {
        multi_x[t] = malloc(sizeof(su3_vector) * SITES_ON_NODE);
        fread(multi_x[t], sizeof(su3_vector), SITES_ON_NODE, f);
    }
    fclose(f);
    printf("read multi_x (nterms = %d, sites = %lu)\n",
           NTERMS, (unsigned long)SITES_ON_NODE);
    return multi_x;
}

//qpaths and num_q_paths
Q_path *read_qpaths(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    int tmp;
    fread(&tmp, sizeof(int), 1, f);
    if(tmp != NUM_Q_PATHS) {
        printf("Warning: expected %d q_paths, file has %d\n",
               NUM_Q_PATHS, tmp);
    }
    Q_path *paths = malloc(sizeof(Q_path) * NUM_Q_PATHS);
    fread(paths, sizeof(Q_path), NUM_Q_PATHS, f);
    fclose(f);
    printf("read qpaths (num_q_paths = %d)\n", NUM_Q_PATHS);
    return paths;
}

//links data
su3_matrix (*read_links(const char *fname))[4] {
    FILE *f = fopen(fname, "rb");
    if (!f) { perror(fname); exit(1); }

    size_t tmp_sites;
    fread(&tmp_sites, sizeof(size_t), 1, f);
    if (tmp_sites != SITES_ON_NODE) {
        printf("Warning: expected %lu sites, file has %lu\n",
               (unsigned long)SITES_ON_NODE, (unsigned long)tmp_sites);
    }

    su3_matrix (*links)[4] = malloc(SITES_ON_NODE * sizeof(*links));
    fread(links, sizeof(su3_matrix), SITES_ON_NODE * 4, f);
    fclose(f);

    printf("read links (sites = %lu)\n", (unsigned long)SITES_ON_NODE);
    return links;
}


/*su3_matrix **read_links(const char *fname) {
    FILE *f = fopen(fname, "rb");
    if(!f){ perror(fname); exit(1); }
    size_t tmp_sites;
    fread(&tmp_sites, sizeof(size_t), 1, f);
    if(tmp_sites != SITES_ON_NODE) {
        printf("Warning: expected %lu sites, file has %lu\n",
               (unsigned long)SITES_ON_NODE, (unsigned long)tmp_sites);
    }
    su3_matrix **links = malloc(sizeof(su3_matrix*) * SITES_ON_NODE);
    for(size_t i = 0; i < SITES_ON_NODE; i++) {
        links[i] = malloc(sizeof(su3_matrix) * 4); // 4 directions
        fread(links[i], sizeof(su3_matrix), 4, f);
    }
    fclose(f);
    printf("read links (sites = %lu)\n", (unsigned long)SITES_ON_NODE);
    return links;
}*/


//momenta read
anti_hermitmat (*read_mom(const char *fname))[4] {
    FILE *f = fopen(fname, "rb");
    if (!f) { perror(fname); exit(1); }

    size_t tmp_sites;
    if (fread(&tmp_sites, sizeof(size_t), 1, f) != 1) { perror("read_mom header"); exit(1); }
    if (tmp_sites != SITES_ON_NODE) {
        printf("Warning: expected %lu sites, file has %lu\n",
               (unsigned long)SITES_ON_NODE, (unsigned long)tmp_sites);
    }

    anti_hermitmat (*mom)[4] = malloc(SITES_ON_NODE * sizeof(*mom));
    if (!mom) { perror("malloc mom"); exit(1); }

    // File is dir-major: [dir0][all sites], [dir1][all sites], ...
    for (int dir = XUP; dir <= TUP; ++dir) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            if (fread(&mom[i][dir], sizeof(anti_hermitmat), 1, f) != 1) {
                perror("read_mom data");
                exit(1);
            }
        }
    }

    fclose(f);
    printf("read mom (sites = %lu)\n", (unsigned long)SITES_ON_NODE);
    return mom;
}

/////////////////////////////////////////
// small funcs used for calculations

//clear su3 matrix
void clear_su3mat( su3_matrix *dest ){
    int i,j;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
	        dest->e[i][j].real = dest->e[i][j].imag = 0.0;
        }
    }
    
}

//projector
void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c ){
    int i,j;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
	        CMUL_J( a->c[i], b->c[j], c->e[i][j] );
        }
    }   
}

//scalar multiply add su3 matrix
void scalar_mult_add_su3_matrix(su3_matrix *a,su3_matrix *b,Real s,
	su3_matrix *c){
    int i,j;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            c->e[i][j].real = a->e[i][j].real + s*b->e[i][j].real;
            c->e[i][j].imag = a->e[i][j].imag + s*b->e[i][j].imag;
        }
    }
}

//multiply two su3 matrices normal*normal
void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j,k;
    fcomplex x,y;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
	        x.real=x.imag=0.0;
	        for(k=0;k<3;k++){
	            CMUL( a->e[i][k] , b->e[k][j] , y );
	            CSUM( x , y );
	        }
	        c->e[i][j] = x;
        }
    }
}

//multiply two su3 matrices normal*adj
void mult_su3_na(  su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j,k;
    fcomplex x,y;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
	        x.real=x.imag=0.0;
	        for(k=0;k<3;k++){
	            CMUL_J( a->e[i][k] , b->e[j][k] , y );
	            CSUM( x , y );
	        }   
	        c->e[i][j] = x;
        }
    }
}

//multiply two su3 matrices adj*normal
void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j,k;
    fcomplex x,y;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
	        x.real=x.imag=0.0;
	        for(k=0;k<3;k++){
	            CMULJ_( a->e[k][i] , b->e[k][j], y );
	            CSUM( x , y );
	        }
	        c->e[i][j] = x;
        }
    }
}

//adjoint of a matrix
void su3_adjoint( su3_matrix *a, su3_matrix *b ){
    int i,j;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
	        CONJG( a->e[j][i], b->e[i][j] );
        }
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
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
	        CADD( a->e[i][j], b->e[i][j], c->e[i][j] );
        }
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


///////////helpers for path directions and indexing///////////////////////

//returns memory index of a given (x,y,z,t) coordinate
int site_index_from_coords(int x, int y, int z, int t)
{
    int ir = x + NX*(y + NY*(z + NZ*t));       // lexicographic, x-fastest

    int parity = (x + y + z + t) & 1;          // 0 = EVEN, 1 = ODD

    if(parity == 0)
        return ir >> 1;                        // EVEN block: ir/2
    else
        return (ir + SITES_ON_NODE) >> 1;      // ODD block: (ir+N)/2
}

 void coords_from_site_index(int idx, int *x, int *y, int *z, int *t)
{
    int neven = (int)(SITES_ON_NODE >> 1);
    int p_block = (idx < neven) ? 0 : 1;          // 0: even block, 1: odd block
    int ib = (idx < neven) ? idx : (idx - neven);       // index within the parity block

    int nX2   = NX >> 1;                          // NX/2
    int slab  = nX2;                              // per y
    int row   = nX2 * NY;                         // per z
    int plane = nX2 * NY * NZ;                    // per t

    *t = ib / plane;    ib -= (*t) * plane;
    *z = ib / row;      ib -= (*z) * row;
    *y = ib / slab;
    int x2 = ib - (*y) * slab;                    // == ib % slab

    int parity_yzt = ((*y + *z + *t) & 1);
    int add = (p_block ^ parity_yzt);             // 0 -> even x, 1 -> odd x
    *x = (x2 << 1) + add;                               // x = 2*x2 + add
}

int walk_dir(int start_idx, int dir)
{
    int x,y,z,t;
    coords_from_site_index(start_idx, &x, &y, &z, &t);

    switch(dir) {
        case XUP:   x = (x + 1) % NX; break;
        case XDOWN: x = (x - 1 + NX) % NX; break;
        case YUP:   y = (y + 1) % NY; break;
        case YDOWN: y = (y - 1 + NY) % NY; break;
        case ZUP:   z = (z + 1) % NZ; break;
        case ZDOWN: z = (z - 1 + NZ) % NZ; break;
        case TUP:   t = (t + 1) % NT; break;
        case TDOWN: t = (t - 1 + NT) % NT; break;
        
        case X3UP:    x = (x + 3) % NX; break;
        case X3DOWN:  x = (x - 3 + NX) % NX; break;
        case Y3UP:    y = (y + 3) % NY; break;
        case Y3DOWN:  y = (y - 3 + NY) % NY; break;
        case Z3UP:    z = (z + 3) % NZ; break;
        case Z3DOWN:  z = (z - 3 + NZ) % NZ; break;
        case T3UP:    t = (t + 3) % NT; break;
        case T3DOWN:  t = (t - 3 + NT) % NT; break;
    }

    return site_index_from_coords(x,y,z,t);
}


int find_backwards_gather_hw(const Q_path *path)
{
    int disp[4], i;
    /* compute total displacement of path */
    for(i=XUP;i<=TUP;i++)disp[i]=0;
    for( i=0; i<path->length; i++){
	if( GOES_FORWARDS(path->dir[i]) )
	    disp[        path->dir[i]  ]++;
	else
	    disp[OPP_DIR(path->dir[i]) ]--;
    }

   // There must be an elegant way??
   if( disp[XUP]==+1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(XDOWN);
   if( disp[XUP]==-1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(XUP);
   if( disp[XUP]== 0 && disp[YUP]==+1 && disp[ZUP]== 0 && disp[TUP]== 0 )return(YDOWN);
   if( disp[XUP]== 0 && disp[YUP]==-1 && disp[ZUP]== 0 && disp[TUP]== 0 )return(YUP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+1 && disp[TUP]== 0 )return(ZDOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-1 && disp[TUP]== 0 )return(ZUP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+1 )return(TDOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-1 )return(TUP);

   if( disp[XUP]==+3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(X3DOWN);
   if( disp[XUP]==-3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(X3UP);
   if( disp[XUP]== 0 && disp[YUP]==+3 && disp[ZUP]== 0 && disp[TUP]== 0 )return(Y3DOWN);
   if( disp[XUP]== 0 && disp[YUP]==-3 && disp[ZUP]== 0 && disp[TUP]== 0 )return(Y3UP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+3 && disp[TUP]== 0 )return(Z3DOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-3 && disp[TUP]== 0 )return(Z3UP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+3 )return(T3DOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-3 )return(T3UP);
    printf("OOOPS: NODIR\n");
   return( NODIR );
}

int sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths ){

    int netdir,dir0,dir1,dir1tmp,num_new,i,j;
    int net_back_dirs[16] = 
      { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
	X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };

    num_new=0; // number of paths in sorted table
    for( i=0; i<16; i++ ){ // loop over net_back_dirs
        netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
	for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
	  for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	    if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	    for( j=0; j<npaths; j++ ){ // pick out paths with right net displacement
	      //	    thislength = src_table[j].length;
	        if( find_backwards_gather_hw( &(src_table[j]) ) == netdir && 
			src_table[j].dir[0]==dir0 &&
			src_table[j].dir[1]==dir1tmp ){
		    dest_table[num_new] = src_table[j];
		    num_new++;
	        }
	    } // loop over paths
	  } //dir1
	} //dir0
    }
    if( num_new!=npaths){ printf("OOPS: path table error\n");}
    return 0;
} /* sort_quark_paths */


/////////////////////////////////////////////////////////////////////////


//link_transport_connection
/*void link_transport_connection(su3_matrix *src, su3_matrix *dest, su3_matrix *work, int dir, su3_matrix (*links)[4]){

    if (GOES_FORWARDS(dir)) {
        //forw: dest[i] = U_dir(i) * src[i + dir]
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            int nbr = neighbor_index_axis((int)i, dir, 1, +1);   // i + 1 along axis
            mult_su3_nn(&links[i][dir], &src[nbr], &dest[i]);
        }
    } 
    else {
        //backwards movement has two steps
        //step 1: work[i]
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            mult_su3_an(&links[i][OPP_DIR(dir)], &src[i], &work[i]);
        }
        //step 2: dest[i] = work[i - 1 along axis]  (i + (-1) step)
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            int nbr = neighbor_index_axis((int)i, OPP_DIR(dir), 1, -1);  // i - 1 along axis
            dest[i] = work[nbr];
        }
    }
}*/
