#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ff_headers.h"

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
su3_vector (*read_multi_x(const char *fname))[SITES_ON_NODE]
{
    FILE *f = fopen(fname, "rb");
    if (!f) { perror(fname); exit(1); }

    su3_vector (*multi_x)[SITES_ON_NODE] =
        malloc(NTERMS * sizeof(*multi_x));
    if (!multi_x) { perror("malloc multi_x"); exit(1); }

    for (int t = 0; t < NTERMS; t++) {
        if (fread(multi_x[t], sizeof(su3_vector),
                  SITES_ON_NODE, f) != SITES_ON_NODE)
        {
            perror("read_multi_x");
            exit(1);
        }
    }

    fclose(f);
    printf("read multi_x (nterms=%d, sites=%lu)\n",
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
        case TUP:   t = (t + 1) % NT_SITES; break;
        case TDOWN: t = (t - 1 + NT_SITES) % NT_SITES; break;

        case X3UP:    x = (x + 3) % NX; break;
        case X3DOWN:  x = (x - 3 + NX) % NX; break;
        case Y3UP:    y = (y + 3) % NY; break;
        case Y3DOWN:  y = (y - 3 + NY) % NY; break;
        case Z3UP:    z = (z + 3) % NZ; break;
        case Z3DOWN:  z = (z - 3 + NZ) % NZ; break;
        case T3UP:    t = (t + 3) % NT_SITES; break;
        case T3DOWN:  t = (t - 3 + NT_SITES) % NT_SITES; break;
    }

    return site_index_from_coords(x,y,z,t);
}

/*int (*create_nbr_table(void))[16]
{
    int (*nbr_table)[16] = malloc(SITES_ON_NODE * sizeof *nbr_table);
    if (!nbr_table) { perror("malloc nbr_table"); exit(1); }

    for (size_t i = 0; i < SITES_ON_NODE; i++) {
        for (int d = 0; d < 16; d++) {
            nbr_table[i][d] = walk_dir((int)i, d);
        }
    }

    printf("Neighbor table created. Sites = %lu, directions = 16\n",
           (unsigned long)SITES_ON_NODE);

    return nbr_table;
}*/

lookup_t *create_lookup_table(void)
{
    lookup_t *lookups = malloc(SITES_ON_NODE * sizeof(lookup_t));
    if (!lookups) { perror("malloc lookups"); exit(1); }

    for (size_t i = 0; i < SITES_ON_NODE; i++) {

        /* compute parity */
        int x, y, z, t;
        coords_from_site_index((int)i, &x, &y, &z, &t);
        lookups[i].parity = (x + y + z + t) & 1;

        /* compute neighbors */
        for (int d = 0; d < 16; d++) {
            lookups[i].nbr[d] = walk_dir((int)i, d);
        }
    }

    printf("lookup table created sites = %lu, struct size = %zu bytes\n", SITES_ON_NODE, sizeof(lookup_t));

    return lookups;
}


int main(void) {

    //inputs reads
    int *netbackdir_table = (int *)malloc( FORW_Q_PATHS * sizeof(int) );
    Real *residues = NULL;
    //su3_vector **multi_x = NULL;
    su3_vector (*multi_x)[SITES_ON_NODE] = NULL;
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
    int j, ipath;
    for (j = 0; j < NUM_Q_PATHS; j++) {
        if (qpaths_sorted[j].forwback == 1) {
            qpaths_forward[idx] = qpaths_sorted[j];
            netbackdir_table[idx] = find_backwards_gather_hw( &(qpaths_forward[idx]) );
            idx++;
        }
    }

    //int (*nbr_table)[16] = create_nbr_table();
    lookup_t *lookups = create_lookup_table();

    //call the hw_friendly function
    fermion_force_fn_multi_hw_friendly(netbackdir_table, residues, multi_x, qpaths_forward, links, mom_main, lookups);
    printf("Finished with the hw_friendly call!\n");


    int mismatch = 0;
    const float eps = 1e-5;
    size_t i;
    int dir;

    for (dir = 0; dir < 4; dir++) {
        for (i = 0; i < SITES_ON_NODE; i++) {

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

    int t;
    free(residues);
    free(multi_x);
    free(q_paths);
    free(qpaths_sorted);
    free(qpaths_forward);
    free(links);
    free(mom_main);
    free(mom_after_tb);

    return 0;
}
