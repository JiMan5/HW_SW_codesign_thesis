/****** fermion_force_fn_multi.c  -- ******************/
/* MIMD version 7 */
/* Multisource fermion force for an FN_type action.  Includes some
 * optimization choices.

 * External entry points in this file:

 * fermion_force_fn_multi
 * fermion_force_fn_multi_reverse
 * fermion_force_fn_multi_june05 (disused)

 * For a general FN action
 * compile with fermion_force_multi.c and fermion_force_eo_milc.c

 * For the Asqtad action
 * compile with fermion_force_asqtad.c

 * 
 * 1. General force for any FN-type action. Based on fermion_forcee_eo_milc.c
 *    and optimized to transport only one set of SU(3) matrices.
 *
 * 2. Same as 1 but with indices on CG solutions reversed
 *    to maybe improve cache hits
 *
 * Select these options with a compiler flag KS_MULTIFF and
 * compile with fermion_force_multi.c
 *
 * D.T. 12/05 Version  3. created for Asqtad RHMC.
 * D.T.  6/06 Versions 1. and 2. created.
 * CD   10/06 Collected versions in this file.
 * CD    5/07 Moved Asqtad-specific procedures to fermion_force_asqtad.c
 */


#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>   // for ptrdiff_t



/* All routines in this file require the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
//fixed consts for run
#define EPS          0.010000f
#define NTERMS       8
#define SITES_ON_NODE 131072UL
#define NUM_Q_PATHS  688
#define FORW_Q_PATHS 344
//#define MAX_PATH_LENGTH 7
#define NX 16
#define NY 16
#define NZ 16
#define NT 32

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif


/**********************************************************************/
/*   General FN Version for "nterms" sources                          */
/**********************************************************************/

/* update the  momenta with the fermion force */
/* Assumes that the multimass conjugate gradient has been run, with the answer in
   multi_x, and dslash_site(multi_x,multi_x,ODD) has been run. (fills in multi_x on odd sites) */
/* SEE LONG COMMENTS AT END */

static Q_path *q_paths_sorted = NULL;	// Quark paths sorted by net displacement and last directions
static int *netbackdir_table = NULL; // table of net path displacements (backwards from usual convention)
// table of gather directions to bring start of path to end, for "FN" = "fat-Naik" actions
static int net_back_dirs[16] = 
      { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
	X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };

static int 
sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths );

static int 
find_backwards_gather( Q_path *path );

static int first_force=1;	// 1 if force hasn't been called yet

static int ff_first_call_logged = 0;


static void write_all_mom(const char *fname) {
    FILE *f = fopen(fname, "wb");
    if(!f){ perror(fname); return; }

    fwrite(&sites_on_node, sizeof(size_t), 1, f);
    for(int dir = XUP; dir <= TUP; ++dir) {
      for(size_t i = 0; i < sites_on_node; ++i) {
          fwrite(&lattice[i].mom[dir], sizeof(anti_hermitmat), 1, f);
      }
  }
  fclose(f);
}

void dump_matrix_array_cpu(const char *fname, su3_matrix *arr) {
    FILE *f = fopen(fname, "wb");
    if (!f) { perror(fname); return; }
    fwrite(arr, sizeof(su3_matrix), sites_on_node, f);
    fclose(f);
}
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

//opposite but with correct EVEN, ODD assumptions
/*void coords_from_site_index(int idx, int *x, int *y, int *z, int *t)
{
    const int neven = SITES_ON_NODE / 2;
    int lex;

    int expected_parity = (idx < neven) ? 0 : 1;

    if (idx < neven)
        lex = idx * 2;           //even block
    else
        lex = (idx - neven) * 2 + 1;  //odd block

    //decode
    *t = lex / (NX * NY * NZ);
    lex -= (*t) * (NX * NY * NZ);

    *z = lex / (NX * NY);
    lex -= (*z) * (NX * NY);

    *y = lex / NX;
    *x = lex - (*y) * NX;

    /*if (((*x + *y + *z + *t) & 1) != expected_parity) {
        (*x)++;
        if (*x >= NX) { *x = 0; (*y)++; }
        if (*y >= NY) { *y = 0; (*z)++; }
        if (*z >= NZ) { *z = 0; (*t)++; }
    }
}

void coords_from_site_index(int idx, int *x, int *y, int *z, int *t)
{
    if(idx==10)printf("something\n");
    const int neven = SITES_ON_NODE / 2;
    int parity = (idx < neven ? 0 : 1);

    int local = (idx < neven ? idx : idx - neven);

    // How many even (or odd) sites per (t,z,y) hyperplane
    int per_hyperplane = (NX * NY * NZ) / 2;
    int per_row        = NX / 2;

    *t = local / (per_hyperplane);
    local -= (*t) * per_hyperplane;

    *z = local / (NY * per_row);
    local -= (*z) * (NY * per_row);

    *y = local / per_row;
    local -= (*y) * per_row;

    // Reconstruct x from parity and local index
    *x = local * 2 + parity;
}*/

void coords_from_site_index(int idx, int *x, int *y, int *z, int *t)
{
    int parity_group_size = (NX * NY * NZ * NT) / 2;
    int parity = (idx < parity_group_size) ? 0 : 1;
    int offset = (parity == 0) ? idx : (idx - parity_group_size);
    if(idx==1)printf("ksana edw!\n");
    // offset is now lexicographic linear index *within that parity class*
    // iterate through lexicographic (t,z,y,x) space skipping wrong parity
    int count = 0;

    for (int tt = 0; tt < NT; tt++)
        for (int zz = 0; zz < NZ; zz++)
            for (int yy = 0; yy < NY; yy++)
                for (int xx = 0; xx < NX; xx++)
                    if (((xx + yy + zz + tt) & 1) == parity) {
                        if (count == offset) {
                            *x = xx;
                            *y = yy;
                            *z = zz;
                            *t = tt;
                            return;
                        }
                        count++;
                    }
}

int walk_netbackdir(int start_idx, int netbackdir)
{
    int x,y,z,t;
    coords_from_site_index(start_idx, &x, &y, &z, &t);

    switch(netbackdir) {
        case XUP:   x = (x + 1) % NX; break;
        case XDOWN: x = (x - 1 + NX) % NX; break;
        case YUP:   y = (y + 1) % NY; break;
        case YDOWN: y = (y - 1 + NY) % NY; break;
        case ZUP:   z = (z + 1) % NZ; break;
        case ZDOWN: z = (z - 1 + NZ) % NZ; break;
        case TUP:   t = (t + 1) % NT; break;
        case TDOWN: t = (t - 1 + NT) % NT; break;
    }

    return site_index_from_coords(x,y,z,t);
}

static int milc_find_site_index_by_xyzt(int x,int y,int z,int t){
  for(size_t ii=0; ii<sites_on_node; ++ii){
    if(lattice[ii].x==x && lattice[ii].y==y && lattice[ii].z==z && lattice[ii].t==t)
      return (int)ii;
  }
  return -1; // not found (should never happen)
}

// Step coords by a single dir with periodic BC
static inline void step_xyzt(int *x,int *y,int *z,int *t,int dir){
  switch(dir){
    case XUP:   *x = (*x + 1) % NX; break;
    case XDOWN: *x = (*x - 1 + NX) % NX; break;
    case YUP:   *y = (*y + 1) % NY; break;
    case YDOWN: *y = (*y - 1 + NY) % NY; break;
    case ZUP:   *z = (*z + 1) % NZ; break;
    case ZDOWN: *z = (*z - 1 + NZ) % NZ; break;
    case TUP:   *t = (*t + 1) % NT; break;
    case TDOWN: *t = (*t - 1 + NT) % NT; break;
  }
}

void 
fermion_force_fn_multi( Real eps, Real *residues, 
			su3_vector **multi_x, int nterms, int prec,
			fermion_links_t *fl ){
  /* prec is ignored for now */
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  ks_action_paths *ap = get_action_paths(fl);
  int term;
  register int i,j,k,lastdir=-99,ipath,ilink;
  register site *s;
  int length,dir,odir;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
  su3_matrix *force_accum[4];  // accumulate force
  int netbackdir, last_netbackdir;	// backwards direction for entire path
//int tempflops = 0; //TEMP

#ifdef FFTIME
  int nflop = 966456 + 1440*nterms; // Asqtad action 11/3/06 version of code;
  double dtime;
#endif
  /* node0_printf("STARTING fermion_force_fn_multi() nterms = %d\n",nterms); */
  if( nterms==0 )return;

  /*if(!ff_first_call_logged) {
    FILE *f;

    // Scalar parameters
    f = fopen("ff_eps.bin","wb");
    fwrite(&eps, sizeof(Real), 1, f);
    fclose(f);

    f = fopen("ff_nterms.bin","wb");
    fwrite(&nterms, sizeof(int), 1, f);
    fclose(f);

    f = fopen("ff_residues.bin","wb");
    fwrite(residues, sizeof(Real), nterms, f);
    fclose(f);

    // multi_x data
    f = fopen("ff_multi_x.bin","wb");
    for(int term=0; term<nterms; term++) {
        fwrite(multi_x[term], sizeof(su3_vector), sites_on_node, f);
    }
    fclose(f);

    // Write q_paths
    ks_action_paths *ap_loc = get_action_paths(fl);
    int num_q_paths_loc = ap_loc->p.num_q_paths;
    f = fopen("ff_qpaths.bin","wb");
    fwrite(&num_q_paths_loc, sizeof(int), 1, f);
    fwrite(ap_loc->p.q_paths, sizeof(Q_path), num_q_paths_loc, f);
    fclose(f);

    // Write gauge links (inputs read during computation)
    f = fopen("ff_links.bin","wb");
    fwrite(&sites_on_node, sizeof(size_t), 1, f);
    FORALLSITES(i, s) {
        for(int dir = XUP; dir <= TUP; dir++) {
            fwrite(&(s->link[dir]), sizeof(su3_matrix), 1, f);
        }
    }
    fclose(f);

    // Momentums before
    write_all_mom("ff_mom_before.bin");

    // Meta info
    f = fopen("ff_meta.txt","w");
    fprintf(f,"sites_on_node %zu\n", sites_on_node);
    fprintf(f,"sizeof_su3_vector %zu\n", sizeof(su3_vector));
    fprintf(f,"sizeof_su3_matrix %zu\n", sizeof(su3_matrix));
    fprintf(f,"sizeof_anti_hermitmat %zu\n", sizeof(anti_hermitmat));
    fprintf(f,"nterms %d\n", nterms);
    fprintf(f,"num_q_paths %d\n", num_q_paths_loc);
    fclose(f);
  } */

  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ // 0 element is never used (it's unit matrix)
     mats_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=XUP;i<=TUP;i++){
     force_accum[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon = 2.0*eps; // we only do forward paths, gives factor of 2
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    if(netbackdir_table==NULL ) netbackdir_table = (int *)malloc( num_q_paths*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    for( ipath=0; ipath<num_q_paths; ipath++ ){
      netbackdir_table[ipath] = find_backwards_gather( &(q_paths_sorted[ipath]) );
    }
    printf("=== PATH AXIS CHECK ===\n");
    first_force=0;
  }

  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)FORALLSITES(i,s)clear_su3mat( &(force_accum[dir][i]) );

  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<num_q_paths; ipath++ ){
    this_path = &(q_paths_sorted[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash */

    length = this_path->length;
    // find gather to bring multi_x[term] from "this site" to end of path
    //netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
    netbackdir = netbackdir_table[ipath];
    // and bring multi_x to end - no gauge transformation !!
    // resulting outer product matrix has gauge transformation properties of a connection
    // from start to end of path
    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	k=0; // which gather we are using
        mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
           netbackdir, EVENANDODD, gen_pt[k] );
        FORALLSITES(i,s){
	  clear_su3mat(  &oprod_along_path[0][i] ); // actually last site in path
        }
        for(term=0;term<nterms;term++){
          if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
		sizeof(su3_vector), netbackdir, EVENANDODD, gen_pt[1-k] );
          wait_gather(mtag[k]);
          printf("\nDEBUG VECTOR COMPARE (CPU original MILC)\n");
          /*
          // ---- CHECK OUR HW INDEX MAPPING ---- //
if (ipath == 0 && term == 0) {
    int mism = 0, printed = 0, maxp = 300;

    for (int i = 0; i < sites_on_node; i++) {
        int x0 = lattice[i].x, y0 = lattice[i].y, z0 = lattice[i].z, t0 = lattice[i].t;

        // forward: coords -> idx
        int idx_fwd = site_index_from_coords(x0,y0,z0,t0);
        if (idx_fwd != i) {
            if (printed < maxp) {
                printf("FWD MISM: i=%d  idx_fwd=%d  coords=(%d,%d,%d,%d)\n", i, idx_fwd, x0,y0,z0,t0);
                printed++;
            }
            mism++;
            continue;
        }

        // inverse: idx -> coords
        int x1,y1,z1,t1;
        coords_from_site_index(i, &x1,&y1,&z1,&t1);

        if (x1!=x0 || y1!=y0 || z1!=z0 || t1!=t0) {
            if (printed < maxp) {
                printf("INV MISM: i=%d  lattice(%d,%d,%d,%d) vs mine(%d,%d,%d,%d)\n",
                       i, x0,y0,z0,t0, x1,y1,z1,t1);
                printed++;
            }
            mism++;
        }
    }

    if (mism == 0) printf("[OK] coords_from_site_index is a true inverse.\n");
    else printf("[FAIL] coords_from_site_index mismatches: %d\n", mism);
}*/
          
          if (ipath == 0 && term == 0) {  // or remove this guard to check every forward path
  int mismatches = 0;
  int printed = 0;
  const int max_print = 1000;   // print at most 5 detailed lines
  const float eps = 1e-6f;

  // netbackdir is already computed above
  FORALLSITES(i,s) {
    su3_vector *milc_vec = (su3_vector*)gen_pt[k][i];
    int real_idx = (int)(milc_vec - multi_x[term]);   // index into multi_x[term]
    int calc_idx = walk_netbackdir(i, netbackdir);

    // quick exact index check
    int idx_ok = (real_idx == calc_idx);

    // optional numeric check on vector contents
    su3_vector *hw_vec = &multi_x[term][calc_idx];
    float d0r = fabsf(milc_vec->c[0].real - hw_vec->c[0].real);
    float d0i = fabsf(milc_vec->c[0].imag - hw_vec->c[0].imag);
    float d1r = fabsf(milc_vec->c[1].real - hw_vec->c[1].real);
    float d1i = fabsf(milc_vec->c[1].imag - hw_vec->c[1].imag);
    float d2r = fabsf(milc_vec->c[2].real - hw_vec->c[2].real);
    float d2i = fabsf(milc_vec->c[2].imag - hw_vec->c[2].imag);
    int vec_ok = (d0r < eps && d0i < eps && d1r < eps && d1i < eps && d2r < eps && d2i < eps);

    if (!idx_ok || !vec_ok) {
      mismatches++;
      if (printed < max_print) {
        printf("MISMATCH: ipath=%d term=%d site=%d netbackdir=%d  real_idx=%d  calc_idx=%d\n",
               ipath, term, i, netbackdir, real_idx, calc_idx);
        printf("  milc: (%g,%g) (%g,%g) (%g,%g)\n",
               milc_vec->c[0].real, milc_vec->c[0].imag,
               milc_vec->c[1].real, milc_vec->c[1].imag,
               milc_vec->c[2].real, milc_vec->c[2].imag);
        printf("   calc: (%g,%g) (%g,%g) (%g,%g)\n",
               hw_vec->c[0].real, hw_vec->c[0].imag,
               hw_vec->c[1].real, hw_vec->c[1].imag,
               hw_vec->c[2].real, hw_vec->c[2].imag);
        printed++;
      }
    }
  }
  if (mismatches == 0)
    printf("[OK] ipath=%d term=%d: gen_pt == walk_netbackdir for all sites\n", ipath, term);
  else
    printf("[WARN] ipath=%d term=%d: mismatches=%d (showing up to %d)\n",
           ipath, term, mismatches, max_print);
}
          FORALLSITES(i,s){
	    su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
	    scalar_mult_add_su3_matrix( &oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i] );
          }
          printf("debug ended, all is good <3\n");
          cleanup_gather(mtag[k]);
	  k=1-k; // swap 0 and 1
        } /* end loop over terms in rational function expansion */

        char fname[128];
        snprintf(fname, sizeof(fname), "cpu_oprod_path_%03d.bin", ipath);
        dump_matrix_array_cpu(fname, oprod_along_path[0]);
        printf("[CPU] dumped %s\n", fname);
//tempflops+=54*nterms;
//tempflops+=36*nterms;
    }

    /* path transport the outer product, or projection matrix, of multi_x[term]
       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.

       maintain a matrix of the outer product transported backwards
	along the path to all sites on the path.
	If new "net displacement", need to completely recreate it.
	Otherwise, use as much of the previous path as possible 

	Note this array is indexed backwards - the outer product transported
	to site number n along the path is in oprod_along_path[length-n].
	This makes reusing it for later paths easier.

	Sometimes we need this at the start point of the path, and sometimes
	one link into the path, so don't always have to do the last link. */

    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( netbackdir == last_netbackdir )
      while ( j>0 && this_path->dir[j] == last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;

    for(ilink=j;ilink>=k;ilink--){
      link_transport_connection( oprod_along_path[length-ilink-1],
      oprod_along_path[length-ilink], mat_tmp0, this_path->dir[ilink]  );
//tempflops+=9*22;
    }

   /* maintain an array of transports "to this point" along the path.
	Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
    if( last_path != NULL )while( this_path->dir[ilink] == last_path->dir[ilink] ) ilink++ ;
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){
            mtag[1] = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLSITES(i,s){ su3_adjoint( (su3_matrix *)gen_pt[1][i], &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag[1]);
          }
          else{
            FORALLSITES(i,s){ mats_along_path[1][i] = s->link[OPP_DIR(dir)]; }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection( mats_along_path[ilink],
        mats_along_path[ilink+1], mat_tmp0, dir  );
//tempflops+=9*22;
      }
    } // end loop over links

    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;
      coeff = ferm_epsilon*this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;

      if(ilink==0 && GOES_FORWARDS(dir) ) FORALLSITES(i,s){
	mat_tmp0[i] = oprod_along_path[length][i]; 
      }
      else if( ilink>0) FORALLSITES(i,s){
        mult_su3_na( &(oprod_along_path[length-ilink][i]),  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
      }
//if(ilink>0)tempflops+=9*22;

      /* add in contribution to the force */
      /* Put antihermitian traceless part into momentum */
      if( ilink<length && GOES_FORWARDS(dir) ){
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[dir][i]) );
	}
	FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum[dir][i]) );
        }
//tempflops+=36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum[odir][i]) );
	}
        FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[odir][i]) );
	}
//tempflops+=36;
      }

      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted[ipath]);
  } /* end loop over paths */

  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++)FORALLSITES(i,s){
     uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
     add_su3_matrix( &tmat2, &(force_accum[dir][i]), &tmat2 );
     make_anti_hermitian( &tmat2, &(s->mom[dir]) );
  }
//tempflops+=4*18;
//tempflops+=4*18;

/*if(!ff_first_call_logged) {
  write_all_mom("ff_mom_after.bin");
  FILE *f = fopen("ff_done.txt","w");
  fprintf(f,"fermion_force_fn_multi first-call logged\n");
  fclose(f);
  ff_first_call_logged = 1;
}*/
	
  free( mat_tmp0 );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     free( mats_along_path[i] );
  }
  for(i=XUP;i<=TUP;i++){
     free( force_accum[i] );
  }
#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (FNMAT) terms = %d mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
//printf("FF flops = %d\n",tempflops);
} /* fermion_force_fn_multi */


//version with "X" vectors in "site major" order.  Since "nterms" is variable, use
// single index array
void 
fermion_force_fn_multi_reverse( Real eps, Real *residues, 
				su3_vector **multi_x, int nterms,
				fermion_links_t *fl ){
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  ks_action_paths *ap = get_action_paths(fl);
  int term;
  register int i,j,k,lastdir=-99,ipath,ilink;
  register site *s;
  int length,dir,odir;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  msg_tag *mtag;
  su3_matrix *mat_tmp0;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
  su3_matrix *force_accum[4];  // accumulate force
  int netbackdir, last_netbackdir;	// backwards direction for entire path
  su3_vector *multi_x_rev;
//int tempflops = 0; //TEMP

#ifdef FFTIME
  int nflop = 966456 + 1440*nterms; // Asqtad action 11/3/06 version of code;
  double dtime;
#endif
  /* node0_printf("STARTING fermion_force_fn_multi_reverse() nterms = %d\n",nterms); */
  if( nterms==0 )return;

  multi_x_rev = (su3_vector *)special_alloc( nterms*sites_on_node*sizeof(su3_vector) );
  if( multi_x_rev==NULL ){printf("NO ROOM\n"); exit(0);}
  FORALLSITES(i,s){
        for( j=0; j<nterms; j++){
            multi_x_rev[i*nterms+j] = multi_x[j][i];
        }
  }
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ // 0 element is never used (it's unit matrix)
     mats_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=XUP;i<=TUP;i++){
     force_accum[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  mat_tmp0 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon = 2.0*eps; // we only do forward paths, gives factor of 2
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    if(netbackdir_table==NULL ) netbackdir_table = (int *)malloc( num_q_paths*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    for( ipath=0; ipath<num_q_paths; ipath++ )
	netbackdir_table[ipath] = find_backwards_gather( &(q_paths_sorted[ipath]) );
    first_force=0;
  }

  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)FORALLSITES(i,s)clear_su3mat( &(force_accum[dir][i]) );

  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<num_q_paths; ipath++ ){
    this_path = &(q_paths_sorted[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash */

    length = this_path->length;
    // find gather to bring multi_x[term] from "this site" to end of path
    //netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
    netbackdir = netbackdir_table[ipath];
    // and bring multi_x to end - no gauge transformation !!
    // resulting outer product matrix has gauge transformation properties of a connection
    // from start to end of path
    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
        mtag = start_gather_field( multi_x_rev, nterms*sizeof(su3_vector),
           netbackdir, EVENANDODD, gen_pt[0] );
        FORALLSITES(i,s){
	  clear_su3mat(  &oprod_along_path[0][i] ); // actually last site in path
        }
        wait_gather(mtag);
        FORALLSITES(i,s){
          for(term=0;term<nterms;term++){
	    su3_projector( &(multi_x_rev[nterms*i+term]), ((su3_vector *)gen_pt[0][i])+term, &tmat );
	    scalar_mult_add_su3_matrix( &oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i] );
          } /* end loop over terms in rational function expansion */
        }
//tempflops+=54*nterms;
//tempflops+=36*nterms;
        cleanup_gather(mtag);
    }

    /* path transport the outer product, or projection matrix, of multi_x[term]
       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.

       maintain a matrix of the outer product transported backwards
	along the path to all sites on the path.
	If new "net displacement", need to completely recreate it.
	Otherwise, use as much of the previous path as possible 

	Note this array is indexed backwards - the outer product transported
	to site number n along the path is in oprod_along_path[length-n].
	This makes reusing it for later paths easier.

	Sometimes we need this at the start point of the path, and sometimes
	one link into the path, so don't always have to do the last link. */

    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( netbackdir == last_netbackdir )
      while ( j>0 && this_path->dir[j] == last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;

    for(ilink=j;ilink>=k;ilink--){
      link_transport_connection( oprod_along_path[length-ilink-1],
      oprod_along_path[length-ilink], mat_tmp0, this_path->dir[ilink]  );
//tempflops+=9*22;
    }

   /* maintain an array of transports "to this point" along the path.
	Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
    if( last_path != NULL )while( this_path->dir[ilink] == last_path->dir[ilink] ) ilink++ ;
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){
            mtag = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag);
            FORALLSITES(i,s){ su3_adjoint( (su3_matrix *)gen_pt[1][i], &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag);
          }
          else{
            FORALLSITES(i,s){ mats_along_path[1][i] = s->link[OPP_DIR(dir)]; }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection( mats_along_path[ilink],
        mats_along_path[ilink+1], mat_tmp0, dir  );
//tempflops+=9*22;
      }
    } // end loop over links

    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;
      coeff = ferm_epsilon*this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;

      if(ilink==0 && GOES_FORWARDS(dir) ) FORALLSITES(i,s){
	mat_tmp0[i] = oprod_along_path[length][i]; 
      }
      else if( ilink>0) FORALLSITES(i,s){
        mult_su3_na( &(oprod_along_path[length-ilink][i]),  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
      }
//if(ilink>0)tempflops+=9*22;

      /* add in contribution to the force */
      /* Put antihermitian traceless part into momentum */
      if( ilink<length && GOES_FORWARDS(dir) ){
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[dir][i]) );
	}
	FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum[dir][i]) );
        }
//tempflops+=36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum[odir][i]) );
	}
        FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[odir][i]) );
	}
//tempflops+=36;
      }

      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted[ipath]);
  } /* end loop over paths */

  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++)FORALLSITES(i,s){
     uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
     add_su3_matrix( &tmat2, &(force_accum[dir][i]), &tmat2 );
     make_anti_hermitian( &tmat2, &(s->mom[dir]) );
  }
//tempflops+=4*18;
//tempflops+=4*18;
	
  free( multi_x_rev );
  free( mat_tmp0 );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     free( mats_along_path[i] );
  }
  for(i=XUP;i<=TUP;i++){
     free( force_accum[i] );
  }
#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (FNMATREV) terms = %d mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
//printf("FF flops = %d\n",tempflops);
} /* fermion_force_fn_multi_reverse */


#if 0
// OLDER VERSION BEFORE SOME OPTIMIZATIONS.  THIS IS PROBABLY CLEARER AS TO
// WHAT IS GOING ON

void 
fermion_force_fn_multi_june05( Real eps, Real *residues, 
			       su3_vector **multi_x, int nterms,
			       fermion_links_t *fl ){
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  ks_action_paths *ap = get_action_paths(fl);
  int term;
  register int i,dir,lastdir=-99,ipath,ilink;
  register site *s;
  int length;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;
  Q_path *this_path;	// pointer to current path
  msg_tag *mtag1;
  su3_matrix *mat_outerprod,*mat_tmp0,*mat_tmp1,*tmp_matpt;
  int netbackdir, last_netbackdir;

#ifdef FFTIME
  int nflop = 0;
  double dtime;
#endif
  /*node0_printf("STARTING fermion_force_fn_multi() nterms = %d\n",nterms); */
  if( nterms==0 )return;

  mat_outerprod = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon = 2.0*eps;
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    first_force=0;
  }

	  /* loop over paths, and loop over links in path */
	  last_netbackdir = NODIR;
	  for( ipath=0; ipath<num_q_paths; ipath++ ){
    	    this_path = &(q_paths_sorted[ipath]);
	    if(this_path->forwback== -1)continue;	/* skip backwards dslash */

	    length = this_path->length;
	    // find gather to bring multi_x[term] from "this site" to end of path
	    netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
	    // and bring multi_x to end - no gauge transformation !!
	    // resulting outer product matrix has gauge transformation properties of a connection
	    // from start to end of path
	    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	        FORALLSITES(i,s){
		  clear_su3mat(  &mat_outerprod[i] );
	        }
	        for(term=0;term<nterms;term++){
	          mtag1 = start_gather_field( multi_x[term], sizeof(su3_vector),
	             netbackdir, EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		    su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[1][i], &tmat );
		    scalar_mult_add_su3_matrix( &mat_outerprod[i], &tmat, residues[term], &mat_outerprod[i] );
	          }
	          cleanup_gather(mtag1);
	        } /* end loop over terms in rational function expansion */
	    }
	
	    /* path transport the outer product, or projection matrix, of multi_x[term]
	       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.  Sometimes
		we need them at the start point of the path, and sometimes
		one link into the path - an optimization for later */
	    path_transport_connection( mat_outerprod, mat_tmp0, EVENANDODD, this_path->dir, length );

	    /* A path has (length+1) points, counting the ends.  At first
		 point, no "down" direction links have their momenta "at this
		 point". At last, no "up" ... */
	    for( ilink=0; ilink<=length; ilink++ ){
	      if(ilink<length)dir = this_path->dir[ilink];
	      else dir=NODIR;
	      coeff = ferm_epsilon*this_path->coeff;
	      if( (ilink%2)==1 )coeff = -coeff;
	
	      /* add in contribution to the force */
	      /* Put antihermitian traceless part into momentum */
	      if( ilink<length && GOES_FORWARDS(dir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[dir]) );
	        }
	      }
	      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[OPP_DIR(lastdir)]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[OPP_DIR(lastdir)]) );
	        }
	      }

	      /* path transport multi_x[term] and Dslash*multi_x[term] to next point */
	      /* sometimes we don't need them at last point */
	      if( ilink<length-1 || 
		(ilink==length-1 && GOES_BACKWARDS(dir)) ){
	
		if( GOES_FORWARDS(dir) ){
		  FORALLSITES(i,s){
	            mult_su3_an( &(s->link[dir]), &(mat_tmp0[i]), &(tmat) );
	            mult_su3_nn( &(tmat), &(s->link[dir]), &(mat_tmp1[i]) );
		  }
		  mtag1 = start_gather_field( mat_tmp1, sizeof(su3_matrix),
	             OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		     mat_tmp0[i] = *(su3_matrix *)gen_pt[1][i];
		  }
	          cleanup_gather(mtag1);
		}
		else{   /* GOES_BACKWARDS(dir) */
	          mtag1 = start_gather_field( mat_tmp0, sizeof(su3_matrix),
	               OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
	            mult_su3_nn( &(s->link[OPP_DIR(dir)]), (su3_matrix *)(gen_pt[1][i]), &(tmat) );
	            mult_su3_na( &(tmat), &(s->link[OPP_DIR(dir)]), &(mat_tmp1[i]) );
		  }
		  tmp_matpt = mat_tmp0; mat_tmp0 = mat_tmp1; mat_tmp1 = tmp_matpt;
	          cleanup_gather(mtag1);
		}
	      }
	      lastdir = dir;
	
	    } /* end loop over links in path */
	    last_netbackdir = netbackdir;
	  } /* end loop over paths */
	
  free( mat_outerprod ); free( mat_tmp0 ); free( mat_tmp1 );
#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (JUN05) terms = %d mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume*nterms/(1e6*dtime*numnodes()) );
#endif
}

#endif


static int 
find_backwards_gather( Q_path *path ){

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
node0_printf("OOOPS: NODIR\n"); exit(0);
   return( NODIR );
} //find_backwards_gather

// Make a new path table.  Sorted principally by total displacement of path.
// Below that, sort by direction of first link
// Below that, sort by direction of second link - note special case of one link paths

static int 
sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths ){

    int netdir,dir0,dir1,dir1tmp,num_new,i,j;

    num_new=0; // number of paths in sorted table
    for( i=0; i<16; i++ ){ // loop over net_back_dirs
        netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
	for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
	  for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	    if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	    for( j=0; j<npaths; j++ ){ // pick out paths with right net displacement
	      //	    thislength = src_table[j].length;
	        if( find_backwards_gather( &(src_table[j]) ) == netdir && 
			src_table[j].dir[0]==dir0 &&
			src_table[j].dir[1]==dir1tmp ){
		    dest_table[num_new] = src_table[j];
		    num_new++;
	        }
	    } // loop over paths
	  } //dir1
	} //dir0
    }
    if( num_new!=npaths){ node0_printf("OOPS: path table error\n"); exit(0); }
    return 0;
} /* sort_quark_paths */


/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X

 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.

Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:

 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 

	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.

The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:


	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)



	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)



	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)

As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.

Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.

Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.

*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.

Examples:
  Fat Link:
    Positive:	X-------X

    Negative     --------
	 	|	|
		|	|
		X	X

  Naik connection, smeared
    Positive:	X-------x-------x-------X

    Negative:	---------
		|	|
		|	|
		X	x-------x-------X

    Positive:	--------x--------
		|		|
		|		|
		X		x-------X

    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/



/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y

		X
		|
		|
		X
		|
		|
	X---->--X

  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/
