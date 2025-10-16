#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "ff_headers.h"


//neighbor indexing algo
static inline int neighbor_index(int i, int dir)
{
    int x = i % NX;
    int y = (i / NX) % NY;
    int z = (i / (NX * NY)) % NZ;
    int t = i / (NX * NY * NZ);

    switch(dir)
    {
        case XUP:   x = (x + 1) % NX; break;
        case XDOWN: x = (x + NX - 1) % NX; break;
        case YUP:   y = (y + 1) % NY; break;
        case YDOWN: y = (y + NY - 1) % NY; break;
        case ZUP:   z = (z + 1) % NZ; break;
        case ZDOWN: z = (z + NZ - 1) % NZ; break;
        case TUP:   t = (t + 1) % NT; break;
        case TDOWN: t = (t + NT - 1) % NT; break;
        default: break; // NODIR or invalid
    }

    return x + NX * (y + NY * (z + NZ * t));
}


void 
fermion_force_fn_multi(
  Real *residues,
	su3_vector **multi_x,
  Q_path *q_paths,
  su3_matrix (*links)[4],
	anti_hermitmat (*mom)[4]
  ){

  int term;
  int i,j,k,lastdir=-99,ipath,ilink;
  int length,dir,odir;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
  su3_matrix *force_accum[4];  // accumulate force
  int netbackdir, last_netbackdir;	// backwards direction for entire path

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
          FORALLSITES(i,s){
	    su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
	    scalar_mult_add_su3_matrix( &oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i] );
          }
          cleanup_gather(mtag[k]);
	  k=1-k; // swap 0 and 1
        } /* end loop over terms in rational function expansion */
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









//kalei ayth
// special case to transport a "connection" by one link, does both parities
void link_transport_connection( su3_matrix *src, su3_matrix *dest,
    su3_matrix *work, int dir ){
      register int i;
      register site *s;
      msg_tag *mtag0;
  
      if( GOES_FORWARDS(dir) ) {
      mtag0 = start_gather_field( src, sizeof(su3_matrix),
          dir, EVENANDODD, gen_pt[0] );
      wait_gather(mtag0);
      FORALLSITES_OMP(i,s,){
          mult_su3_nn( &(s->link[dir]), (su3_matrix *)(gen_pt[0][i]),
          &(dest[i]) );
      } END_LOOP_OMP;
      cleanup_gather(mtag0);
      }
  
      else{ /* GOES_BACKWARDS(dir) */
        FORALLSITES_OMP(i,s,){
          mult_su3_an( &(s->link[OPP_DIR(dir)]),
          &(src[i]), &(work[i]) );
        } END_LOOP_OMP;
        mtag0 = start_gather_field( work, sizeof(su3_matrix),
                    dir, EVENANDODD, gen_pt[0] );
        wait_gather(mtag0);
        FORALLSITES_OMP(i,s,){
      dest[i] = *(su3_matrix *)gen_pt[0][i];
        } END_LOOP_OMP;
        cleanup_gather(mtag0);
      }
  } /* link_transport_connection */