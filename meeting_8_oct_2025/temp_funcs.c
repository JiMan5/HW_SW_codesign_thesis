#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"

#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#endif

#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

#ifdef MILC_GLOBAL_DEBUG
#include "debug.h"
#endif /* MILC_GLOBAL_DEBUG */

/* For information */
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int i,meascount,traj_done, naik_index;
  int prompt;
  int s_iters, avs_iters, avbcorr_iters;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif
  
  initialize_machine(&argc,&argv);
  
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  
  starttime = dclock();
  
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");
  
  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
#ifdef MILC_GLOBAL_DEBUG
    global_current_time_step = 0;
#endif /* MILC_GLOBAL_DEBUG */
    
    for( traj_done=0; traj_done < warms; traj_done++ ){
      update();
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=0;		/* number of measurements 		*/
    avs_iters = avbcorr_iters = 0;
    
    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
      {
	int isite, idir;
	site *s;
	FORALLSITES(isite,s) {
	  for( idir=XUP;idir<=TUP;idir++ ) {
	    lattice[isite].on_step_Y[idir] = 0;
	    lattice[isite].on_step_W[idir] = 0;
	    lattice[isite].on_step_V[idir] = 0;
	  }
	}
      }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
      /* do the trajectories */
      STARTTIME;
      s_iters=update();
      ENDTIME("do one trajectory");
      
      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	STARTTIME;
	g_measure_ks( );
	ENDTIME("do gauge measurement");
#ifdef MILC_GLOBAL_DEBUG
#if FERM_ACTION == HISQ
        g_measure_plaq( );
#endif
#ifdef MEASURE_AND_TUNE_HISQ
        g_measure_tune( );
#endif /* MEASURE_AND_TUNE_HISQ */
#endif /* MILC_GLOBAL_DEBUG */
	
	
	/**************************************************************/
	/* Compute chiral condensate and related quantities           */
	
	/* Make fermion links if not already done */
	
	STARTTIME;
	invalidate_fermion_links(fn_links);
	restore_fermion_links_from_site(fn_links, param.prec_pbp);
	for(i = 0; i < param.num_pbp_masses; i++){
#if FERM_ACTION == HISQ
	  naik_index = param.ksp_pbp[i].naik_term_epsilon_index;
#else
	  naik_index = 0;
#endif
 	  f_meas_imp_field( param.npbp_reps, &param.qic_pbp[i], 
 			    param.ksp_pbp[i].mass, naik_index, fn_links);
	  
#ifdef D_CHEM_POT
	  Deriv_O6_field( param.npbp_reps, &param.qic_pbp[i],
			  param.ksp_pbp[i].mass, naik_index, fn_links);
#endif
	}
	ENDTIME("do pbp measurements");
	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
    }	/* end loop over trajectories */
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %e\n",
		   (double)avs_iters/meascount);
    }
    
    endtime = dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",(double)(endtime-starttime));
      printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
      printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
      
#ifdef HISQ_FORCE_FILTER_COUNTER
      printf("hisq_force_filter_counter = %d\n",hisq_force_filter_counter);
#endif
    }
    fflush(stdout);
    starttime = endtime;
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );
    }
    
    /* Destroy fermion links (created in readin() */
    
#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#endif
    fn_links = NULL;

  }
  free_lattice();
  
#ifdef HAVE_QUDA
  finalize_quda();
#endif
  
#ifdef HAVE_QPHIX
  finalize_qphix();
#endif
  
  normal_exit(0);
  return 0;
}







int update() {
    int step, iters=0;
    double startaction,endaction;
  #ifdef MILC_GLOBAL_DEBUG
    double tempaction;
  #endif
  #ifdef HMC
    Real xrandom;
  #endif
    int i,j;
    su3_vector **multi_x;
    int n_multi_x;	// number of vectors in multi_x
    su3_vector *sumvec;
    int iphi, int_alg, inaik, jphi, n;
    Real lambda, alpha, beta; // parameters in integration algorithms
    imp_ferm_links_t* fn;
  
    int_alg = INT_ALG;
    switch(int_alg){
      case INT_LEAPFROG:
        node0_printf("Leapfrog integration, steps= %d eps= %e\n",steps,epsilon);
        n_multi_x = max_rat_order;
        for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
        if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
      break;
      case INT_OMELYAN:
        lambda = 0.8;
        node0_printf("Omelyan integration, steps= %d eps= %e lambda= %e\n",steps,epsilon,lambda);
        if (steps %2 != 0 ){
          node0_printf("BONEHEAD! need even number of steps\n");
          terminate(1);
        }
        n_multi_x = max_rat_order;
        for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
        if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
      break;
      case INT_2G1F:
        alpha = 0.1; beta = 0.1;
        node0_printf("Omelyan integration, 2 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
            steps,epsilon,alpha,beta);
        if (steps %2 != 0 ){
            node0_printf("BONEHEAD! need even number of steps\n"); terminate(0);
        }
        n_multi_x = max_rat_order;
        for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
        if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
      break;
      case INT_3G1F:
        alpha = 0.1; beta = 0.1;
        node0_printf("Omelyan integration, 3 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
            steps,epsilon,alpha,beta);
        if (steps %2 != 0 ){
            node0_printf("BONEHEAD! need even number of steps\n"); terminate(0);
        }
        n_multi_x = max_rat_order;
        for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
        if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
      break;
      case INT_5G1F:
        alpha = 0.1; beta = 0.1;
        node0_printf("Omelyan integration, 5 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
            steps,epsilon,alpha,beta);
        if (steps %2 != 0 ){
            node0_printf("BONEHEAD! need even number of steps\n"); terminate(0);
        }
        n_multi_x = max_rat_order;
        for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
        if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
      break;
      case INT_6G1F:
        alpha = 0.1; beta = 0.1;
        node0_printf("Omelyan integration, 6 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
            steps,epsilon,alpha,beta);
        if (steps %2 != 0 ){
            node0_printf("BONEHEAD! need even number of steps\n"); terminate(0);
        }
        n_multi_x = max_rat_order;
        for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
        if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
      break;
      case INT_2EPS_3TO1:
        lambda = 0.8;
        node0_printf("Omelyan integration, steps= %d eps= %e lambda= %e\n",steps,epsilon,lambda);
        node0_printf("3X step size leapfrog for factor 1 in determinant\n");
        if ( n_pseudo != 2 ){
            node0_printf("BONEHEAD! need exactly two pseudofermions\n"); terminate(0);
        }
        if (steps %6 != 0 ){
            node0_printf("BONEHEAD! need 6N number of steps\n"); terminate(0);
        }
        n_multi_x = max_rat_order;
      break;
      default:
        n_multi_x = 0;  /* Humor the compiler error checker */
        node0_printf("No integration algorithm, or unknown one\n");
        terminate(1);
      break;
    }
    
    /* allocate space for multimass solution vectors */
  
    multi_x = (su3_vector **)malloc(n_multi_x*sizeof(su3_vector *));
    if(multi_x == NULL){
      printf("update: No room for multi_x\n");
      terminate(1);
    }
    for(i=0;i<n_multi_x;i++){
      multi_x[i]=(su3_vector *)special_alloc( sizeof(su3_vector)*sites_on_node );
      if(multi_x[i] == NULL){
        printf("update: No room for multi_x\n");
        terminate(1);
      }
    }
  
    sumvec = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
    if( sumvec==NULL ){
      printf("update: No room for sumvec\n"); 
        terminate(1);
    }
    
    /* refresh the momenta */
    ranmom();
    
    /* generate a pseudofermion configuration only at start*/
    // NOTE used to clear xxx here.  May want to clear all solutions for reversibility
    iphi=0;
  #if FERM_ACTION == HISQ
    n = fermion_links_get_n_naiks(fn_links);
  #else
    n = 1;
  #endif
    for( inaik=0; inaik<n; inaik++ ) {
      restore_fermion_links_from_site(fn_links, prec_gr[0]);
      fn = get_fm_links(fn_links, inaik);
      for( jphi=0; jphi<n_pseudo_naik[inaik]; jphi++ ) {
        grsource_imp_rhmc( F_OFFSET(phi[iphi]), &(rparam[iphi].GR), EVEN,
               multi_x, sumvec, rsqmin_gr[iphi], niter_gr[iphi],
               prec_gr[iphi], fn, inaik, 
               rparam[iphi].naik_term_epsilon);
        iphi++;
      }
      destroy_fn_links(fn);
    }
    
    /* find action */
    startaction=d_action_rhmc(multi_x,sumvec);
  #ifdef HMC
    /* copy link field to old_link */
    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
  #endif
    
    switch(int_alg){
      case INT_LEAPFROG:
        /* do "steps" microcanonical steps"  */
        for(step=1; step <= steps; step++){
          /* update U's to middle of interval */
          update_u(0.5*epsilon);
          /* now update H by full time interval */
          iters += update_h_rhmc( epsilon, multi_x);
          /* update U's by half time step to get to even time */
          update_u(epsilon*0.5);
          /* reunitarize the gauge field */
          reunitarize_ks();
          /*TEMP - monitor action*/if(step%4==0)d_action_rhmc(multi_x,sumvec);
        }	/* end loop over microcanonical steps */
      break;
      case INT_OMELYAN:
        /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
        for(step=2; step <= steps; step+=2){
          /* update U's and H's - see header comment */
          update_u(0.5*epsilon*lambda);
          iters += update_h_rhmc( epsilon, multi_x);
          update_u(epsilon*(2.0-lambda));
          iters += update_h_rhmc( epsilon, multi_x);
          update_u(0.5*epsilon*lambda);
          /* reunitarize the gauge field */
          reunitarize_ks();
          /*TEMP - monitor action*/ //if(step%4==0)d_action_rhmc(multi_x,sumvec);
        }	/* end loop over microcanonical steps */
      break;
      case INT_2G1F:
          /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
          for(step=2; step <= steps; step+=2){
          /* update U's and H's - see header comment */
               update_u( epsilon*( (0.25-0.5*alpha) ) );
          update_h_gauge( 0.5*epsilon);
               update_u( epsilon*( (0.5-beta)-(0.25-0.5*alpha) ) );
          iters += update_h_fermion( epsilon, multi_x);
               update_u( epsilon*( (0.75+0.5*alpha)-(0.5-beta) ) );
          update_h_gauge( 0.5*epsilon);
  
               update_u( epsilon*( (1.25-0.5*alpha)-(0.75+0.5*alpha) ) );
          update_h_gauge( 0.5*epsilon);
               update_u( epsilon*( (1.5+beta)-(1.25-0.5*alpha) ) );
          iters += update_h_fermion( epsilon, multi_x);
               update_u( epsilon*( (1.75+0.5*alpha)-(1.5+beta) ) );
          update_h_gauge( 0.5*epsilon);
               update_u( epsilon*( (2.0)-(1.75+0.5*alpha) ) );
  
              /* reunitarize the gauge field */
              reunitarize_ks();
              /*TEMP - monitor action*/ //if(step%6==0)d_action_rhmc(multi_x,sumvec);
          }	/* end loop over microcanonical steps */
      break;
      case INT_3G1F:
          /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
          for(step=2; step <= steps; step+=2){
  #ifdef MILC_GLOBAL_DEBUG
              global_current_time_step = step-1;
              node0_printf( "Current time step: %d\n", global_current_time_step );
  #endif /* MILC_GLOBAL_DEBUG */
          /* update U's and H's - see header comment */
              update_u( epsilon*( (1.0/6.0-alpha/3.0) ) );
          update_h_gauge( epsilon/3.0);
               update_u( epsilon*( (0.5-beta)-(1.0/6.0-alpha/3.0) ) );
          iters += update_h_fermion( epsilon, multi_x);
               update_u( epsilon*( (3.0/6.0+alpha/3.0)-(0.5-beta) ) );
          update_h_gauge( epsilon/3.0);
               update_u( epsilon*( (5.0/6.0-alpha/3.0)-(3.0/6.0+alpha/3.0) ) );
          update_h_gauge( epsilon/3.0);
               update_u( epsilon*( (7.0/6.0+alpha/3.0)-(5.0/6.0-alpha/3.0) ) );
          update_h_gauge( epsilon/3.0);
  
  #ifdef MILC_GLOBAL_DEBUG
  #ifdef HISQ_REUNITARIZATION_DEBUG
              {
              double max_delta_phase = 0.0;
              double max_delta_norm = 0.0;
              site *s;
              int idir;
              double delta_phase,delta_norm;
              su3_matrix Wdiff;
              double min_det_V, max_det_V;
              double min_eigen,max_eigen,min_denom;
              double Wm1unit; /* norm of (W^+W - I) */
              double flag_detV=1;
              FORALLSITES(i,s) {
                for( idir=XUP;idir<=TUP;idir++ ) {
                  /* phase deviation */
                  delta_phase = fabs( lattice[i].phase_Y[idir] -
                                      lattice[i].phase_Y_previous[idir] );
                  if( delta_phase > max_delta_phase ) max_delta_phase = delta_phase;
                  /* Wlink norm deviation */
                  sub_su3_matrix( &(lattice[i].Wlink[idir]),
                                  &(lattice[i].Wlink_previous[idir]), &Wdiff );
                  delta_norm = su3_norm_frob( &Wdiff );
                  if( delta_norm > max_delta_norm ) max_delta_norm = delta_norm;
  
                  if( 1==flag_detV ) {
                    min_det_V = lattice[i].Vdet[idir];
                    max_det_V = min_det_V;
                    min_eigen = lattice[i].gmin[idir];
                    max_eigen = lattice[i].gmax[idir];
                    min_denom = lattice[i].denom[idir];
                    Wm1unit = lattice[i].unitW1[idir];
                    flag_detV = 0;
                  }
                  else {
                    if( min_det_V > lattice[i].Vdet[idir] )
                      min_det_V = lattice[i].Vdet[idir];
                    if( max_det_V < lattice[i].Vdet[idir] )
                      max_det_V = lattice[i].Vdet[idir];
                    if( min_eigen > lattice[i].gmin[idir] )
                      min_eigen = lattice[i].gmin[idir];
                    if( max_eigen < lattice[i].gmax[idir] )
                      max_eigen = lattice[i].gmax[idir];
                    if( min_denom > lattice[i].denom[idir] )
                      min_denom = lattice[i].denom[idir];
                    if( Wm1unit < lattice[i].unitW1[idir] )
                      Wm1unit = lattice[i].unitW1[idir];
                  }
                }
              }
              g_doublemax( &max_delta_phase );
              g_doublemax( &max_delta_norm );
              min_det_V = -min_det_V;
              g_doublemax( &min_det_V );
              g_doublemax( &max_det_V );
              min_det_V = -min_det_V;
              min_eigen = -min_eigen;
              g_doublemax( &min_eigen );
              g_doublemax( &max_eigen );
              min_eigen = -min_eigen;
              min_denom = -min_denom;
              g_doublemax( &min_denom );
              min_denom = -min_denom;
              node0_printf("PHASE_Y maximum jump: %28.14g\n", max_delta_phase );
              node0_printf("NORM_W maximum jump: %28.14g\n", max_delta_norm );
              node0_printf("DET_V minimum: %28.14g\n", min_det_V);
              node0_printf("DET_V maximum: %28.14g\n", max_det_V);
              node0_printf("(V^+V) eigenvalue minimum: %28.14g\n", min_eigen);
              node0_printf("(V^+V) eigenvalue maximum: %28.14g\n", max_eigen);
              node0_printf("denom=ws*(us*vs-ws)  minimum: %28.14g\n", min_denom);
              node0_printf("Deviation from unitary |W^+W-1|  maximum: %28.14g\n", Wm1unit);
              }
  #endif /* HISQ_REUNITARIZATION_DEBUG */
              global_current_time_step = step;
              node0_printf( "Current time step: %d\n", global_current_time_step );
  #endif /* MILC_GLOBAL_DEBUG */
  
               update_u( epsilon*( (9.0/6.0-alpha/3.0)-(7.0/6.0+alpha/3.0) ) );
          update_h_gauge( epsilon/3.0);
               update_u( epsilon*( (1.5+beta)-(9.0/6.0-alpha/3.0) ) );
          iters += update_h_fermion( epsilon, multi_x);
               update_u( epsilon*( (11.0/6.0+alpha/3.0)-(1.5+beta) ) );
          update_h_gauge( epsilon/3.0);
               update_u( epsilon*( (2.0)-(11.0/6.0+alpha/3.0) ) );
  
              /* reunitarize the gauge field */
              reunitarize_ks();
  #ifdef MILC_GLOBAL_DEBUG
  #ifdef HISQ_REUNITARIZATION_DEBUG
              {
              double max_delta_phase = 0.0;
              double max_delta_norm = 0.0;
              site *s;
              int idir;
              double delta_phase,delta_norm;
              su3_matrix Wdiff;
              double min_det_V, max_det_V;
              double min_eigen,max_eigen,min_denom;
              double Wm1unit; /* norm of (W^+W - I) */
              double flag_detV=1;
              FORALLSITES(i,s) {
                for( idir=XUP;idir<=TUP;idir++ ) {
                  /* phase deviation */
                  delta_phase = fabs( lattice[i].phase_Y[idir] -
                                      lattice[i].phase_Y_previous[idir] );
                  if( delta_phase > max_delta_phase ) max_delta_phase = delta_phase;
                  /* Wlink norm deviation */
                  sub_su3_matrix( &(lattice[i].Wlink[idir]),
                                  &(lattice[i].Wlink_previous[idir]), &Wdiff );
                  delta_norm = su3_norm_frob( &Wdiff );
                  if( delta_norm > max_delta_norm ) max_delta_norm = delta_norm;
  
                  if( 1==flag_detV ) {
                    min_det_V = lattice[i].Vdet[idir];
                    max_det_V = min_det_V;
                    min_eigen = lattice[i].gmin[idir];
                    max_eigen = lattice[i].gmax[idir];
                    min_denom = lattice[i].denom[idir];
                    Wm1unit = lattice[i].unitW1[idir];
                    flag_detV = 0;
                  }
                  else {
                    if( min_det_V > lattice[i].Vdet[idir] )
                      min_det_V = lattice[i].Vdet[idir];
                    if( max_det_V < lattice[i].Vdet[idir] )
                      max_det_V = lattice[i].Vdet[idir];
                    if( min_eigen > lattice[i].gmin[idir] )
                      min_eigen = lattice[i].gmin[idir];
                    if( max_eigen < lattice[i].gmax[idir] )
                      max_eigen = lattice[i].gmax[idir];
                    if( min_denom > lattice[i].denom[idir] )
                      min_denom = lattice[i].denom[idir];
                    if( Wm1unit < lattice[i].unitW1[idir] )
                      Wm1unit = lattice[i].unitW1[idir];
                  }
                }
              }
              g_doublemax( &max_delta_phase );
              g_doublemax( &max_delta_norm );
              min_det_V = -min_det_V;
              g_doublemax( &min_det_V );
              g_doublemax( &max_det_V );
              min_det_V = -min_det_V;
              min_eigen = -min_eigen;
              g_doublemax( &min_eigen );
              g_doublemax( &max_eigen );
              min_eigen = -min_eigen;
              min_denom = -min_denom;
              g_doublemax( &min_denom );
              min_denom = -min_denom;
              node0_printf("PHASE_Y maximum jump: %28.14g\n", max_delta_phase );
              node0_printf("NORM_W maximum jump: %28.14g\n", max_delta_norm );
              node0_printf("DET_V minimum: %28.14g\n", min_det_V);
              node0_printf("DET_V maximum: %28.14g\n", max_det_V);
              node0_printf("(V^+V) eigenvalue minimum: %28.14g\n", min_eigen);
              node0_printf("(V^+V) eigenvalue maximum: %28.14g\n", max_eigen);
              node0_printf("denom=ws*(us*vs-ws)  minimum: %28.14g\n", min_denom);
              node0_printf("Deviation from unitary |W^+W-1|  maximum: %28.14g\n", Wm1unit);
              }
  #endif /* HISQ_REUNITARIZATION_DEBUG */
              /*TEMP - monitor action*/ tempaction = d_action_rhmc(multi_x,sumvec);
              /*TEMP - monitor action*/ node0_printf("DELTA_SO_FAR %e\n",tempaction-startaction);
  #endif /* MILC_GLOBAL_DEBUG */
          }	/* end loop over microcanonical steps */
      break;
      case INT_5G1F:
          /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
          for(step=2; step <= steps; step+=2){
          /* update U's and H's - see header comment */
          /* do 2 gauge updates */
               update_u( epsilon*( (1.0/10.0-alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
              update_u( epsilon*( (3.0/10.0+alpha/5.0)-(1.0/10.0-alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
          /* do 1 ferm update */
               update_u( epsilon*( (0.5-beta)-(3.0/10.0+alpha/5.0) ) );
          iters += update_h_fermion( epsilon, multi_x);
          /* do 3 gauge updates */
               update_u( epsilon*( (5.0/10.0-alpha/5.0)-(0.5-beta) ) );
          update_h_gauge( epsilon/5.0);
               update_u( epsilon*( (7.0/10.0+alpha/5.0)-(5.0/10.0-alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
               update_u( epsilon*( (9.0/10.0-alpha/5.0)-(7.0/10.0+alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
          /* next step */
          /* do 3 gauge updates */
               update_u( epsilon*( (11.0/10.0+alpha/5.0)-(9.0/10.0-alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
               update_u( epsilon*( (13.0/10.0-alpha/5.0)-(11.0/10.0+alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
               update_u( epsilon*( (15.0/10.0+alpha/5.0)-(13.0/10.0-alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
          /* do 1 ferm update */
               update_u( epsilon*( (1.5+beta)-(15.0/10.0+alpha/5.0) ) );
          iters += update_h_fermion( epsilon, multi_x);
          /* do 2 gauge updates */
               update_u( epsilon*( (17.0/10.0-alpha/5.0)-(1.5+beta) ) );
          update_h_gauge( epsilon/5.0);
               update_u( epsilon*( (19.0/10.0+alpha/5.0)-(17.0/10.0-alpha/5.0) ) );
          update_h_gauge( epsilon/5.0);
          /* finish */
               update_u( epsilon*( (2.0)-(19.0/10.0+alpha/5.0) ) );
          /* reunitarize the gauge field */
              reunitarize_ks();
          }
      break;
      case INT_6G1F:
          /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
          for(step=2; step <= steps; step+=2){
          /* update U's and H's - see header comment */
          /* do 3 gauge updates */
               update_u( epsilon*( (1.0/12.0-alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
              update_u( epsilon*( (3.0/12.0+alpha/6.0)-(1.0/12.0-alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (5.0/12.0-alpha/6.0)-(3.0/12.0+alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
          /* do 1 ferm update */
               update_u( epsilon*( (0.5-beta)-(5.0/12.0-alpha/6.0) ) );
          iters += update_h_fermion( epsilon, multi_x);
          /* do 3 gauge updates */
               update_u( epsilon*( (7.0/12.0+alpha/6.0)-(0.5-beta) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (9.0/12.0-alpha/6.0)-(7.0/12.0+alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (11.0/12.0+alpha/6.0)-(9.0/12.0-alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
          /* next step */
          /* do 3 gauge updates */
               update_u( epsilon*( (13.0/12.0-alpha/6.0)-(11.0/12.0+alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (15.0/12.0+alpha/6.0)-(13.0/12.0-alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (17.0/12.0-alpha/6.0)-(15.0/12.0+alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
          /* do 1 ferm update */
               update_u( epsilon*( (1.5+beta)-(17.0/12.0-alpha/6.0) ) );
          iters += update_h_fermion( epsilon, multi_x);
          /* do 3 gauge updates */
               update_u( epsilon*( (19.0/12.0+alpha/6.0)-(1.5+beta) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (21.0/12.0-alpha/6.0)-(19.0/12.0+alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
               update_u( epsilon*( (23.0/12.0+alpha/6.0)-(21.0/12.0-alpha/6.0) ) );
          update_h_gauge( epsilon/6.0);
          /* finish */
               update_u( epsilon*( (2.0)-(23.0/12.0+alpha/6.0) ) );
          /* reunitarize the gauge field */
              reunitarize_ks();
          }
      break;
      case INT_2EPS_3TO1:
  #if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
        printf("update(%d): INT_2EPS_3TO1 is not supported for HISQ or HYPISQ\n",
           this_node);
        terminate(1);
  #endif
          /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
          for(step=6; step <= steps; step+=6){
          /* update U's and H's - first Omelyan step */
               update_u(0.5*epsilon*lambda);
          imp_gauge_force_ks(epsilon,F_OFFSET(mom));
              eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
                     multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
                     niter_md[1], prec_md[1], prec_ff,
                     fn_links );
               update_u(epsilon*( 1.0 + 0.5*(1-lambda) )); // to time = (3/2)*epsilon
              eo_fermion_force_rhmc( 3.0*epsilon,  &rparam[0].MD,
                     multi_x, F_OFFSET(phi[0]), rsqmin_md[0], 
                     niter_md[0], prec_md[0], prec_ff,
                     fn_links );
  
               update_u(epsilon*( 0.5*(1.0-lambda) ));
  
          imp_gauge_force_ks(epsilon,F_OFFSET(mom));
              eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
                     multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
                     niter_md[1], prec_md[1], prec_ff,
                     fn_links );
  
              update_u(0.5*epsilon*lambda);
  
          /* update U's and H's - second Omelyan step */
               update_u(0.5*epsilon*lambda);
  
          imp_gauge_force_ks(epsilon,F_OFFSET(mom));
              eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
                     multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
                     niter_md[1], prec_md[1], prec_ff,
                     fn_links );
  
               update_u(epsilon*(2.0-lambda));
  
          imp_gauge_force_ks(epsilon,F_OFFSET(mom));
              eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
                     multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
                     niter_md[1], prec_md[1], prec_ff,
                     fn_links );
  
              update_u(0.5*epsilon*lambda);
  
          /* update U's and H's - third Omelyan step */
               update_u(0.5*epsilon*lambda);
  
          imp_gauge_force_ks(epsilon,F_OFFSET(mom));
              eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
                     multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
                     niter_md[1], prec_md[1], prec_ff,
                     fn_links );
  
               update_u(epsilon*(  0.5*(1.0-lambda) )); // to time 2*epsilon + epsilon/2
  
              eo_fermion_force_rhmc( 3.0*epsilon,  &rparam[0].MD,
                     multi_x, F_OFFSET(phi[0]), rsqmin_md[0], 
                     niter_md[0], prec_md[0], prec_ff,
                     fn_links );
  
               update_u(epsilon*( 1.0 + 0.5*(1.0-lambda) ));
  
          imp_gauge_force_ks(epsilon,F_OFFSET(mom));
              eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
                     multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
                     niter_md[1], prec_md[1], prec_ff,
                     fn_links );
  
              update_u(0.5*epsilon*lambda);
  
              /* reunitarize the gauge field */
          reunitarize_ks();
              /*TEMP - monitor action*/ //if(step%6==0)d_action_rhmc(multi_x,sumvec);
  
          }	/* end loop over microcanonical steps */
      break;
      default:
        node0_printf("No integration algorithm, or unknown one\n");
        terminate(1);
      break;
    }
  
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    endaction=d_action_rhmc(multi_x,sumvec);
    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
  #ifdef HMC
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( (double)(startaction-endaction) ) < xrandom ){
      if(steps > 0)
        gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
  #ifdef FN
      invalidate_fermion_links(fn_links);
  #endif
      node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
    }
    else {
      node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
    }
  #else  // not HMC
    node0_printf("CHECK: delta S = %e\n", (double)(endaction-startaction));
  #endif // HMC
    
    /* free multimass solution vector storage */
    for(i=0;i<n_multi_x;i++)special_free(multi_x[i]);
    free(sumvec);
    
    if(steps > 0)return (iters/steps);
    else return(-99);
  }










#include "ks_imp_includes.h"	/* definitions files and prototypes */

int update_h_rhmc( Real eps, su3_vector **multi_x ){
  int iters;
#ifdef FN
  invalidate_fermion_links(fn_links);
#endif
  /*  node0_printf("update_h_rhmc:\n"); */
  /* gauge field force */
  imp_gauge_force_ks(eps,F_OFFSET(mom));
  /* fermionic force */
  
  iters = update_h_fermion( eps,  multi_x );
  return iters;
} /* update_h_rhmc */

// gauge and fermion force parts separately, for algorithms that use
// different time steps for them
void update_h_gauge( Real eps ){
  /* node0_printf("update_h_gauge:\n");*/
  /* gauge field force */
  imp_gauge_force_ks(eps,F_OFFSET(mom));
} /* update_h_gauge */

// fermion force update grouping pseudofermions with the same path coeffs
int update_h_fermion( Real eps, su3_vector **multi_x ){
  int iphi,jphi;
  Real final_rsq;
  int i,j,n;
  int order, tmporder;
  Real *residues,*allresidues;
  Real *roots;
  int iters = 0;
  imp_ferm_links_t *fn;

  /* Algorithm sketch: assemble multi_x with all |X> fields,
     then call force routine for each part (so far we have to parts:
     zero correction to Naik and non-zero correction to Naik */

  allresidues = (Real *)malloc(n_order_naik_total*sizeof(Real));

  // Group the fermion force calculation according to sets of like
  // path coefficients.
  tmporder = 0;
  iphi = 0;
  restore_fermion_links_from_site(fn_links, prec_md[0]);
#if FERM_ACTION == HISQ
  n = fermion_links_get_n_naiks(fn_links);
#else
  n = 1;
#endif
  for( i=0; i<n; i++ ) {
    /* Assume prec_md is the same for all pseudo-fermions
       with the same Naik epsilon */
    fn = get_fm_links(fn_links, i);
    for( jphi=0; jphi<n_pseudo_naik[i]; jphi++ ) {
      
      // Add the current pseudofermion to the current set
      order = rparam[iphi].MD.order;
      residues = rparam[iphi].MD.res;
      roots = rparam[iphi].MD.pole;

      // Compute ( M^\dagger M)^{-1} in xxx_even
      // Then compute M*xxx in temporary vector xxx_odd 
      /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
      iters += ks_ratinv( F_OFFSET(phi[iphi]), multi_x+tmporder, roots, residues,
                          order, niter_md[iphi], rsqmin_md[iphi], prec_md[iphi], EVEN,
			  &final_rsq, fn, i, rparam[iphi].naik_term_epsilon );

      for(j=0;j<order;j++){
	dslash_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD, fn);
	allresidues[tmporder+j] = residues[j+1];
	// remember that residues[0] is constant, no force contribution.
      }
      tmporder += order;
      iphi++;
    }
    destroy_fn_links(fn);
  }

#ifdef MILC_GLOBAL_DEBUG
  node0_printf("update_h_rhmc: MULTI_X ASSEMBLED\n");fflush(stdout);
  node0_printf("update_h_rhmc: n_distinct_Naik=%d\n",n);
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: orders[%d]=%d\n",j,n_orders_naik[j]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: masses_Naik[%d]=%f\n",j,fn_links.hl.eps_naik[j]);
#endif
  fflush(stdout);
#endif /* MILC_GLOBAL_DEBUG */

  restore_fermion_links_from_site(fn_links, prec_ff);
  eo_fermion_force_multi( eps, allresidues, multi_x,
			  n_order_naik_total, prec_ff, fn_links );
  free(allresidues);
  return iters;
} /* update_h_fermion */














// kalei ayth
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