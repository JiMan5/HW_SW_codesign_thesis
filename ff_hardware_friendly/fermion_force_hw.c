#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

void dump_matrix_array(const char *fname, su3_matrix *arr) {
    FILE *f = fopen(fname, "wb");
    if (!f) { perror(fname); return; }
    fwrite(arr, sizeof(su3_matrix), SITES_ON_NODE, f);
    fclose(f);
}

    /*
    char fname[128];
    snprintf(fname, sizeof(fname), "oprod_bins/hw_oprod_path_newwalk_%03d.bin", ipath);
    dump_matrix_array(fname, oprod_along_path[0]);
    printf("Dumped oprod_along_path[0] for path %d\n", ipath);
    */

    //FOR DEBUG
      /*if(ipath==158){
            printf("path ---> (");
            for(int j=0; j<length; j++){
                printf("%d ", this_path->dir[j]);
            }
            printf(")\n");
        int countertemp = 0;
        for (size_t i = 0; i < 10; ++i) {
          countertemp++;
          int test = i * 12345 % SITES_ON_NODE; // pseudo-random
          //int nbrs = walk_dir(test, dir);
          int x,y,z,t;
          coords_from_site_index(test, &x,&y,&z,&t);
          //coords_from_site_index(nbrs,&xn,&yn,&zn,&tn);
          printf("i=%d site=(%d,%d,%d,%d) "
              "\n",
              test, x,y,z,t);

          printf("oprod_along_path[%d][%d] = ", length, test);print_su3(&oprod_along_path[length][test]);
          if(countertemp==10){
              printf("countertemp = %d\n", countertemp); 
          }
        }
        exit(0);
      }*/

void fermion_force_fn_multi_hw_friendly(
    int *netbackdirs_table,
    Real *residues,               //size NTERMS
    su3_vector **multi_x,         //[NTERMS][SITES_ON_NODE]
    Q_path *q_paths_forward,      //[FORW_Q_PATHS] (only forward paths)
    su3_matrix (*links)[4],       //[SITES_ON_NODE][4]
    anti_hermitmat (*mom)[4]      //[SITES_ON_NODE][4] (packed)
)
{
    //static arrays
    static su3_matrix oprod_along_path[MAX_PATH_LENGTH+1][SITES_ON_NODE];
    static su3_matrix mats_along_path[MAX_PATH_LENGTH+1][SITES_ON_NODE];
    static su3_matrix force_accum[4][SITES_ON_NODE];
    static su3_matrix mat_tmp_work[SITES_ON_NODE]; // work buffer used by transports/mults

    Real ferm_epsilon = 2.0f * EPS;
    su3_matrix tmat;
    su3_matrix tmat2;

    //clear junk data loop
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&force_accum[d][i]);
        }
    }

    //big loop over paths

    for (int ipath = 0; ipath < FORW_Q_PATHS; ++ipath) {
        const Q_path *this_path = &q_paths_forward[ipath];
        int length = this_path->length;
        int dir0  = this_path->dir[0]; //for first link of path later
        int odir0 = OPP_DIR(dir0);
        int dirlast = this_path->dir[length-1]; //for last link of path later
        Real coeff = ferm_epsilon * this_path->coeff;

        //clear junk data
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&oprod_along_path[0][i]);
        }
        
        //loop over terms
        for (int term = 0; term < NTERMS; term++) {
            for (size_t i = 0; i < SITES_ON_NODE; i++) {
                int nbr = walk_dir((int)i, netbackdirs_table[ipath]);
                su3_projector(&multi_x[term][i], &multi_x[term][nbr], &tmat);
                scalar_mult_add_su3_matrix(&oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i]);
            }
        }        


        int j = length - 1; 
        int k = GOES_BACKWARDS(dir0) ? 1 : 0;
        for (int ilink = MAX_PATH_LENGTH - 1; ilink >= 0; --ilink) {
            if (ilink > j) continue;
            if (ilink < k) continue;

            int src_layer = length - ilink - 1;
            int dst_layer = length - ilink;
            int dir = this_path->dir[ilink];
            link_transport_connection(oprod_along_path[src_layer], oprod_along_path[dst_layer], mat_tmp_work, dir, links);
        }

    
        //first link of path
        if (GOES_FORWARDS(dir0)){
            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                int nbr = walk_dir((int)i, OPP_DIR(dir0));
                su3_adjoint(&links[nbr][dir0], &mats_along_path[1][i]);
            }
        } 
        else { //backward means take opposite link directly
            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                mats_along_path[1][i] = links[i][odir0];
            }
        }

        //remaining links 
        k = GOES_FORWARDS(dirlast) ? (length - 1) : length;
        for (int ilink = 1; ilink < MAX_PATH_LENGTH; ++ilink) {
            if (ilink >= k || k == 0) continue;
            int dir = OPP_DIR(this_path->dir[ilink]);
            printf("mpika%d???\n", ipath);
            link_transport_connection( mats_along_path[ilink], mats_along_path[ilink + 1], mat_tmp_work, dir, links);
        }

        /*if(ipath == 98){
            int countertemp = 0;
            printf("path ---> (");
            for(int j = 0; j < length; j++){
                printf("%d ", this_path->dir[j]);
            }
            printf(")\n");

            for(size_t i = 0; i < 10; ++i){
                countertemp++;
                int test = i * 12345 % SITES_ON_NODE;
                int x, y, z, t;
                coords_from_site_index(test, &x, &y, &z, &t);

                printf("i=%d site=(%d,%d,%d,%d)\n", test, x, y, z, t);

                printf("mats_along_path[%d][%d] = ", length, test);
                print_su3(&mats_along_path[length][test]);

                if(countertemp == 10){
                    printf("countertemp = %d\n", countertemp);
                }
            }
            exit(0);
        }*/

        
        //edw we walk along each link of the path and we accumulate contributions into force_accum

        int lastdir_local = NODIR; // no previous direction initially
        int lastdir_saved = (length > 0) ? this_path->dir[length - 1] : NODIR;

        for (int ilink = 0; ilink < MAX_PATH_LENGTH; ++ilink) {

            //for active iterations so that the loop bounds are const
            int active = (ilink < length);
            int is_extra = (!active) && (ilink == length) && GOES_BACKWARDS(lastdir_saved);

             // Skip all iterations beyond the path end (past the extra)
            if (!active && !is_extra)
                continue;

            int local_dir = active ? this_path->dir[ilink] : NODIR;
            Real local_coeff = (ilink & 1) ? -coeff : coeff;

            int need_A = active && GOES_FORWARDS(local_dir);
            int need_B = ( (active && ilink > 0 && GOES_BACKWARDS(lastdir_local)) || is_extra );

            int need_mat = need_A || need_B;

            // ---- Compute mat_tmp_work when needed ----
            if (need_mat) {
                if (ilink == 0) {
                    for (size_t i = 0; i < SITES_ON_NODE; ++i)
                        mat_tmp_work[i] = oprod_along_path[length][i];
                } 
                else {
                    int mat_idx = (active) ? ilink : length;
                    int op_idx = (active) ? (length - ilink) : 0;   // for the extra iteration, use 0 (same as MILC)
                    for (size_t i = 0; i < SITES_ON_NODE; ++i)
                        mult_su3_na(&oprod_along_path[op_idx][i], &mats_along_path[mat_idx][i], &mat_tmp_work[i]);
                }
            }

            // ---- Forward-link contribution ----
            if (need_A) {
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int parity = ((i % NX) + ((i / NX) % NY) + ((i / (NX * NY)) % NZ) + (i / (NX * NY * NZ))) & 1;
                    Real signed_coeff = (parity == 0) ? local_coeff : -local_coeff;
                    scalar_mult_add_su3_matrix(&force_accum[local_dir][i], &mat_tmp_work[i], signed_coeff, &force_accum[local_dir][i]);
                }
            }

            // ---- Backward-link / odir contribution ----
            if (need_B) {
                int odir = (is_extra) ? OPP_DIR(lastdir_saved) : OPP_DIR(lastdir_local);
                //when we’re in the extra final iteration, use parity-flipped sign as usual
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int parity = ((i % NX) + ((i / NX) % NY) + ((i / (NX * NY)) % NZ) + (i / (NX * NY * NZ))) & 1;
                    Real signed_coeff = (parity == 0) ? -local_coeff : local_coeff;
                    scalar_mult_add_su3_matrix(&force_accum[odir][i], &mat_tmp_work[i], signed_coeff, &force_accum[odir][i]);
                }
            }

            // ---- Update lastdir_local only while inside the real path ----
            if (active)
                lastdir_local = local_dir;
        }

    } // end ipath loop

    /*
    // 2) Add force_accum into mom (uncompress -> add -> compress)
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            uncompress_anti_hermitian(&mom[i][d], &tmat2);
            add_su3_matrix(&tmat2, &force_accum[d][i], &tmat2);
            make_anti_hermitian(&tmat2, &mom[i][d]);
        }
    }*/
    }
}

