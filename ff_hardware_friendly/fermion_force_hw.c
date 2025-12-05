#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

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
        int dir_last_in_path = this_path->dir[length-1]; //for last link of path later
        int lastdir = NODIR; //for logic for force accum (milc has it too, can't skip it, den einai software trick, einai physics necessary)
        Real base_coeff = ferm_epsilon * this_path->coeff;

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
        k = GOES_FORWARDS(dir_last_in_path) ? (length - 1) : length;
        for (int ilink = 1; ilink < MAX_PATH_LENGTH; ++ilink) {
            if (ilink >= k || k == 0) continue;
            int dir = OPP_DIR(this_path->dir[ilink]);
            link_transport_connection( mats_along_path[ilink], mats_along_path[ilink + 1], mat_tmp_work, dir, links);
        }


        
        //gia to prwto
        if (GOES_FORWARDS(dir0)) {
            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                mat_tmp_work[i] = oprod_along_path[length][i];
            }

            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                int x,y,z,t;
                coords_from_site_index((int)i,&x,&y,&z,&t);
                int parity = (x+y+z+t)&1;

                Real sign = (parity==0 ? base_coeff : -base_coeff);
                scalar_mult_add_su3_matrix(&force_accum[dir0][i], &mat_tmp_work[i], sign, &force_accum[dir0][i]);
            }
        }

        lastdir = dir0;

        //gia ta alla
        k = GOES_FORWARDS(dir_last_in_path) ? (length - 1) : length;
        for (int ilink = 1; ilink <= MAX_PATH_LENGTH; ++ilink) {
            if (ilink > k) continue;
            
            int dir_val = 0;
            int dir = NODIR;
            if(ilink < length){
                dir = this_path->dir[ilink];
                dir_val = 1;
            } 

            Real coeff = (ilink & 1) ? -base_coeff : base_coeff;

            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                mult_su3_na(&oprod_along_path[length - ilink][i], &mats_along_path[ilink][i], &mat_tmp_work[i]);
            }

            if (dir_val && GOES_FORWARDS(dir)) {
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int x,y,z,t;
                    coords_from_site_index((int)i,&x,&y,&z,&t);
                    int parity=(x+y+z+t)&1;

                    Real sign = (parity==0 ? coeff : -coeff);
                    scalar_mult_add_su3_matrix(&force_accum[dir][i], &mat_tmp_work[i], sign, &force_accum[dir][i]);
                }
            }

            if (GOES_BACKWARDS(lastdir)) {
                int odir = OPP_DIR(lastdir);

                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int x,y,z,t;
                    coords_from_site_index((int)i,&x,&y,&z,&t);
                    int parity=(x+y+z+t)&1;

                    Real sign = (parity==0 ? -coeff : coeff);
                    scalar_mult_add_su3_matrix(&force_accum[odir][i], &mat_tmp_work[i], sign, &force_accum[odir][i]);
                }
            }

            lastdir = dir;
        }

    } // end ipath loop

    
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            uncompress_anti_hermitian(&mom[i][d], &tmat2);
            add_su3_matrix(&tmat2, &force_accum[d][i], &tmat2);
            make_anti_hermitian(&tmat2, &mom[i][d]);
        }
    }
}

