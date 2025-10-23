#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

void fermion_force_fn_multi_hw_friendly(
    Real *residues,               //size NTERMS
    su3_vector **multi_x,         //[NTERMS][SITES_ON_NODE]
    Q_path *q_paths_forward,      //[FORW_Q_PATHS] (only forward paths)
    const int path_axis[FORW_Q_PATHS],
    const int path_steps[FORW_Q_PATHS],
    const int path_sign[FORW_Q_PATHS],
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
        int dir0  = this_path->dir[0]; //for first link of path later
        int axis  = path_axis[ipath];
        int steps = path_steps[ipath];
        int sign  = path_sign[ipath];
        int length = this_path->length;
        Real coeff = ferm_epsilon * this_path->coeff;

        
        //clear junk data
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&oprod_along_path[0][i]);
        }
        
        //loop over terms
        for (int term = 0; term < NTERMS; term++) {
            for (size_t i = 0; i < SITES_ON_NODE; i++) {
                int nbr = neighbor_index_axis((int)i, axis, steps, sign);
                su3_projector(&multi_x[term][i], &multi_x[term][nbr], &tmat);
                scalar_mult_add_su3_matrix(&oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i]);
            }
        }
        
        //hw friendly oso ginetai static for loop path. Tha mporoysa isws na kanw presort ta paths me to length toys akrivws kai na treksw 3 diaforetikes loopes. TBD
        for (int ilink = 0; ilink < MAX_PATH_LENGTH; ++ilink) {
            if (ilink < length) {
                int dir = this_path->dir[ilink];
                link_transport_connection(oprod_along_path[ilink], oprod_along_path[ilink+1], mat_tmp_work, dir, links);
            }
        }

        
        //first link of path
        if (GOES_FORWARDS(dir0)){
            //forward means use adjoint of link from neighbor in opposite direction
            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                //step one site backward along this axis
                int nbr = neighbor_index_axis((int)i, dir0, 1, -1);
                su3_adjoint(&links[nbr][dir0], &mats_along_path[1][i]);
            }
        } 
        else { //backward means take opposite link directly
            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                mats_along_path[1][i] = links[i][OPP_DIR(dir0)];
            }
        }

        //remaining links 
        for (int ilink = 1; ilink < MAX_PATH_LENGTH; ++ilink) {
            if (ilink < length) {
                int dir = OPP_DIR(this_path->dir[ilink]);
                link_transport_connection( mats_along_path[ilink], mats_along_path[ilink + 1], mat_tmp_work, dir, links);
            }
        }

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

    // 2) Add force_accum into mom (uncompress -> add -> compress)
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            uncompress_anti_hermitian(&mom[i][d], &tmat2);
            add_su3_matrix(&tmat2, &force_accum[d][i], &tmat2);
            make_anti_hermitian(&tmat2, &mom[i][d]);
        }
    }
}

