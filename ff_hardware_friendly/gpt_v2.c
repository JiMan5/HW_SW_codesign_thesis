// ---------- Assumes these macros/constants are defined in ff_headers.h ----------
/*
#define EPS          0.010000f
#define NTERMS       8
#define SITES_ON_NODE 131072UL
#define NUM_Q_PATHS  688
#define FORW_Q_PATHS 344
#define NX 16
#define NY 16
#define NZ 16
#define NT 32

#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

#define NODIR -1
#define OPP_DIR(dir) (7 - (dir))
#define GOES_FORWARDS(dir) ((dir) <= TUP)
#define GOES_BACKWARDS(dir) ((dir) > TUP)
*/

#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

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
        default: break;
    }
    return x + NX * (y + NY * (z + NZ * t));
}

/* compute net displacement for path: returns axis (0=X,1=Y,2=Z,3=T), steps>0 and sign (+1 forward, -1 backward).
   If path net displacement is not purely along a single axis, returns axis=-1 and steps=0.
   This mirrors the original find_backwards_gather logic but returns (axis,steps,sign) so the kernel can compute neighbor.
*/
static void compute_net_disp(const Q_path *path, int *axis_out, int *steps_out, int *sign_out) {
    int disp[4] = {0,0,0,0};
    for (int i = 0; i < path->length; ++i) {
        int d = path->dir[i];
        if (GOES_FORWARDS(d)) {
            disp[d] += 1;
        } else {
            // map backward codes 4..7 to axes 0..3 via OPP_DIR
            int opp = OPP_DIR(d);
            disp[opp] -= 1;
        }
    }
    // find the non-zero axis (assume only one axis non-zero in MILC's FN paths)
    int axis = -1;
    int steps = 0;
    int sign = 0;
    for (int a = 0; a < 4; ++a) {
        if (disp[a] != 0) {
            if (axis != -1) { axis = -1; steps = 0; sign = 0; break; } // more than one axis non-zero
            axis = a;
            steps = abs(disp[a]);
            sign = (disp[a] > 0) ? +1 : -1;
        }
    }
    if (axis == -1) {
        *axis_out = -1; *steps_out = 0; *sign_out = 0;
    } else {
        *axis_out = axis; *steps_out = steps; *sign_out = sign;
    }
}

// External helper; DO NOT IMPLEMENT HERE: provide single-node version separately.
// void link_transport_connection(su3_matrix *src, su3_matrix *dest, su3_matrix *work, int dir);

void fermion_force_fn_multi_hw_friendly(
    Real *residues,               // size NTERMS
    su3_vector **multi_x,         // [NTERMS][SITES_ON_NODE]
    Q_path *q_paths_forward,      // [FORW_Q_PATHS] (only forward paths passed)
    su3_matrix (*links)[4],       // [SITES_ON_NODE][4]
    anti_hermitmat (*mom)[4]      // [SITES_ON_NODE][4] (packed)
)
{
    // static scratch arrays to avoid malloc/free
    static su3_matrix oprod_along_path[MAX_PATH_LENGTH+1][SITES_ON_NODE];
    static su3_matrix mats_along_path[MAX_PATH_LENGTH+1][SITES_ON_NODE];
    static su3_matrix force_accum[4][SITES_ON_NODE];
    static su3_matrix mat_tmp_work[SITES_ON_NODE]; // work buffer used by transports/mults
    // note: using SITES_ON_NODE (macro) to make arrays static sized

    Real ferm_epsilon = 2.0f * EPS;
    su3_matrix tmat; (void)tmat; // temp when needed locally
    su3_matrix tmat2;

    // 0) Initialize force_accum to zero matrices
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&force_accum[d][i]);
        }
    }

    // 1) Loop over forward q_paths (assumed q_paths_forward length = FORW_Q_PATHS)
    for (int ipath = 0; ipath < FORW_Q_PATHS; ++ipath) {
        const Q_path *this_path = &q_paths_forward[ipath];
        int length = this_path->length;

        // 1.a) compute net displacement (axis, steps, sign)
        int axis, steps, sign;
        compute_net_disp(this_path, &axis, &steps, &sign);
        // If axis == -1 then path net displacement is not a single-axis simple displacement.
        // In MILC find_backwards_gather would have exited; here we require axis != -1.
        if (axis == -1) {
            // fallback: skip (should not happen for FN paths produced by MILC)
            continue;
        }

        // Determine direction code to step *from end to start*:
        // If net disp sign>0 (path moves forward), the gather should follow the opposite direction (down)
        int axis_up = axis; // 0->XUP,1->YUP,2->ZUP,3->TUP
        int gather_dir;
        if (sign > 0) gather_dir = OPP_DIR(axis_up); // e.g., XUP(+1) -> XDOWN gather
        else gather_dir = axis_up; // net negative -> gather the UP direction

        // 1.b) Build outer-product at the "end" of path for every site:
        // oprod_along_path[0][i] = sum_term residues[term] * projector(multi_x[term][i], multi_x[term][startSite])
        // where startSite is the site at distance 'steps' along gather_dir from i
        // compute startSite by stepping `steps` times in gather_dir from i
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&oprod_along_path[0][i]);
        }
        for (int term = 0; term < NTERMS; ++term) {
            for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                int nbr = (int)i;
                for (int s = 0; s < steps; ++s) {
                    nbr = neighbor_index(nbr, gather_dir);
                }
                // projector: produces an su3_matrix from two vectors
                su3_projector(&multi_x[term][i], &multi_x[term][nbr], &tmat);
                scalar_mult_add_su3_matrix(&oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i]);
            }
        }

        // 1.c) Compute transported outer-products along path:
        // We'll compute oprod_along_path[m] for m=1..length (transported backwards along path)
        // We follow the simple "do transport per link" approach using link_transport_connection() external helper.
        for (int m = 1; m <= length; ++m) {
            // src = oprod_along_path[m-1], dest = oprod_along_path[m]
            // link dir for this step is this_path->dir[length - m] (same indexing as original)
            int step_dir = this_path->dir[length - m];
            // Call external single-node link transport helper:
            // link_transport_connection(src, dest, work, dir)
            link_transport_connection(oprod_along_path[m-1], oprod_along_path[m], mat_tmp_work, step_dir);
        }

        // 1.d) Build mats_along_path: sequence of link products to be used by mult_su3_na
        // mats_along_path[1] (first link) special-case: if first step is forwards, we need adjoint of neighbor link
        {
            int dir0 = this_path->dir[0];
            if (GOES_FORWARDS(dir0)) {
                int opp = OPP_DIR(dir0);
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int nbr = neighbor_index((int)i, opp);
                    su3_adjoint(&links[nbr][opp], &mats_along_path[1][i]);
                }
            } else {
                int opp = OPP_DIR(dir0);
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    mats_along_path[1][i] = links[i][opp];
                }
            }
        }
        // For ilink = 1 .. length-1 compute mats_along_path[ilink+1] via link_transport_connection-like steps
        for (int ilink = 1; ilink < length; ++ilink) {
            int dir_step = OPP_DIR(this_path->dir[ilink]); // same pattern as original logic
            link_transport_connection(mats_along_path[ilink], mats_along_path[ilink+1], mat_tmp_work, dir_step);
        }

        // 1.e) Walk along points on path and accumulate contributions into force_accum
        int lastdir_local = -99;
        for (int ilink = 0; ilink <= length; ++ilink) {
            int local_dir = (ilink < length) ? this_path->dir[ilink] : NODIR;
            Real coeff = ferm_epsilon * this_path->coeff;
            if ((ilink & 1) == 1) coeff = -coeff;

            // compute mat_tmp_work[i] for all sites:
            if (ilink == 0 && local_dir != NODIR && GOES_FORWARDS(local_dir)) {
                // mat_tmp = oprod_along_path[length]
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    mat_tmp_work[i] = oprod_along_path[length][i];
                }
            } else if (ilink > 0) {
                // mat_tmp = oprod_along_path[length - ilink] * mats_along_path[ilink]
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    mult_su3_na(&oprod_along_path[length - ilink][i], &mats_along_path[ilink][i], &mat_tmp_work[i]);
                }
            }

            // add contribution to force accumulators
            if (ilink < length && local_dir != NODIR && GOES_FORWARDS(local_dir)) {
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int p = ((i % NX) + ((i / NX) % NY) + ((i / (NX*NY)) % NZ) + (i / (NX*NY*NZ))) & 1; // parity
                    Real sign = (p == 0) ? 1.0f : -1.0f;
                    scalar_mult_add_su3_matrix(&force_accum[local_dir][i], &mat_tmp_work[i], coeff * sign, &force_accum[local_dir][i]);
                }
            }

            // update opposite-direction contribution if lastdir_local was backwards
            if (ilink > 0 && GOES_BACKWARDS(lastdir_local)) {
                int odir = OPP_DIR(lastdir_local);
                for (size_t i = 0; i < SITES_ON_NODE; ++i) {
                    int p = ((i % NX) + ((i / NX) % NY) + ((i / (NX*NY)) % NZ) + (i / (NX*NY*NZ))) & 1;
                    Real sign = (p == 0) ? -1.0f : 1.0f;
                    scalar_mult_add_su3_matrix(&force_accum[odir][i], &mat_tmp_work[i], coeff * sign, &force_accum[odir][i]);
                }
            }

            lastdir_local = local_dir;
        } // end ilink loop

    } // end ipath loop

    // 2) Add force_accum into mom (uncompress -> add -> compress)
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE; ++i) {
            uncompress_anti_hermitian(&mom[i][d], &tmat2);
            add_su3_matrix(&tmat2, &force_accum[d][i], &tmat2);
            make_anti_hermitian(&tmat2, &mom[i][d]);
        }
    }

    // done (no frees, static arrays)
}
