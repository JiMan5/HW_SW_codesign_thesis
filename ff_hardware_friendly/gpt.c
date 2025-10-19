// --------------------------- hw-friendly fermion force ---------------------------
// Assumes the following compile-time defines exist in ff_headers.h:
//   #define NTERMS_CONST 8
//   #define SITES_ON_NODE_CONST 131072UL
//   #define NUM_Q_PATHS_CONST 688
//   #define NX 16
//   #define NY 16
//   #define NZ 16
//   #define NT 32
//
// Also assumes presence of SU(3) helper prototypes in ff_general_funcs.c:
//   clear_su3mat, su3_projector, scalar_mult_add_su3_matrix,
//   mult_su3_nn, mult_su3_na, mult_su3_an, su3_adjoint,
//   uncompress_anti_hermitian, add_su3_matrix, make_anti_hermitian
//
// And utility: find_backwards_gather(Q_path *path) (or we inline simple equivalent).
// ---------------------------------------------------------------------------------

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

// parity helper: 0 for even, 1 for odd
static inline int site_parity(int i)
{
    int x = i % NX;
    int y = (i / NX) % NY;
    int z = (i / (NX * NY)) % NZ;
    int t = i / (NX * NY * NZ);
    return ( (x + y + z + t) & 1 );
}

/*
 * Hardware-friendly version of fermion_force_fn_multi.
 *
 * Inputs:
 *  - residues: array[ NTERMS_CONST ] of Real
 *  - multi_x:  array[ NTERMS_CONST ][ SITES_ON_NODE_CONST ] (provided as su3_vector **)
 *  - q_paths:  array[ NUM_Q_PATHS_CONST ] of Q_path (path definitions)
 *  - links:    array[ SITES_ON_NODE_CONST ][4] of su3_matrix (links[i][dir])
 *  - mom:      array[ SITES_ON_NODE_CONST ][4] of anti_hermitmat (momenta in anti-hermitian packed form)
 *
 * Notes:
 *  - This function performs the same arithmetic as the original MILC first-call,
 *    but with deterministic loops and static storage suitable for HLS/prototyping.
 *  - It expects the multi_x and links arrays to be indexed by the same site ordering
 *    (x-fastest), i.e. index = x + NX*(y + NY*(z + NZ*t)).
 */
void fermion_force_fn_multi_hw_friendly(
    Real *residues,
    su3_vector **multi_x,
    Q_path *q_paths,
    su3_matrix (*links)[4],
    anti_hermitmat (*mom)[4]
){
    // static scratch buffers (big — placed in static/global memory)
    static su3_matrix oprod_along_path[MAX_PATH_LENGTH+1][SITES_ON_NODE_CONST];
    static su3_matrix mats_along_path[MAX_PATH_LENGTH+1][SITES_ON_NODE_CONST];
    static su3_matrix force_accum[4][SITES_ON_NODE_CONST];
    static su3_matrix mat_tmp0[SITES_ON_NODE_CONST];

    // locals
    Real ferm_epsilon = 2.0f * EPS_CONST; // keep original scaling
    int ipath, term, ilink;
    int length, dir, odir;
    Q_path *this_path;
    int lastdir = -99;
    int last_netbackdir = NODIR;

    // --- initialize force_accum to zero ---
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
            clear_su3mat(&force_accum[d][i]);
        }
    }

    // Main loop over paths (no sorting, deterministic order contained in q_paths)
    for (ipath = 0; ipath < NUM_Q_PATHS_CONST; ++ipath) {
        this_path = &q_paths[ipath];
        if (this_path->forwback == -1) continue; // skip backwards Dslash entries

        length = this_path->length;

        // compute netbackdir (use existing function or inline logic)
        int netbackdir = find_backwards_gather(this_path);

        // Build outer-product transported to end-of-path for EVERY site if net displacement changed
        if (netbackdir != last_netbackdir) {
            // zero oprod_along_path[0]
            for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                clear_su3mat(&oprod_along_path[0][i]);
            }

            // for each term, accumulate residues[term] * projector( multi_x[term][i], multi_x[term][i+netback] )
            for (term = 0; term < NTERMS_CONST; ++term) {
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int nbr = neighbor_index((int)i, netbackdir);
                    su3_matrix tmat;
                    // su3_projector(a, b, &tmat) computes outer product-like projector between vectors
                    su3_projector(&multi_x[term][i], &multi_x[term][nbr], &tmat);
                    // oprod_along_path[0][i] += residues[term] * tmat
                    scalar_mult_add_su3_matrix(&oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i]);
                }
            }
            last_netbackdir = netbackdir;
        }

        // Recompute transported outer products along the path (backwards transport)
        // Original logic reused parts when consecutive paths shared prefix; we recompute fully for clarity.
        // oprod_along_path[length] is the outer product at the path start end (indexing reversed in original)
        // We iterate ilink from length-1 down to 0 and compute transport steps.
        // We'll compute: for m = 0..length : oprod_along_path[m] corresponds to transported quantity
        // Start with oprod_along_path[0] already computed (at end); now compute oprod_along_path[1..length]
        for (int m = 1; m <= length; ++m) {
            // compute transport for link index = (length - m)  (this corresponds to original ordering)
            int step_dir = this_path->dir[length - m];
            // transport: if GOES_FORWARDS(dir) dest[i] = links[i][dir] * src[neighbor_index(i,dir)]
            // else dest[i] = work gathered from adjoint( links[i][OPP_DIR(dir)] ) * src[i] then shifted by neighbor
            if (GOES_FORWARDS(step_dir)) {
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int nbr = neighbor_index((int)i, step_dir);
                    // dest = links[i][step_dir] * src[nbr]
                    mult_su3_nn(&links[i][step_dir], &oprod_along_path[m-1][nbr], &oprod_along_path[m][i]);
                }
            } else {
                int opp = OPP_DIR(step_dir);
                // compute work[i] = adjoint( links[i][opp] ) * src[i]
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    mult_su3_an(&links[i][opp], &oprod_along_path[m-1][i], &mat_tmp0[i]);
                }
                // now dest[i] = work[ neighbor_index(i, step_dir) ]
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int nbr = neighbor_index((int)i, step_dir);
                    oprod_along_path[m][i] = mat_tmp0[nbr];
                }
            }
        } // end building oprod_along_path up to [length]

        // Build mats_along_path: product of links along path up to each ilink (we will use OPP_DIR in transport)
        // mats_along_path[1] is special case (first link)
        // We follow MILC's semantics:
        // If first link goes forwards, then mats_along_path[1][i] = adjoint( link at neighbor in OPP_DIR(dir) )[i]
        // else mats_along_path[1][i] = links[i][OPP_DIR(dir)]
        {
            int dir0 = this_path->dir[0];
            if (GOES_FORWARDS(dir0)) {
                int opp = OPP_DIR(dir0);
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int nbr = neighbor_index((int)i, opp);
                    // adjoint of links[nbr][opp]
                    su3_adjoint(&links[nbr][opp], &mats_along_path[1][i]);
                }
            } else {
                int opp = OPP_DIR(dir0);
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    mats_along_path[1][i] = links[i][opp]; // copy
                }
            }
        }

        // For ilink > 1, transport mats_along_path forward along remaining path steps.
        for (ilink = 1; ilink < length; ++ilink) {
            int step_dir = OPP_DIR(this_path->dir[ilink]); // as in original
            if (GOES_FORWARDS(step_dir)) {
                // mats_along_path[ilink+1][i] = links[i][step_dir] * mats_along_path[ilink][ neighbor ]
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int nbr = neighbor_index((int)i, step_dir);
                    mult_su3_nn(&links[i][step_dir], &mats_along_path[ilink][nbr], &mats_along_path[ilink+1][i]);
                }
            } else {
                int opp = OPP_DIR(step_dir);
                // work[i] = adjoint( links[i][opp] ) * mats_along_path[ilink][i]
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    mult_su3_an(&links[i][opp], &mats_along_path[ilink][i], &mat_tmp0[i]);
                }
                // mats_along_path[ilink+1][i] = work[ neighbor_index(i, step_dir) ]
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int nbr = neighbor_index((int)i, step_dir);
                    mats_along_path[ilink+1][i] = mat_tmp0[nbr];
                }
            }
        }

        // Now loop over links (and points along path) to compute mat_tmp0 and accumulate forces
        int k_for_lastdir = (GOES_FORWARDS(this_path->dir[0]) ? 0 : 1); // original code's k initial guess (not needed now)
        int lastdir_local = -99; // track lastdir inside this path
        for (ilink = 0; ilink <= length; ++ilink) {
            if (ilink < length) dir = this_path->dir[ilink];
            else dir = NODIR;

            Real coeff = ferm_epsilon * this_path->coeff;
            if ((ilink & 1) == 1) coeff = -coeff;

            // compute mat_tmp0 per site
            if (ilink == 0 && dir != NODIR && GOES_FORWARDS(dir)) {
                // mat_tmp0[i] = oprod_along_path[length][i]
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    mat_tmp0[i] = oprod_along_path[length][i];
                }
            } else if (ilink > 0) {
                // mat_tmp0[i] = oprod_along_path[length - ilink][i] * mats_along_path[ilink][i]^dag ?
                // original used mult_su3_na( &(oprod_along_path[length-ilink][i]),  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    mult_su3_na(&oprod_along_path[length - ilink][i], &mats_along_path[ilink][i], &mat_tmp0[i]);
                }
            }

            // add contribution to the force accumulators
            if (ilink < length && dir != NODIR && GOES_FORWARDS(dir)) {
                // even sites add +coeff, odd sites add -coeff
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int p = site_parity((int)i);
                    Real sign = (p == 0) ? 1.0f : -1.0f;
                    scalar_mult_add_su3_matrix(&force_accum[dir][i], &mat_tmp0[i], coeff * sign, &force_accum[dir][i]);
                }
            }

            if (ilink > 0 && GOES_BACKWARDS(lastdir_local)) {
                odir = OPP_DIR(lastdir_local);
                // even: -coeff ; odd: +coeff
                for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
                    int p = site_parity((int)i);
                    Real sign = (p == 0) ? -1.0f : 1.0f;
                    scalar_mult_add_su3_matrix(&force_accum[odir][i], &mat_tmp0[i], coeff * sign, &force_accum[odir][i]);
                }
            }

            lastdir_local = dir;
        } // end ilink loop

    } // end ipath loop

    // Final: add force accumulators into momenta
    for (int d = XUP; d <= TUP; ++d) {
        for (size_t i = 0; i < SITES_ON_NODE_CONST; ++i) {
            su3_matrix tmat2;
            uncompress_anti_hermitian(&mom[i][d], &tmat2);
            add_su3_matrix(&tmat2, &force_accum[d][i], &tmat2);
            make_anti_hermitian(&tmat2, &mom[i][d]);
        }
    }

    // No mallocs to free (static arrays)
} // fermion_force_fn_multi_hw_friendly
