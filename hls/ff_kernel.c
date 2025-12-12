#include <stdio.h>
#include <stdlib.h>
#include "ff_headers.h"

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

//link_transport_connection
void link_transport_connection(su3_matrix *src, su3_matrix *dest, su3_matrix *work, int dir, su3_matrix (*links)[4], lookup_t lookups[SITES_ON_NODE]){

	size_t i;
    if (GOES_FORWARDS(dir)) {
        for (i = 0; i < SITES_ON_NODE; ++i) {
        	int nbr = lookups[i].nbr[dir];
            //int nbr = nbr_table[i][dir];
            mult_su3_nn(&links[i][dir], &src[nbr], &dest[i]);
        }
    }
    else {
        int odir = OPP_DIR(dir);
        for (i = 0; i < SITES_ON_NODE; ++i) {
            mult_su3_an(&links[i][odir], &src[i], &work[i]);
        }

        for (i = 0; i < SITES_ON_NODE; ++i) {
        	int nbr = lookups[i].nbr[dir];
            //int nbr = nbr_table[(int)i][dir];
            dest[i] = work[nbr];
        }
    }
}


void fermion_force_fn_multi_hw_friendly(
    int netbackdirs_table[FORW_Q_PATHS],
    Real residues[NTERMS],
    su3_vector multi_x[NTERMS][SITES_ON_NODE],
    Q_path q_paths_forward[FORW_Q_PATHS],
    su3_matrix links[SITES_ON_NODE][4],
    anti_hermitmat mom[SITES_ON_NODE][4],
    lookup_t lookups[SITES_ON_NODE]
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

    int d, ipath, term, ilink;
    size_t i;

    //clear junk data loop
    for (d = XUP; d <= TUP; ++d) {
        for (i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&force_accum[d][i]);
        }
    }

    //big loop over paths

    for (ipath = 0; ipath < FORW_Q_PATHS; ++ipath) {
    	printf("ipath = %d\n", ipath);
        const Q_path *this_path = &q_paths_forward[ipath];
        int length = this_path->length;
        int dir0  = this_path->dir[0]; //for first link of path later
        int odir0 = OPP_DIR(dir0);
        int dir_last_in_path = this_path->dir[length-1]; //for last link of path later
        int lastdir = NODIR; //for logic for force accum (milc has it too, can't skip it, den einai software trick, einai physics necessary)
        Real base_coeff = ferm_epsilon * this_path->coeff;

        //clear junk data
        for (i = 0; i < SITES_ON_NODE; ++i) {
            clear_su3mat(&oprod_along_path[0][i]);
        }

        //loop over terms
        for (term = 0; term < NTERMS; term++) {
            for (i = 0; i < SITES_ON_NODE; i++) {
            	int nbr = lookups[(int)i].nbr[netbackdirs_table[ipath]];
                //int nbr = nbr_table[(int)i][netbackdirs_table[ipath]];
                su3_projector(&multi_x[term][i], &multi_x[term][nbr], &tmat);
                scalar_mult_add_su3_matrix(&oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i]);
            }
        }


        int j = length - 1;
        int k = GOES_BACKWARDS(dir0) ? 1 : 0;
        for (ilink = MAX_PATH_LENGTH - 1; ilink >= 0; --ilink) {
            if (ilink > j) continue;
            if (ilink < k) continue;

            int src_layer = length - ilink - 1;
            int dst_layer = length - ilink;
            int dir = this_path->dir[ilink];
            link_transport_connection(oprod_along_path[src_layer], oprod_along_path[dst_layer], mat_tmp_work, dir, links, lookups);
        }


        //first link of path
        if (GOES_FORWARDS(dir0)){
            for (i = 0; i < SITES_ON_NODE; ++i) {
            	int nbr = lookups[(int)i].nbr[OPP_DIR(dir0)];
                //int nbr = nbr_table[(int)i][OPP_DIR(dir0)];
                su3_adjoint(&links[nbr][dir0], &mats_along_path[1][i]);
            }
        }
        else { //backward means take opposite link directly
            for (i = 0; i < SITES_ON_NODE; ++i) {
                mats_along_path[1][i] = links[i][odir0];
            }
        }
        //remaining links
        k = GOES_FORWARDS(dir_last_in_path) ? (length - 1) : length;
        for (ilink = 1; ilink < MAX_PATH_LENGTH; ++ilink) {
            if (ilink >= k || k == 0) continue;
            int dir = OPP_DIR(this_path->dir[ilink]);
            link_transport_connection( mats_along_path[ilink], mats_along_path[ilink + 1], mat_tmp_work, dir, links, lookups);
        }



        //gia to prwto
        if (GOES_FORWARDS(dir0)) {
            for (i = 0; i < SITES_ON_NODE; ++i) {
                mat_tmp_work[i] = oprod_along_path[length][i];
            }

            for (i = 0; i < SITES_ON_NODE; ++i) {
            	int parity = lookups[i].parity;
                /*int x,y,z,t;
                coords_from_site_index((int)i,&x,&y,&z,&t);
                int parity = (x+y+z+t)&1;*/

                Real sign = (parity==0 ? base_coeff : -base_coeff);
                scalar_mult_add_su3_matrix(&force_accum[dir0][i], &mat_tmp_work[i], sign, &force_accum[dir0][i]);
            }
        }

        lastdir = dir0;

        //gia ta alla
        k = GOES_FORWARDS(dir_last_in_path) ? (length - 1) : length;
        for (ilink = 1; ilink <= MAX_PATH_LENGTH; ++ilink) {
            if (ilink > k) continue;

            int dir_val = 0;
            int dir = NODIR;
            if(ilink < length){
                dir = this_path->dir[ilink];
                dir_val = 1;
            }

            Real coeff = (ilink & 1) ? -base_coeff : base_coeff;

            for (i = 0; i < SITES_ON_NODE; ++i) {
                mult_su3_na(&oprod_along_path[length - ilink][i], &mats_along_path[ilink][i], &mat_tmp_work[i]);
            }

            if (dir_val && GOES_FORWARDS(dir)) {
                for (i = 0; i < SITES_ON_NODE; ++i) {
                	int parity = lookups[i].parity;
                	/*int x,y,z,t;
                	coords_from_site_index((int)i,&x,&y,&z,&t);
                	int parity = (x+y+z+t)&1;*/

                    Real sign = (parity==0 ? coeff : -coeff);
                    scalar_mult_add_su3_matrix(&force_accum[dir][i], &mat_tmp_work[i], sign, &force_accum[dir][i]);
                }
            }

            if (GOES_BACKWARDS(lastdir)) {
                int odir = OPP_DIR(lastdir);

                for (i = 0; i < SITES_ON_NODE; ++i) {
                	int parity = lookups[i].parity;
                	/*int x,y,z,t;
                	coords_from_site_index((int)i,&x,&y,&z,&t);
                	int parity = (x+y+z+t)&1;*/

                    Real sign = (parity==0 ? -coeff : coeff);
                    scalar_mult_add_su3_matrix(&force_accum[odir][i], &mat_tmp_work[i], sign, &force_accum[odir][i]);
                }
            }

            lastdir = dir;
        }

    } // end ipath loop


    for (d = XUP; d <= TUP; ++d) {
        for (i = 0; i < SITES_ON_NODE; ++i) {
            uncompress_anti_hermitian(&mom[i][d], &tmat2);
            add_su3_matrix(&tmat2, &force_accum[d][i], &tmat2);
            make_anti_hermitian(&tmat2, &mom[i][d]);
        }
    }
}
