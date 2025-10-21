#include <stdlib.h>
#include <stdint.h>

//fixed consts for run
#define EPS          0.010000f
#define NTERMS       8
#define SITES_ON_NODE 131072UL
#define NUM_Q_PATHS  688
#define FORW_Q_PATHS 344
#define MAX_PATH_LENGTH 7
#define NX 16
#define NY 16
#define NZ 16
#define NT 32

//directions
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

#define NODIR -1  /* not a direction */

#define OPP_DIR(dir)	(7-(dir))	/* Opposite direction */
#define GOES_FORWARDS(dir) ((dir)<=TUP)
#define GOES_BACKWARDS(dir) ((dir)>TUP)
#define NDIRS 8				/* number of directions */

#define CMUL_J(a,b,c) { (c).real = (a).real*(b).real + (a).imag*(b).imag; \
	  	        (c).imag = (a).imag*(b).real - (a).real*(b).imag; }

#define CMUL(a,b,c) { (c).real = (a).real*(b).real - (a).imag*(b).imag; \
		      (c).imag = (a).real*(b).imag + (a).imag*(b).real; }

#define CADD(a,b,c) { (c).real = (a).real + (b).real;  \
		      (c).imag = (a).imag + (b).imag; }

#define CSUM(a,b) { (a).real += (b).real; (a).imag += (b).imag; }

#define CMULJ_(a,b,c) { (c).real = (a).real*(b).real + (a).imag*(b).imag; \
		        (c).imag = (a).real*(b).imag - (a).imag*(b).real; }

#define CONJG(a,b) { (b).real = (a).real; (b).imag = -(a).imag; }

typedef float Real;

typedef struct {   
  float real;	   
  float imag; 
} fcomplex; 

typedef struct { fcomplex e[3][3]; } fsu3_matrix;
typedef struct { fcomplex c[3]; } fsu3_vector;

typedef struct { 
  fcomplex m01,m02,m12; 
  float m00im,m11im,m22im; 
  float space; } fanti_hermitmat;


#define su3_matrix      fsu3_matrix
#define su3_vector      fsu3_vector
#define anti_hermitmat  fanti_hermitmat


#define MAX_PATH_LENGTH 16
typedef struct {
  int dir[MAX_PATH_LENGTH];	/* directions in path */
  int length;		/* length of path */
  Real coeff;	        /* coefficient, including minus sign if backwards */
  Real forwback;	/* +1 if in forward Dslash, -1 if in backward */
} Q_path;

typedef struct {
    int axis;
    int steps;
    int sign;
} PathDisp;

//read files
Real *read_residues(const char *fname);
su3_vector **read_multi_x(const char *fname);
Q_path *read_qpaths(const char *fname);
su3_matrix **read_links(const char *fname);
anti_hermitmat **read_mom(const char *fname);
//void *read_lattice(const char *fname, size_t *sites_on_node, size_t site_size);

//calculations
void clear_su3mat( su3_matrix *dest );
void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c );
void scalar_mult_add_su3_matrix(su3_matrix *a,su3_matrix *b,Real s, su3_matrix *c);
void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_na(  su3_matrix *a, su3_matrix *b, su3_matrix *c );
void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void su3_adjoint( su3_matrix *a, su3_matrix *b );
void uncompress_anti_hermitian( const anti_hermitmat * const mat_antihermit, su3_matrix *mat_su3 );
void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c );
void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3 );

//helpers for path directions and indexing///////////////////////

void compute_net_disp(const Q_path *path, int *axis_out, int *steps_out, int *sign_out);
int neighbor_index_disp(int i, int axis, int steps, int sign);

//fermion_force_hw.c
fermion_force_fn_multi_hw_friendly(residues, multi_x, qpaths_forward, axis, steps, sign, links, mom_before);

//link_transport_connection
void link_transport_connection(const su3_matrix *src, su3_matrix *dest, su3_matrix *work, int dir, su3_matrix (*links)[4]);