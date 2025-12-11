#include <stdlib.h>
#include <stdint.h>

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
#define NT_SITES 32

//directions
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

//helper directions for netbackdir
#define X3UP 8
#define Y3UP 9
#define Z3UP 10
#define T3UP 11
#define T3DOWN 12
#define Z3DOWN 13
#define Y3DOWN 14
#define X3DOWN 15

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
