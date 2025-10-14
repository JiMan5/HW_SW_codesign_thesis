#include <stdlib.h>
#include <stdint.h>


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
#define NDIRS 8				/* number of directions */

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
  Real one_link;
  Real naik;
  Real three_staple;
  Real five_staple;
  Real seven_staple;
  Real lepage;
} asqtad_coeffs_t;

typedef struct {
  asqtad_coeffs_t act_path_coeff;    /* For optimized Asqtad action */
  int num_q_paths;                   /* For all actions */
  Q_path *q_paths;                   /* For all actions */
} ks_component_paths;


typedef struct {
  ks_component_paths p;
  int constructed;         /* Boolean */
} ks_action_paths;

typedef struct {
  int twist_in;       /* Whether boundary-twist phases are in */
  Real bdry_phase[4]; /* The boundary-twist phases in current use */
  int r0[4];          /* The effective lattice origin */
} link_phase_info_t;

typedef struct {
  int preserve;
  link_phase_info_t *phase;
  su3_matrix *fat;
  su3_matrix *lng;
  su3_matrix *fatback;  // NULL if unused
  su3_matrix *lngback;  // NULL if unused
  double eps_naik;
  int notify_quda_new_links;
} fn_links_t;

#define imp_ferm_links_t fn_links_t

typedef struct {
  ks_action_paths *ap;     // Paths defining the action
  imp_ferm_links_t *fm;    // The usual links
} fm_ap_links_t;


typedef struct {
  fm_ap_links_t *fm_ap;
  fm_ap_links_t *fm_ap_du0; /* Derivative of links wrto u0 */
} milc_fm_links_t;


typedef milc_fm_links_t ferm_links_generic_t;
typedef struct {
  int want_du0;          /* Do we need to calculate the derivative wrto u0? */
  int want_deps;         /* Do we need to calculate the derivative wrto eps? */
  int want_aux;          /* Do we need to keep the auxiliary HISQ
			    links?  (They are needed for the HISQ
			    fermion force, but not if we just
			    computing propagators.) */
  int want_back;         /* Do we want backward links as well? */
} ferm_links_options_t;

typedef struct {
  ferm_links_options_t options;
  ferm_links_generic_t *flg;
} fermion_links_t;

/*static ks_action_paths *get_fm_ap_links_ap(fm_ap_links_t *al){
  return al->ap;
}

static ks_action_paths *get_milc_fm_ap_links_ap(milc_fm_links_t *al){
  return get_fm_ap_links_ap(al->fm_ap);
}

ks_action_paths *get_action_paths(fermion_links_t *fl){
  return get_milc_fm_ap_links_ap(fl->flg);
}*/


/* msg_tag structure used for gathers */
/* actual structure defined in individual com_*.c files */
typedef struct msg_tag msg_tag;

/* Structure to keep track of outstanding sends and receives */
struct msg_tag {
 int dummy;
};  /* don't need anything */

//io_bin.c
Real read_eps(const char *fname);
int read_nterms(const char *fname);
Real *read_residues(const char *fname, int nterms);
su3_vector **read_multi_x(const char *fname, int nterms, size_t sites_on_node);
Q_path *read_qpaths(const char *fname, int *num_q_paths);
su3_matrix **read_links(const char *fname, size_t *sites_on_node);
anti_hermitmat **read_mom(const char *fname, size_t *sites_on_node);
//void *read_lattice(const char *fname, size_t *sites_on_node, size_t site_size);

//fermion_force_hw.c
void fermion_force_fn_multi_hw(Real eps, Real *residues, su3_vector **multi_x, int nterms, int prec, fermion_links_t *fl, su3_matrix *links, anti_hermitmat *mom, size_t sites_on_node);
