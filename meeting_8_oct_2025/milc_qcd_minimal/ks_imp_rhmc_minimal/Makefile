# Compiler and Flags
CC = gcc
CFLAGS = -O3 -g -std=c99 -Wall -Wno-unused-variable -Wno-unused-but-set-variable \
-DVTUNE -g -fopt-info-vec -fopt-info-optimized \
-DGAUSS_SMEAR_KS_TWOLINK -DSINGLE -DNERSC_TIME -DCGTIME -DFFTIME -DFLTIME -DGFTIME \
-DREMAP -DPRTIME -DIOTIME -DWMTIME -DCG_OK -DREMAP_STDIO_APPEND -DDBLSTORE_FN \
-DFEWSUMS -DD_FN_GATHER13 -DKS_MULTICG=HYBRID -DKS_MULTIFF=FNMAT \
-DMILC_CODE_VERSION="\"6be2\"" -DC_GLOBAL_INLINE -DMILC_PRECISION=1 \
-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DMAX_MIXED -DFN -DHMC -DINT_ALG=INT_OMELYAN \
-I. -I../include -I/opt/intel/oneapi/vtune/2025.1/include

CFLAGS += -MD

# List of object files from build log
OBJS = \
  d_action_rhmc.o eo_fermion_force_rhmc.o gauge_info.o setup.o update_h_rhmc.o update_rhmc.o update_u.o control.o \
  gauge_action_imp_ks.o gauge_force_imp_ks.o gauge_measure_ks.o grsource_rhmc.o ks_ratinv.o load_rhmc_params.o reunitarize_ks.o ploop3_ks.o \
  ape_smear.o check_unitarity.o d_plaq4.o gaugefix2.o io_lat4.o momentum_twist.o nersc_cksum.o path_product.o project_su3_hit.o \
  reunitarize2.o show_generic_opts.o show_scidac_opts.o layout_hyper_prime.o blind_data.o field_translation.o field_utilities.o \
  gauge_utilities.o io_detect.o io_helpers.o io_lat_utils.o make_lattice.o mmap_cache.o ranstuff.o remap_stdio_from_args.o io_ansi.o \
  com_vanilla.o ploop3.o show_generic_ks_md_opts.o \
  charge_utilities.o fermion_links_from_site.o f_meas.o gauss_smear_ks.o gauss_smear_ks_cpu.o gauss_smear_ks_QUDA.o grsource_imp.o \
  naik_eps_utilities.o path_transport.o rephase.o show_generic_ks_opts.o show_hisq_links_opts.o \
  fermion_links_milc.o fermion_links.o fermion_links_fn_load_milc.o fermion_links_fn_twist_milc.o general_staple.o fn_links_milc.o ks_action_paths.o \
  d_congrad5_two_src.o d_congrad5_fn_milc.o mat_invert.o ks_invert.o d_congrad5_fn.o d_congrad_opt.o report_invert_status.o dslash_fn_dblstore.o \
  ks_multicg.o ks_multicg_offset.o fermion_force_asqtad.o fermion_force_multi.o fermion_force_fn_multi.o ff_opt.o \
  gauge_force_imp.o gauge_stuff.o gauge_action_imp.o ranmom.o

# Source directories
VPATH = . ../generic ../generic_ks ../generic_pg

# Libraries
LIBS = ../libraries/su3.1.a ../libraries/complex.1.a -L/opt/intel/oneapi/vtune/2025.1/lib64 -lm

# Default target
all: su3_omelyan_rhmc

# Link the final executable
su3_omelyan_rhmc: $(OBJS)
	$(CC) -o $@ -g -littnotify $(OBJS) $(LIBS)

# Compilation rule
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean
clean:
	rm -f *.o *.d su3_omelyan_rhmc

