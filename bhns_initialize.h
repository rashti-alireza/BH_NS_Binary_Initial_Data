#include "bhns_header.h"
#include "maths_equation_solvings_lib.h"

/* index of bh */
#define Ibh (0)
/* index of ns */
#define Ins (1)

/* encapsulate use last NS surface for multiple uses */
#define USE_LAST_NS_SURFACE() \
  printf(Pretty0"Using the last NS surface.\n");\
  lmax = (Uint)Pgeti("NS_surface_R|lmax");\
  n    = Ncoeffs_Ylm(lmax);\
  double *realClm = alloc_ClmYlm(lmax);/* freed in free_grid_char */\
  double *imagClm = alloc_ClmYlm(lmax);/* freed in free_grid_char */\
  double *coeffs  = 0;\
  coeffs = Pgetdd("NS_surface_R|realClm");\
  for (Uint ij = 0; ij < n; ++ij)\
    realClm[ij] = coeffs[ij];\
  coeffs = Pgetdd("NS_surface_R|imagClm");\
  for (Uint ij = 0; ij < n; ++ij)\
    imagClm[ij] = coeffs[ij];\
  /* might already have values so free them. */\
  Free(grid_char->params[Ins]->relClm);\
  Free(grid_char->params[Ins]->imgClm);\
  grid_char->params[Ins]->relClm = realClm;\
  grid_char->params[Ins]->imgClm = imagClm;\
  grid_char->params[Ins]->lmax   = lmax;\
  grid_char->params[Ins]->r_min  = Pgetd("NS_min_radius");\
  grid_char->params[Ins]->r_max  = Pgetd("NS_max_radius");\
  Pseti("NS_did_NS_surface_change?",0);


Physics_T *bhns_initialize_new_physics(Physics_T *const old_phys);
static Physics_T *guess_new_physics(void);

static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const bhns);


static void update_partial_derivatives(Physics_T *const phys,
                                       const char *const region,
                                       const char *const regex);

static void initial_B0I(Physics_T *const phys,
                       const char *const region);


Physics_T *bhns_read_physics_from_checkpoint(void);
static Physics_T *infer_new_physics(Physics_T *const old_phys);
static void initialize_fields_using_previous_solve
            (Physics_T *const new_phys,Physics_T *const old_phys);
static void update_params(Physics_T *const phys);
static void move_jacobian
            (Physics_T *const new_phys,Physics_T *const old_phys);
