#include "bhns_header.h"
#include "maths_equation_solvings_lib.h"

/* index of bh */
#define Ibh (0)
/* index of ns */
#define Ins (1)

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
