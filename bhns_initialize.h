#include "sbh_header.h"

/* index of bh */
#define Ibh (0)

Physics_T *sbh_initialize_new_physics(Physics_T *const old_phys);
static Physics_T *guess_new_physics(void);

static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const sbh);


static void update_partial_derivatives(Physics_T *const phys,
                                       const char *const region,
                                       const char *const regex);

static void initial_B0I(Physics_T *const phys,
                       const char *const region);


Physics_T *sbh_read_physics_from_checkpoint(void);

