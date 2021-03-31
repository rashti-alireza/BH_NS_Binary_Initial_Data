#include "bhns_header.h"
#include "physics_equation_lib.h"
#include "maths_equation_solvings_lib.h"
#include "physics_adm_lib.h"


/* external functions, defined in "eq_main.h" */
extern fFunc_stop_criteria_T(*eq_stop_criteria);
extern fFunc_source_update_T(*eq_source_update);
extern fFunc_field_update_T(*eq_field_update);
extern fFunc_analyze_solution_T(*eq_analyze_solution);


void bhns_solve_equation(Physics_T *const phys);
static void backtrack_solutions(Grid_T *const grid,const char *const name);
static int stop_criteria (Grid_T *const grid,const char *const name);
static void source_update(Grid_T *const grid,const char *const name);
static void field_update(Patch_T *const patch,const char *const name);
static void prepare_dFdu(Physics_T *const phys);


