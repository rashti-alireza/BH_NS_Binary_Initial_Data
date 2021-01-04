/*
// Alireza Rashti
// December 2020
*/

/* solving physics equations (elliptic solve) */

#include "sbh_solve_eqs.h"

/* setup and issue physics solve */
void sbh_solve_equation(Physics_T *const phys)
{
  FUNC_TIC
  
  /* populate value of inner BC */
  Physics_T *const bh = init_physics(phys,BH);
  physics(bh,BH_UPDATE_INNER_BC);
  free_physics(bh);
  
  /* initialize */
  physics(phys,EQ_SET_PARAMS);
  physics(phys,EQ_ADD_FIELDS);
  
  /* setup external eq functions
  // NOTE: it must be after EQ_SET_PARAMS */
  eq_field_update  = field_update;
  eq_source_update = source_update;
  eq_stop_criteria = stop_criteria;
  eq_analyze_solution = sbh_analyze;

  /* solve */
  physics(phys,EQ_SOLVE);
  
  FUNC_TOC
}


/* stop criteria for solver, namely, if some conditions satisfied, 
// it stops. one can give specific criteria according to the given field's name.
// ->return value: 0 means stop, 1 means continue to solve */
static int stop_criteria(Grid_T *const grid,const char *const name)
{
  int stop = 1;
  int stop_max = 1;
  int stop_res = 0;
  int stop_backtrack = 1;
  int stop_abnormal  = 1;
  const double res_d    = Pgetd("solve_residual");/* desired residual */
  const int max_step    = Pgeti("solve_max_Newton_step");
  const double res_fac  = Pgetd("solve_residual_factor");
  const Uint npatch = grid->np;
  Uint p;
  
  /* if no step should be taken */
  if (max_step  == 0)
  {
    printf("%s equation:\n"
           Pretty0"Newton solver reached maximum step number so existing ...\n",name);
    fflush(stdout);
    return 0;
  }
    
  
  /* NOTE: due to the break command, the order of ifs are important */
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    double res      = patch->solving_man->Frms;/* current residual */
    int solver_step = patch->solving_man->settings->solver_step;/* iteration number */
    
    /* if nan or inf */
    if (!isfinite(res))
    {
      stop_abnormal = 0;
      break;
    }
    
    /* if this is the very first step, don't check the following */
    if (solver_step  == 0)
      continue;
    
    /* note: all patches have same solver_step */
    if (solver_step >= max_step)
    {
      stop_max = 0;
      break;
    }
  }
  
  if (!stop_abnormal)
  {
    printf("%s equation:\n"
           Pretty0"Newton solver got abnormal residual so exit ...\n",name);
    fflush(stdout);
    return stop_abnormal;
  }
  
  if (!stop_backtrack)
  {
    printf("%s equation:\n"
           Pretty0"Newton solver increased the residual so backtrack and exist ...\n",name);
    fflush(stdout);
    backtrack_solutions(grid,name);
    return stop_backtrack;
  }
  
  if (!stop_max)
  {
    printf("%s equation:\n"
           Pretty0"Newton solver reached maximum step number so existing ...\n",name);
    fflush(stdout);
    return stop_max;
  }
  
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    double res      = patch->solving_man->Frms;/* current residual */
    double res_i    = patch->solving_man->settings->Frms_i;/* initial residual */
    
    /* since one of them is enough to continue */
    if (res > res_d && res > res_fac*res_i)
    {
      stop_res = 1;
      break;
    }
    
  }
  
  if (!stop_res)  
  {
    printf("%s equation:\n"
           Pretty0"Newton solver satisfies demanding residual so existing ...\n",name);
    fflush(stdout);
    return stop_res;
  }
  
  return stop;
}


/* how update sourc */
static void source_update(Grid_T *const grid,const char *const name)
{
  Physics_T *const sbh = init_physics(0,SBH);
  sbh->grid            = grid;
  
  physics(sbh,ADM_UPDATE_AConfIJ);
  
  sbh->grid = 0;/* don't free grid */
  free_physics(sbh);
  
  UNUSED(name);
}

/* how update field */
static void field_update(Patch_T *const patch,const char *const name)
{
  if (!strcmp(name,"psi"))
  {
    partial_derivative_with_regex(patch,"^dpsi_D.$,^ddpsi_D.D.$");
  }
  else if (!strcmp(name,"alphaPsi"))
  {
    partial_derivative_with_regex(patch,"^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  }
  else if (!strcmp(name,"B0_U0"))
  {
    partial_derivative_with_regex(patch,"^dB0_U0D.$,^ddB0_U0D.D.$");
    adm_update_beta_U0(patch);
  }
  else if (!strcmp(name,"B0_U1"))
  {
    partial_derivative_with_regex(patch,"^dB0_U1D.$,^ddB0_U1D.D.$");
    adm_update_beta_U1(patch);
  }
  else if (!strcmp(name,"B0_U2"))
  {
    partial_derivative_with_regex(patch,"^dB0_U2D.$,^ddB0_U2D.D.$");
    adm_update_beta_U2(patch);
  }
  else
    Error0(NO_OPTION);

}

/* restore fields to the last solution */
static void backtrack_solutions(Grid_T *const grid,const char *const name)
{
  const Uint npatch = grid->np;
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npatch; ++p)
  {
    Patch_T *patch  = grid->patch[p];
    Field_T *f      = patch->fields[Ind(name)];
    double *v = f->v;
    const double *last_sol = patch->solving_man->settings->last_sol;
    Uint ijk;
    
    free_coeffs(f);
    for(ijk = 0; ijk < patch->nn; ++ijk)
      v[ijk] = last_sol[ijk];
    
    eq_field_update(patch,name);
  }
  
}
