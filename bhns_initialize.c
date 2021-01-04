/*
// Alireza Rashti
// November 2020
*/

/* various functions for start off a new physics for single BH */


#include "sbh_initialize.h"

/* decide how to initialize the new physics */
Physics_T *sbh_initialize_new_physics(Physics_T *const old_phys)
{
  Physics_T *new_phys = 0;
  if (!old_phys)/* if empty, come up with a start off */
  {
    /* if we wanna use checkpoint file */
    if (Pcmps(P_"start_off","checkpoint_file"))
      Error0(NO_OPTION);
      //new_phys = read_physics_from_checkpoint();
    
    /* can we resume from a useful checkpoint file */
    else if (can_we_use_checkpoint(Pgets("top_directory")))
      new_phys = sbh_read_physics_from_checkpoint();
      
    else 
      new_phys = guess_new_physics();
  }
  else/* use old physics, tune it and make new physics */
  {
    Error0(NO_OPTION);
    //new_phys = infer_new_physics(old_phys);
  }
  
  return new_phys;
}


/* use Kerr-Schild solution to initialize the physics */
static Physics_T *guess_new_physics(void)
{
  FUNC_TIC
  
  Physics_T *const sbh = init_physics(0,SBH);/* the whole system */
  Physics_T *const bh  = init_physics(sbh,BH);/* BH part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  
  /* set paramters */
  physics(sbh,FREE_DATA_SET_PARAMS);
  physics(sbh,ADM_SET_PARAMS);
  physics(sbh,BH_SET_PARAMS);
  physics(sbh,SYS_SET_PARAMS);
  physics(sbh,STRESS_ENERGY_SET_PARAMS);
  physics(sbh,OBSERVE_SET_PARAMS);  
  
  /* create grid */
  bh->grid_char = grid_char;
  bh->igc       = Ibh;
  physics(bh,BH_START);
  physics(bh,BH_FIND_SURFACE);
  create_new_grid(grid_char,sbh);
  bh->grid = sbh->grid;
  
  /* add fields */
  physics(sbh,ADM_ADD_FIELDS);
  physics(sbh,BH_ADD_FIELDS);
  physics(sbh,FREE_DATA_ADD_FIELDS);
  physics(sbh,STRESS_ENERGY_ADD_FIELDS);
  physics(sbh,SYS_ADD_FIELDS);
  physics(sbh,OBSERVE_ADD_FIELDS);
  
  /* populate fields */
  physics(sbh,FREE_DATA_POPULATE);
  physics(sbh,SYS_INITIALIZE_FIELDS);
  /* beta = B0+B1 */
  initial_B0I(sbh,".*");
  physics(sbh,ADM_UPDATE_B1I);
  update_partial_derivatives(sbh,".*","^dB0_U.+,^ddB0_U.+");
  physics(sbh,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(sbh,".*","^dpsi_D.+,^ddpsi_D.+,"
                                      "^dalphaPsi_D.+,^ddalphaPsi_D.+");
  
  /* update AConf^{ij} */
  physics(sbh,ADM_UPDATE_AConfIJ);
  
  /* update normal on AH */
  physics(bh,BH_UPDATE_sConf);
  
  /* free */
  free_physics(bh);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return sbh;
}

/* based on grid character, make a new grid for the system. */
static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const sbh)
{
  FUNC_TIC
  
  /* a new grid */
  Grid_T *const grid = alloc_grid();
  
  /* grid for characters */
  grid_char->grid = grid;
  
  if (Pcmps("grid_kind","SplitCubedSpherical(BH)"))
  {
    grid_char->S              = Pgetd("grid_around_box_length");
    grid_char->params[Ibh]->l = Pgetd("grid_central_box_length");
    grid_char->params[Ibh]->w = Pgetd("grid_central_box_length");
    grid_char->params[Ibh]->h = Pgetd("grid_central_box_length");
    
    set_params_of_split_cubed_spherical_grid(grid_char);
  }
  else
    Error0(NO_OPTION);
    
  make_patches(grid);
  realize_interfaces(grid);
  
  sbh->grid = grid;
  
  FUNC_TOC
}

/* update partial derivatives of the given field name regex match 
// at the given region */
static void update_partial_derivatives(Physics_T *const phys,
                                       const char *const region,
                                       const char *const regex)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  printf(Pretty0"%s\n",regex);
  fflush(stdout);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    partial_derivative_with_regex(patch,regex);
  }
  
  FUNC_TOC
}

/* initial B0^i in beta = B0+B1 */
static void initial_B0I(Physics_T *const phys,
                       const char *const region)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    READ_v(beta_U0);
    READ_v(beta_U1);
    READ_v(beta_U2);
    
    REALLOC_v_WRITE_v(B0_U0);
    REALLOC_v_WRITE_v(B0_U1);
    REALLOC_v_WRITE_v(B0_U2);
    
    FOR_ALL_ijk
    {
      B0_U0[ijk] = beta_U0[ijk];
      B0_U1[ijk] = beta_U1[ijk];
      B0_U2[ijk] = beta_U2[ijk];
    }
  }
  
  FUNC_TOC
}

/* loading from checkpoint */
Physics_T *sbh_read_physics_from_checkpoint(void)
{
  FUNC_TIC
  Physics_T *const sbh = init_physics(0,SBH);
  FILE *file = 0;
  
  /* first load grid and parameters */
  file = open_checkpoint_file_then_read_grid_and_params(sbh);
  
  /* make the patches */
  make_patches(sbh->grid);
  
  /* realizing the geometry */
  realize_interfaces(sbh->grid);
  
  /* set parameters, it's important to add paramters 
  // since these call also reposible to set default functions. */
  physics(sbh,FREE_DATA_SET_PARAMS);
  physics(sbh,ADM_SET_PARAMS);
  physics(sbh,BH_SET_PARAMS);
  physics(sbh,SYS_SET_PARAMS);
  physics(sbh,STRESS_ENERGY_SET_PARAMS);
  physics(sbh,OBSERVE_SET_PARAMS);

  /* now add fields */
  physics(sbh,ADM_ADD_FIELDS);
  physics(sbh,BH_ADD_FIELDS);
  physics(sbh,FREE_DATA_ADD_FIELDS);
  physics(sbh,STRESS_ENERGY_ADD_FIELDS);
  physics(sbh,SYS_ADD_FIELDS);
  physics(sbh,OBSERVE_ADD_FIELDS);
  
  /* then read those saved fields */
  read_fields_from_checkpoint_file(sbh,file);
  
  
  {//temp
  Warning("temp");
  physics(sbh,FREE_DATA_POPULATE);
  physics(sbh,SYS_INITIALIZE_FIELDS);
  initial_B0I(sbh,".*");
  physics(sbh,ADM_UPDATE_B1I);
  update_partial_derivatives(sbh,".*","^dB0_U.+,^ddB0_U.+");
  physics(sbh,ADM_UPDATE_beta);
  update_partial_derivatives(sbh,".*","^dpsi_D.+,^ddpsi_D.+,"
                                      "^dalphaPsi_D.+,^ddalphaPsi_D.+");
  physics(sbh,ADM_UPDATE_AConfIJ);
  }
  
  Fclose(file);
  
  FUNC_TOC
  return sbh;
}

