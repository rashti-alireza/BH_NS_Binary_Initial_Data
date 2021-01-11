/*
// Alireza Rashti
// January 2021
*/

/* various functions to start off a new physics for BH-NS */


#include "bhns_initialize.h"

/* decide how to initialize the new physics */
Physics_T *bhns_initialize_new_physics(Physics_T *const old_phys)
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
      new_phys = bhns_read_physics_from_checkpoint();
      
    else 
      new_phys = guess_new_physics();
  }
  else/* use old physics, tune it and make new physics */
  {
    new_phys = infer_new_physics(old_phys);
  }
  
  return new_phys;
}


/* use old physics to infer the new physics */
static Physics_T *infer_new_physics(Physics_T *const old_bhns)
{
  FUNC_TIC
  
  Physics_T *const bhns = init_physics(0,BHNS);/* the whole system */
  Physics_T *const bh   = init_physics(bhns,BH);/* BH part */
  Physics_T *const ns   = init_physics(bhns,NS);/* NS part */
  Physics_T *const old_bh = init_physics(old_bhns,BH);/* BH part */
  Physics_T *const old_ns = init_physics(old_bhns,NS);/* NS part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  old_bh->grid_char = grid_char;
  old_bh->igc       = Ibh;
  old_ns->grid_char = grid_char;
  old_ns->igc       = Ins;
  
  /* adjust and tuning */
  Psets("NS_enthalpy_neat","no");
  physics(old_ns,STRESS_ENERGY_UPDATE);
  physics(old_ns,STAR_TUNE_EULER_CONST);
  physics(old_ns,STRESS_ENERGY_UPDATE);
  physics(old_bh,BH_TUNE_RADIUS);
  physics(old_bhns,SYS_TUNE_P_ADM);
  physics(old_ns,STRESS_ENERGY_UPDATE);
  physics(old_ns,STAR_TUNE_FORCE_BALANCE);
  physics(old_ns,STRESS_ENERGY_UPDATE);
  physics(old_ns,STAR_EXTRAPOLATE_MATTERS);
  physics(old_ns,STAR_TUNE_CENTER);
  physics(old_ns,STAR_FIND_SURFACE);
  
  /* new grid */
  create_new_grid(grid_char,bhns);
  bh->grid = bhns->grid;
  ns->grid = bhns->grid;
  
  /* set paramters */
  physics(bhns,FREE_DATA_SET_PARAMS);
  physics(bhns,ADM_SET_PARAMS);
  physics(bhns,SYS_SET_PARAMS);
  physics(bhns,STRESS_ENERGY_SET_PARAMS);
  physics(bhns,OBSERVE_SET_PARAMS);  
  physics(bhns,BH_SET_PARAMS);
  physics(bhns,STAR_SET_PARAMS);
  
  /* add fields */
  physics(bhns,ADM_ADD_FIELDS);
  physics(bhns,FREE_DATA_ADD_FIELDS);
  physics(bhns,STRESS_ENERGY_ADD_FIELDS);
  physics(bhns,SYS_ADD_FIELDS);
  physics(bhns,OBSERVE_ADD_FIELDS);
  physics(bhns,BH_ADD_FIELDS);
  physics(bhns,STAR_ADD_FIELDS);
  
  /* populate fields */
  physics(bhns,FREE_DATA_POPULATE);
  initialize_fields_using_previous_solve(bhns,old_bhns);
  
  /* beta = B0+B1 */
  physics(bhns,ADM_UPDATE_B1I);
  initial_B0I(bhns,".*");
  update_partial_derivatives(bhns,".*","^dB0_U.+,^ddB0_U.+");
  physics(bhns,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(bhns,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  update_partial_derivatives(ns,"NS","^dphi_D.$,^ddphi_D.D.$");
  
  /* update AConf^{ij} */
  physics(bhns,ADM_UPDATE_AConfIJ);
  
  /* update normal on AH */
  physics(bh,BH_UPDATE_sConf);
  
  /* update matter fields */
  Psets("NS_enthalpy_neat","yes");
  physics(ns,STRESS_ENERGY_UPDATE);
  
  /* free */
  free_physics(bh);
  free_physics(ns);
  free_physics(old_bh);
  free_physics(old_ns);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return bhns;
}  

/* use a known BH and NS solution to initialize the physics */
static Physics_T *guess_new_physics(void)
{
  FUNC_TIC
  
  Physics_T *const bhns = init_physics(0,BHNS);/* the whole system */
  Physics_T *const bh   = init_physics(bhns,BH);/* BH part */
  Physics_T *const ns   = init_physics(bhns,NS);/* NS part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  
  /* set paramters */
  physics(bhns,FREE_DATA_SET_PARAMS);
  physics(bhns,ADM_SET_PARAMS);
  physics(bhns,SYS_SET_PARAMS);
  physics(bhns,STRESS_ENERGY_SET_PARAMS);
  physics(bhns,OBSERVE_SET_PARAMS);  
  physics(bhns,BH_SET_PARAMS);
  physics(bhns,STAR_SET_PARAMS);
  
  /* create grid */
  bh->grid_char = grid_char;
  bh->igc       = Ibh;
  ns->grid_char = grid_char;
  ns->igc       = Ins;
  physics(bh,BH_START);
  physics(bh,BH_FIND_SURFACE);
  physics(ns,STAR_START);
  physics(ns,STAR_FIND_SURFACE);
  create_new_grid(grid_char,bhns);
  bh->grid = bhns->grid;
  ns->grid = bhns->grid;
  
  /* add fields */
  physics(bhns,ADM_ADD_FIELDS);
  physics(bhns,FREE_DATA_ADD_FIELDS);
  physics(bhns,STRESS_ENERGY_ADD_FIELDS);
  physics(bhns,SYS_ADD_FIELDS);
  physics(bhns,OBSERVE_ADD_FIELDS);
  physics(bhns,BH_ADD_FIELDS);
  physics(bhns,STAR_ADD_FIELDS);
  
  /* populate fields */
  physics(bhns,FREE_DATA_POPULATE);
  physics(bhns,SYS_INITIALIZE_FIELDS);
  /* beta = B0+B1 */
  physics(bhns,ADM_UPDATE_B1I);
  initial_B0I(bhns,".*");
  update_partial_derivatives(bhns,".*","^dB0_U.+,^ddB0_U.+");
  physics(bhns,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(bhns,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  update_partial_derivatives(ns,"NS","^dphi_D.$,^ddphi_D.D.$");
  
  /* update AConf^{ij} */
  physics(bhns,ADM_UPDATE_AConfIJ);
  
  /* update normal on AH */
  physics(bh,BH_UPDATE_sConf);
  
  /* update stress energy-tensor */
  Psetd("NS_Euler_equation_constant",
        star_NS_current_Euler_eq_const(ns));
  physics(ns,STRESS_ENERGY_UPDATE);
  
  /* free */
  free_physics(bh);
  free_physics(ns);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return bhns;
}

/* based on grid character, make a new grid for the system. */
static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const bhns)
{
  FUNC_TIC
  
  /* a new grid */
  Grid_T *const grid = alloc_grid();
  
  /* grid for characters */
  grid_char->grid = grid;
  
  if (!Pcmps("grid_kind","SplitCubedSpherical(BH+NS)"))
    Error0(NO_OPTION);
  
  grid_char->S              = Pgetd("BHNS_separation");
  
  grid_char->params[Ibh]->l = Pgetd("grid_BH_central_box_length");
  grid_char->params[Ibh]->w = Pgetd("grid_BH_central_box_length");
  grid_char->params[Ibh]->h = Pgetd("grid_BH_central_box_length");
  
  grid_char->params[Ins]->l = Pgetd("grid_NS_central_box_length");
  grid_char->params[Ins]->w = Pgetd("grid_NS_central_box_length");
  grid_char->params[Ins]->h = Pgetd("grid_NS_central_box_length");
  
  set_params_of_split_cubed_spherical_grid(grid_char);
    
  make_patches(grid);
  realize_interfaces(grid);
  
  bhns->grid = grid;
  
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
    READ_v(B1_U0);
    READ_v(B1_U1);
    READ_v(B1_U2);
    
    REALLOC_v_WRITE_v(B0_U0);
    REALLOC_v_WRITE_v(B0_U1);
    REALLOC_v_WRITE_v(B0_U2);
    
    FOR_ALL_ijk
    {
      B0_U0[ijk] = beta_U0[ijk]-B1_U0[ijk];
      B0_U1[ijk] = beta_U1[ijk]-B1_U1[ijk];
      B0_U2[ijk] = beta_U2[ijk]-B1_U2[ijk];
    }
  }
  
  FUNC_TOC
}

/* loading from checkpoint */
Physics_T *bhns_read_physics_from_checkpoint(void)
{
  FUNC_TIC
  Physics_T *const bhns = init_physics(0,BHNS);
  FILE *file = 0;
  
  /* first load grid and parameters */
  file = open_checkpoint_file_then_read_grid_and_params(bhns);
  
  /* make the patches */
  make_patches(bhns->grid);
  
  /* realizing the geometry */
  realize_interfaces(bhns->grid);
  
  /* set parameters, it's important to add paramters 
  // since these call also reposible to set default functions. */
  physics(bhns,FREE_DATA_SET_PARAMS);
  physics(bhns,ADM_SET_PARAMS);
  physics(bhns,SYS_SET_PARAMS);
  physics(bhns,STRESS_ENERGY_SET_PARAMS);
  physics(bhns,OBSERVE_SET_PARAMS);  
  physics(bhns,BH_SET_PARAMS);
  physics(bhns,STAR_SET_PARAMS);
  
  /* now add fields */
  physics(bhns,ADM_ADD_FIELDS);
  physics(bhns,FREE_DATA_ADD_FIELDS);
  physics(bhns,STRESS_ENERGY_ADD_FIELDS);
  physics(bhns,SYS_ADD_FIELDS);
  physics(bhns,OBSERVE_ADD_FIELDS);
  physics(bhns,BH_ADD_FIELDS);
  physics(bhns,STAR_ADD_FIELDS);
    
  /* then read those saved fields */
  read_fields_from_checkpoint_file(bhns,file);
  
  
  {//temp
  Warning("temp");
  physics(bhns,FREE_DATA_POPULATE);
  physics(bhns,SYS_INITIALIZE_FIELDS);
  initial_B0I(bhns,".*");
  physics(bhns,ADM_UPDATE_B1I);
  update_partial_derivatives(bhns,".*","^dB0_U.+,^ddB0_U.+");
  physics(bhns,ADM_UPDATE_beta);
  update_partial_derivatives(bhns,".*","^dpsi_D.+,^ddpsi_D.+,"
                                      "^dalphaPsi_D.+,^ddalphaPsi_D.+");
  physics(bhns,ADM_UPDATE_AConfIJ);
  }
 
  
  Fclose(file);
  
  FUNC_TOC
  return bhns;
}

/* using copy or interpolation from old physics to 
// initialize fields for new physics */
static void initialize_fields_using_previous_solve
            (Physics_T *const new_phys,Physics_T *const old_phys)
{
  FUNC_TIC
  
  /* matter fields */
  interpolate_fields_from_old_grid_to_new_grid
    (mygrid(old_phys,"NS"),mygrid(new_phys,"NS"),"phi,enthalpy",0);
  
  /* if resolution changed */
  if(Pgeti(P_"did_resolution_change?"))
  {
    interpolate_fields_from_old_grid_to_new_grid
      (old_phys->grid,new_phys->grid,"psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
  }
  else
  {
    const char *region = 0;
    if (new_phys->grid->kind == Grid_SplitCubedSpherical_BHNS)
    {
      /* since filling_box,outermost are fixed, only copy */
      region = "filling_box,outermost";
      interpolate_fields_from_old_grid_to_new_grid
        (mygrid(old_phys,region),mygrid(new_phys,region),
         "psi,alphaPsi,B0_U0,B0_U1,B0_U2",1);
      
      region = "NS,NS_around,BH,BH_around";
      interpolate_fields_from_old_grid_to_new_grid
        (mygrid(old_phys,region),mygrid(new_phys,region),
         "psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
    }
    else
      Error0(NO_OPTION);
  }
  
  FUNC_TOC
}

