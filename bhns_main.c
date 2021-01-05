/*
// Alireza Rashti
// January 2021
*/

#include "bhns_main.h"


/* initial data for BH-NS binary system */
int BH_NS_Binary_Initial_Data(void *vp)
{
  /* if this is a BAM call */
  if (strcmp_i(PgetsEZ(P_"bam_export_id"),"yes"))
    bhns_bam_exporting_initial_data(vp);
  
  /* otherwise construct initial data */
  else
    construct_initial_data(vp);
  
  return EXIT_SUCCESS;
}

/* main algorithm for construction of ID */
static void construct_initial_data(void *vp)
{
  FUNC_TIC
  
  Physics_T *new_phys = 0,*old_phys = 0;
  int Stop = 0;/* if 1 it stops the main loop */
  Uint iter = 0;
  
  /* setting the default parameters */
  set_default_parameters();
  
  /* main iteration loop */
  Stop = update_iteration_params(iter,P_,P_"%s_%ux%ux%u");
  while(!Stop)
  {
     printf("{ Outermost iteration %u ...\n",iter);
     
     new_phys = bhns_initialize_new_physics(old_phys);
     
     bhns_analyze(new_phys,(int)iter);
     
     write_checkpoint(new_phys,Pgets(P_"my_directory"));
     
     bhns_solve_equation(new_phys);
     
     free_physics(old_phys);
     
     old_phys = new_phys;
     
     printf("} Outermost iteration %u ==> Done.\n",iter);
     
     iter++;
     
     Stop = update_iteration_params(iter,P_,P_"%s_%ux%ux%u");
  }
  
  /* save */
  write_checkpoint(new_phys,Pgets(P_"my_directory"));
  
  free_physics(new_phys);
  UNUSED(vp);
  FUNC_TOC
}

/* default parameters for this project */
static void set_default_parameters(void)
{
  /* BH-NS parameters:
  // ================== */
 
  /* how far are BH-NS */
  Pset_default(P_"separation","0.");
  
  /* how fast BH-NS angular velocity */
  Pset_default(P_"angular_velocity","0.");
  
  /* how fast BH-NS infall velocity */
  Pset_default(P_"infall_velocity","0.");
 
  /* how to start off:
  // options:
  // "parameter_file" : it reads parameter file and initialize.
  // "checkpoint_file": it uses checkpoint file and initialize. */
  Pset_default(P_"start_off","parameter_file");
  
  /* total number of iterations that have been executed */
  Pseti(P_"iteration_number",0);
  
  /* system center of mass. */
  Pset_default(P_"x_CM","0."); 
  Pset_default(P_"y_CM","0."); 
  Pset_default(P_"z_CM","0."); 
  
  /* boost velocity for BHNS */
  Pset_default(P_"boost_Vx","0."); 
  Pset_default(P_"boost_Vy","0."); 
  Pset_default(P_"boost_Vz","0."); 
  
  /* what to print for properties of BHNS, add and separate with comma */
  Pset_default(P_"BHNS_properties","x_CM,y_CM,z_CM");
  
  /* NS paramters:
  // ============= */
  
  /* NS baryonic mass */
  Pset_default("NS_baryonic_mass","0.");
  
  /* NS EoS: */
  Pset_default("NS_EoS_description","0");
  /* [polytropic,piecewise_polytropic] */
  Pset_default("NS_EoS_type","0");
  /* unit: [geo] */
  Pset_default("NS_EoS_unit","geo");
  Pset_default("NS_EoS_K","0");
  Pset_default("NS_EoS_Gamma","0");
  
  /* -> central rho0 */
  Pset_default("NS_rho_center","1E-3");
  
  /* geometrical center of NS.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("NS_center_x","0."); 
  Pset_default("NS_center_y","0."); 
  Pset_default("NS_center_z","0."); 
  
  /* box length at the center of NS */
  Pset_default("grid_NS_central_box_length","1.");
  
  /* spin vector to adjust spin for NS */
  Pset_default("NS_Omega_x","0."); 
  Pset_default("NS_Omega_y","0."); 
  Pset_default("NS_Omega_z","0."); 
  
  /* dimensionless spin */
  Pset_default("NS_chi_x","0."); 
  Pset_default("NS_chi_y","0."); 
  Pset_default("NS_chi_z","0."); 
  
  /* ADM momentum */
  Pset_default("NS_Px_ADM","0."); 
  Pset_default("NS_Py_ADM","0."); 
  Pset_default("NS_Pz_ADM","0."); 
  
  /* ADM angular momentum */
  Pset_default("NS_Jx_ADM","0."); 
  Pset_default("NS_Jy_ADM","0."); 
  Pset_default("NS_Jz_ADM","0."); 
  
  /* what to print for properties of NS, add and separate with comma */
  Pset_default(P_"NS_properties",
   "center_x,center_y,center_z,max_radius,min_radius,"
   "ADM_mass,baryonic_mass,Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM");
  
  /* the very first NS approximation
  // options:
  // ========
  // o. TOV (see NS physics) */
  Pset_default("NS_start_off","TOV"); 
  
  /* max l in Ylm expansion */
  Pset_default("NS_surface_Ylm_max_l","1"); 
  
  /* NS surface type */
  Pset_default("NS_surface_type","perfect_s2"); 
  
  
  
  /* BH paramters:
  // ============= */
  
  /* BH irreducible mass */
  Pset_default("BH_irreducible_mass","0.");
  
  /* geometrical center of BH.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("BH_center_x","0."); 
  Pset_default("BH_center_y","0."); 
  Pset_default("BH_center_z","0."); 
  
  /* box length at the center of BH */
  Pset_default("grid_BH_central_box_length","1.");
  
  /* boost velocity for BH */
  Pset_default("BH_boost_Vx","0."); 
  Pset_default("BH_boost_Vy","0."); 
  Pset_default("BH_boost_Vz","0."); 
  
  /* spin vector to adjust spin for BH */
  Pset_default("BH_Omega_x","0."); 
  Pset_default("BH_Omega_y","0."); 
  Pset_default("BH_Omega_z","0."); 
  
  /* dimensionless spin */
  Pset_default("BH_chi_x","0."); 
  Pset_default("BH_chi_y","0."); 
  Pset_default("BH_chi_z","0."); 
  
  /* ADM momentum */
  Pset_default("BH_Px_ADM","0."); 
  Pset_default("BH_Py_ADM","0."); 
  Pset_default("BH_Pz_ADM","0."); 
  
  /* ADM angular momentum */
  Pset_default("BH_Jx_ADM","0."); 
  Pset_default("BH_Jy_ADM","0."); 
  Pset_default("BH_Jz_ADM","0."); 
  
  /* what to print for properties of BH, add and separate with comma */
  Pset_default(P_"BH_properties",
   "center_x,center_y,center_z,max_radius,min_radius,"
   "irreducible_mass,Christodoulou_mass,Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM");
  
  /* the very first BH approximation
  // options:
  // ========
  // o. CloseKerrSchild (see BH physics) */
  Pset_default("BH_start_off","CloseKerrSchild"); 
  
  /* max l in Ylm expansion */
  Pset_default("BH_surface_Ylm_max_l","1"); 
  
  /* BH surface type */
  Pset_default("BH_surface_type","perfect_s2"); 
  
  
  /* BH filler parameters:
  // { */
  /* how to fill: 
  // options:
  // ChebTn_Ylm_perfect_s2:
  //    fill PERFECT S2 surface excised BH with data 
  //    demanding C2 continuity across horizon. extrapolant is:
  //    f(r,th,ph) = C_{ilm}*ChebyshevT(i,r)*Ylm(th,ph). 
  //    thie method is faster than ChebTn_general_s2
  // ChebTn_general_s2: 
  //    fill a GENERAL S2 surface excised BH with data 
  //    demanding C2 continuity across horizon. extrapolant is:
  //    f(r) = C_{i}ChebyshevT(i,r) along radius. 
  // None:
  //    No filling. */
  Pset_default("BH_filler_method","ChebTn_Ylm_perfect_s2");
  /* max l for Y_{lm} expansion in ChebTn_Ylm_perfect_s2 */
  Pset_default("BH_filler_Ylm_expansion_lmax","15");
  /* it won't fill inside the BH if r < Rmin_cutoff.
  // note: this only kicks in for ChebTn_general_s2 method.
  // the ChebTn_Ylm_perfect_s2 is so fast which we don't bother. */
  Pset_default("BH_filler_Rmin_cutoff","-1.");
  /* which fields to be filled */
  Pset_default("BH_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                                  "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                                  "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                                  "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                                  "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  /* value of extrapolated fields at the center of BH r = 0 */
  Pset_default("BH_filler_r0_trK","0.");
  Pset_default("BH_filler_r0_alpha","0.2");
  Pset_default("BH_filler_r0_psi","2.");
  Pset_default("BH_filler_r0_alphaPsi","0.4");
  Pset_default("BH_filler_r0_beta_U0","0.");
  Pset_default("BH_filler_r0_beta_U1","0.");
  Pset_default("BH_filler_r0_beta_U2","0.");
  Pset_default("BH_filler_r0_gConf_D0D0","4.");
  Pset_default("BH_filler_r0_gConf_D0D1","0.");
  Pset_default("BH_filler_r0_gConf_D0D2","0.");
  Pset_default("BH_filler_r0_gConf_D1D1","4.");
  Pset_default("BH_filler_r0_gConf_D1D2","0.");
  Pset_default("BH_filler_r0_gConf_D2D2","4.");
  Pset_default("BH_filler_r0_adm_Kij_D0D0","0.");
  Pset_default("BH_filler_r0_adm_Kij_D0D1","0.");
  Pset_default("BH_filler_r0_adm_Kij_D0D2","0.");
  Pset_default("BH_filler_r0_adm_Kij_D1D1","0.");
  Pset_default("BH_filler_r0_adm_Kij_D1D2","0.");
  Pset_default("BH_filler_r0_adm_Kij_D2D2","0.");
  /* verbos? [1 = yes; 0 = no] */
  Pset_default("BH_filler_verbose","0");
  
  /* test BH filler: interpolate along this line:
  // x = init_x + t*slop_x
  // y = init_y + t*slop_y
  // z = init_z + t*slop_z. */
  Pset_default("BH_filler_test_print_1d","no");
  Pset_default("BH_filler_test_print_1d_slop_x","0.");
  Pset_default("BH_filler_test_print_1d_slop_y","0.");
  Pset_default("BH_filler_test_print_1d_slop_z","0.");
  Pset_default("BH_filler_test_print_1d_init_x","0.");
  Pset_default("BH_filler_test_print_1d_init_y","0.");
  Pset_default("BH_filler_test_print_1d_init_z","0.");
  Pset_default("BH_filler_test_print_1d_length","20.");
  Pset_default("BH_filler_test_print_1d_points","1");
  /* } */
}
