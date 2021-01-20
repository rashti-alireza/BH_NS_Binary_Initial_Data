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
     
     write_checkpoint(new_phys,Pgets(P_"my_directory"));
     
     bhns_solve_equation(new_phys);
     
     bhns_analyze(new_phys,(int)iter);
     
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
  
  /* masses */
  Pset_default(P_"ADM_mass","0.");
  Pset_default(P_"Komar_mass","0.");
  
  /* what to print for properties of BHNS, add and separate with comma */
  Pset_default(P_"BHNS_properties",
    "separation,x_CM,y_CM,z_CM,ADM_mass,Komar_mass,mass_ratio,"
    "angular_velocity,infall_velocity,"
    "Px_ADM,Py_ADM,Pz_ADM,"
    "Jx_ADM,Jy_ADM,Jz_ADM,"
    "binding_energy,Virial_error");
    
  /* how to tune P_ADM */
  Pset_default(P_"P_ADM_control_method","adjust(x_CM,y_CM)");
  Pset_default(P_"P_ADM_control_update_weight","0.");
  Pset_default(P_"P_ADM_control_tolerance","1E-5");
  
  /* observer method */
  Pset_default(P_"Observe_ADM_P","S+V,Rashti");
  Pset_default(P_"Observe_ADM_J","S+V,Ossokine");
  Pset_default(P_"Observe_ADM_M","S+V,default");
  Pset_default(P_"Observe_Komar_M","S_inf,default");
  
  /* NS paramters:
  // ============= */
  
  /* NS masses */
  Pset_default("NS_baryonic_mass_current","1.");
  Pset_default("NS_baryonic_mass","1.");
  Pset_default("NS_ADM_mass","1.");
  Pset_default("NS_Komar_mass","1.");
  Pset_default("NS_TOV_ADM_mass","1.");
  Pset_default("NS_TOV_radius","1.");
  Pset_default("NS_TOV_compactness","0.");
  Pset_default("NS_mass_shedding_indicator","1.");
  
  /* NS EoS: */
  Pset_default("NS_EoS_description","NA");
  /* [polytropic,piecewise_polytropic] */
  Pset_default("NS_EoS_type","NA");
  /* unit: [geo] */
  Pset_default("NS_EoS_unit","NA");
  Pset_default("NS_EoS_K","NA");
  Pset_default("NS_EoS_Gamma","NA");
  Pset_default("NS_EoS_rho_th","NA");
  
  /* -> central matters */
  Pset_default("NS_rho_center","1E-3");
  Pset_default("NS_pressure_center","1E-3");
  Pset_default("NS_energy_density_center","1E-3");
  
  /* geometrical center of NS.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("NS_center_x","0."); 
  Pset_default("NS_center_y","0."); 
  Pset_default("NS_center_z","0."); 
  Pset_default("NS_x_CM","0."); 
  Pset_default("NS_y_CM","0."); 
  Pset_default("NS_z_CM","0."); 
  
  /* box length at the center of NS */
  Pset_default("grid_NS_central_box_length","1.");
  
  /* spin vector to adjust spin for NS */
  Pset_default("NS_Omega_x","0."); 
  Pset_default("NS_Omega_y","0."); 
  Pset_default("NS_Omega_z","0."); 
  
  /* spin */
  Pset_default("NS_chi_x","0."); 
  Pset_default("NS_chi_y","0."); 
  Pset_default("NS_chi_z","0."); 
  Pset_default("NS_spin_x","0."); 
  Pset_default("NS_spin_y","0."); 
  Pset_default("NS_spin_z","0."); 
  
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
   "center_x,center_y,center_z,x_CM,y_CM,z_CM,"
   "max_radius,min_radius,TOV_radius,"
   "ADM_mass,TOV_ADM_mass,Komar_mass,baryonic_mass_current,"
   "baryonic_mass,mass_shedding_indicator,TOV_compactness,"
   "rho_center,pressure_center,energy_density_center,"
   "enthalpy_L2_residual,Euler_equation_constant,"
   "Omega_x,Omega_y,Omega_z,chi_x,chi_y,chi_z,"
   "spin_x,spin_y,spin_z,"
   "Px_ADM,Py_ADM,Pz_ADM,"
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
  
  /* observe method pertinet to NS */
  Pset_default("NS_Observe_ADM_M","V_obj,default");
  Pset_default("NS_Observe_Komar_M","S_obj,default");
  Pset_default("NS_Observe_baryonic_M","V_obj,default");
  Pset_default("NS_Observe_ADM_P","S_obj,default");
  Pset_default("NS_Observe_ADM_J","S_obj,default");
  Pset_default("NS_Observe_CM","S_obj,default");
  Pset_default("NS_Observe_spin","S_obj,JRP");
  
  /* smooth and polish phi equation close to the surface */
  Pset_default("NS_Eq_phi_polish","0.1");
  
  /* tune and adjust: */
  Pset_default("NS_force_balance_equation","adjust(d/dy:Omega)");
  Pset_default("NS_force_balance_update_weight","0.");
  Pset_default("NS_adjust_center_method","interpolation");
  Pset_default("NS_enthalpy_allowed_residual","1E-5");
  Pset_default("NS_enthalpy_L2_residual","0.");
  
  /* extrapolation of matter fields outside NS:
  // options = [inverse_r2,exp2,poly2,inverse_r2_expmr,
                inverse_r2_expmAr,enthalpy_expmr_phi_inverse_r2]. */
  Pset_default("NS_extrapolate_matter_fields","inverse_r2_expmAr");
  
  /* Euler eq. constant */
  Pset_default("NS_Euler_equation_constant","0.");
  Pset_default("NS_Euler_const_update_weight","1.");
  
  /* NS enhtalpy update weight */
  Pset_default("NS_enthalpy_update_weight","0.");
  /* set enthalpy != 1 on surface to 1
  // options: [yes/no]. */
  Pset_default("NS_enthalpy_neat","yes");
  /* residual of enthalpy root finder */
  Pset_default("NS_enthalpy_residual","0.");
  
  /* root finder pertinent to NS */
  Pset_default("NS_RootFinder_method","Steepest_Descent");
  Pset_default("NS_RootFinder_Tolerance","1E-9");
  Pset_default("NS_RootFinder_Iteration","1000");
  Pset_default("NS_RootFinder_verbose","no");
  
  /* BH paramters:
  // ============= */
  
  /* BH irreducible mass */
  Pset_default("BH_irreducible_mass","0.");
  
  /* geometrical center of BH.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("BH_center_x","0."); 
  Pset_default("BH_center_y","0."); 
  Pset_default("BH_center_z","0."); 
  
  Pset_default("BH_x_CM","0."); 
  Pset_default("BH_y_CM","0."); 
  Pset_default("BH_z_CM","0."); 
  
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
  
  /* spin */
  Pset_default("BH_chi_x","0."); 
  Pset_default("BH_chi_y","0."); 
  Pset_default("BH_chi_z","0."); 
  Pset_default("BH_spin_x","0."); 
  Pset_default("BH_spin_y","0."); 
  Pset_default("BH_spin_z","0."); 
  
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
   "center_x,center_y,center_z,x_CM,y_CM,z_CM,max_radius,min_radius,"
   "Komar_mass,irreducible_mass,irreducible_mass_current,"
   "Christodoulou_mass,Christodoulou_mass_current,AH_area,"
   "Omega_x,Omega_y,Omega_z,chi_x,chi_y,chi_z,"
   "chi_x_current,chi_y_current,chi_z_current,"
   "spin_x,spin_y,spin_z,"
   "Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM,"
   "boost_Vx,boost_Vy,boost_Vz");
 
  /* the very first BH approximation
  // options:
  // ========
  // o. CloseKerrSchild (see BH physics) */
  Pset_default("BH_start_off","CloseKerrSchild"); 
  
  /* max l in Ylm expansion */
  Pset_default("BH_surface_Ylm_max_l","1"); 
  
  /* BH surface type */
  Pset_default("BH_surface_type","perfect_s2"); 
  
  /* equation related: */
  /* XCTS means: alpha,beta,psi */
  Pset_default("BH_Eq_inner_BC_fields","XCTS");
  
  /* set alpha on AH to be exp(-r^4)*KerrSchild value */
  Pset_default("BH_Eq_inner_BC_alpha","exp(-r^4)*KerrSchild");
  
  /* set bete^i on AH to be alpha*s^i + Omega x r.  */
  Pset_default("BH_Eq_inner_BC_beta","alpha+Omega*r");
  
  /* observe method pertinet to BH */
  Pset_default("BH_Observe_ADM_M","S_obj,default");
  Pset_default("BH_Observe_Komar_M","S_obj,default");
  Pset_default("BH_Observe_Irreducible_M","S_obj,default");
  Pset_default("BH_Observe_ADM_P","S_obj,default");
  Pset_default("BH_Observe_ADM_J","S_obj,default");
  Pset_default("BH_Observe_CM","S_obj,default");
  Pset_default("BH_Observe_spin","S_obj,Campanelli");
  
  /* how to tune BH Radius */
  Pset_default("BH_tune_BH_radius_criteria","fix_irreducible_mass");
  Pset_default("BH_mass_tolerance","1E-5");
  Pset_default("BH_radius_update_weight","0.");
  Pset_default("BH_did_BH_surface_change","1");
  
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
  
  /* modify parameters if needed */
  
  /* use Kepler's law to set angular velocity */
  if (Pcmps(P_"angular_velocity","auto"))
  {
    const double m1 = Pgetd("BH_irreducible_mass");
    const double m2 = Pgetd("NS_baryonic_mass");
    const double r  = Pgetd(P_"separation");
    const double O  = sqrt((m1+m2)/pow(r,3.));
    
    Psetd(P_"angular_velocity",O);
  }
  
  /* set BH center and NS center */
  if (strstr(Pgets("grid_set_BH"),"right"))
  {
    const double S = Pgetd(P_"separation");
    
    /* BH center in +y */
    Psetd("BH_center_x",0.);
    Psetd("BH_center_y",S/2.);
    Psetd("BH_center_z",0.);
    
    /* NS center in -y */
    Psetd("NS_center_x",0.);
    Psetd("NS_center_y",-S/2.);
    Psetd("NS_center_z",0.);
  }
  else
  {
    const double S = Pgetd(P_"separation");
    
    /* BH center in -y */
    Psetd("BH_center_x",0.);
    Psetd("BH_center_y",-S/2.);
    Psetd("BH_center_z",0.);
    
    /* NS center in +y */
    Psetd("NS_center_x",0.);
    Psetd("NS_center_y",S/2.);
    Psetd("NS_center_z",0.);
  }
  
  /* set center of mass */
  {
    const double irr_mass = Pgetd("BH_irreducible_mass");
    const double bar_mass = Pgetd("NS_baryonic_mass");
    const double y_CM = (bar_mass*Pgetd("NS_center_y") +
                       irr_mass*Pgetd("BH_center_y"))/(bar_mass+irr_mass);
   
    Psetd(P_"x_CM",0.);
    Psetd(P_"y_CM",y_CM);
    Psetd(P_"z_CM",0.);
  }
  
  /* boost speed for BH in x direction */
  {
    const double Omega = Pgetd(P_"angular_velocity");
    const double y_CM  = Pgetd(P_"y_CM");
    const double BH_center_y = Pgetd("BH_center_y");
    
    Psetd("BH_boost_Vx",-Omega*(BH_center_y-y_CM));
  }
}
