/*
// Alireza Rashti
// January 2021
*/

/* analyzing initial data such as mass, momentum, constraints etc.  */

#include "bhns_analyze.h"

/* analyzing physics properties, constraints etc */
void bhns_analyze(Physics_T *const phys,const int iteration)
{
  if (!phys) return;

  FUNC_TIC
  
  const char *const properties_file_name = P_"properties.txt";
  FILE *file = 0;
  char str[MAX_STR_LEN];
   
  /* compute properties and constraints */ 
  physics(phys,ADM_COMPUTE_CONSTRAINTS);
  
  /* compute various properties */
  compute_properties(phys);
  
  /* open properties file in "my_directory" and save */
  sprintf(str,"%s/%s",Pgets(P_"my_directory"),properties_file_name);
  file = Fopen(str,"w");
  bhns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);

  /* open properties file in "Diagnostics" and save */
  sprintf(str,"%s/%s",Pgets(P_"Diagnostics"),properties_file_name);
  file = Fopen(str,"a");
  bhns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);
  
  /* prints */
  print_fields_0D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_1D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_2D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_3D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  
  FUNC_TOC
}

/* print physical system properties such as mass, spin etc in the given
// file, if pr_screen is 1, it also prints in stdout */
void bhns_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen)
{
  Physics_T *const bh = init_physics(phys,BH);
  Physics_T *const ns = init_physics(phys,NS);

  if (pr_screen)
  {
    printf(Pretty0"iteration = %d:\n",iteration);
  }
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %d\n",iteration);
  fprintf(file,"\n");
  
  bh_print_properties(bh,Pgets(P_"BH_properties"),file,pr_screen);
  star_print_properties(ns,Pgets(P_"NS_properties"),file,pr_screen);
  sys_print_properties(phys,Pgets(P_"BHNS_properties"),file,pr_screen);
  
  free_physics(bh);
  free_physics(ns);
}

/* compute variety of properties.
// NOTE: order of parameter calculations matter. 
// NOTE: if there is a confusion between target params and current
//       params, "current" suffix added to the latter. */
static void compute_properties(Physics_T *const phys/* bhns */)
{
  Physics_T *const ns = init_physics(phys,NS);
  Physics_T *const bh = init_physics(phys,BH);
  const double x_CM   = Pgetd(P_"x_CM");
  const double y_CM   = Pgetd(P_"y_CM");
  const double z_CM   = Pgetd(P_"z_CM");
  double im[2] = {0.};
  double p[3]  = {0.};
  double j[3]  = {0.};
  double s[3]  = {0.};
  double cm[3] = {0.};
  double m     = 0.;
  
  /* NS: */
  observe(ns,"ADM(M)",Pgets("NS_Observe_ADM_M"),&m);
  Psetd("NS_ADM_mass",m);
  
  observe(ns,"Komar(M)",Pgets("NS_Observe_Komar_M"),&m);
  Psetd("NS_Komar_mass",m);
  
  observe(ns,"Baryonic(M)",Pgets("NS_Observe_baryonic_M"),&m);
  Psetd("NS_baryonic_mass_current",m);
  
  TOV_T *tov = TOV_init();
  tov->exit_if_error = 0;
  tov->phys  = ns;
  tov->bar_m = Pgetd("NS_baryonic_mass_current");
  tov = TOV_solution(tov);
  if (tov->status == 0)
  {
    Psetd("NS_TOV_ADM_mass",tov->ADM_m);
    /* Note: compactness = adm_mass/Schwarzschild_radius 
      (not isotropic radius) */
    Psetd("NS_TOV_compactness",tov->ADM_m/tov->r[tov->N-1]);
    Psetd("NS_TOV_radius",tov->rbar[tov->N-1]);
  }
  TOV_free(tov);
  
  Psetd("NS_mass_shedding_indicator",star_NS_mass_shedding_indicator(ns));
  
  observe(ns,"CM",Pgets("NS_Observe_CM"),cm);
  Psetd("NS_x_CM",cm[0]+x_CM);
  Psetd("NS_y_CM",cm[1]+y_CM);
  Psetd("NS_z_CM",cm[2]+z_CM);

  observe(ns,"ADM(P)",Pgets("NS_Observe_ADM_P"),p);
  Psetd("NS_Px_ADM",p[0]);
  Psetd("NS_Py_ADM",p[1]);
  Psetd("NS_Pz_ADM",p[2]);

  observe(ns,"ADM(J)",Pgets("NS_Observe_ADM_J"),j);
  Psetd("NS_Jx_ADM",j[0]);
  Psetd("NS_Jy_ADM",j[1]);
  Psetd("NS_Jz_ADM",j[2]);
  
  observe(ns,"spin",Pgets("NS_Observe_spin"),s);
  Psetd("NS_Spin_x",s[0]);
  Psetd("NS_Spin_y",s[1]);
  Psetd("NS_Spin_z",s[2]);
  
  m = Pgetd("NS_adm_mass");
  Psetd("NS_chi_x",s[0]/Pow2(m));
  Psetd("NS_chi_y",s[1]/Pow2(m));
  Psetd("NS_chi_z",s[2]/Pow2(m));
  
  /* BH: */
  observe(bh,"Komar(M)",Pgets("BH_Observe_Komar_M"),&m);
  Psetd("BH_Komar_mass",m);
  
  observe(bh,"Irreducible(M)",Pgets("BH_Observe_irreducible_M"),im);
  Psetd("BH_irreducible_mass_current",im[0]);
  Psetd("BH_AH_area",im[1]);
  
  observe(bh,"CM",Pgets("BH_Observe_CM"),cm);
  Psetd("BH_x_CM",cm[0]+x_CM);
  Psetd("BH_y_CM",cm[1]+y_CM);
  Psetd("BH_z_CM",cm[2]+z_CM);

  observe(bh,"ADM(P)",Pgets("BH_Observe_ADM_P"),p);
  Psetd("BH_Px_ADM",p[0]);
  Psetd("BH_Py_ADM",p[1]);
  Psetd("BH_Pz_ADM",p[2]);

  observe(bh,"ADM(J)",Pgets("BH_Observe_ADM_J"),j);
  Psetd("BH_Jx_ADM",j[0]);
  Psetd("BH_Jy_ADM",j[1]);
  Psetd("BH_Jz_ADM",j[2]);
  
  observe(bh,"spin",Pgets("BH_Observe_spin"),s);
  Psetd("BH_Spin_x",s[0]);
  Psetd("BH_Spin_y",s[1]);
  Psetd("BH_Spin_z",s[2]);
  
  /* calculate BH current Christodoulou mass.
  // NOTE: s[?] depend on BH spin on top */
  double irr_mass = Pgetd("BH_irreducible_mass_current");
  double net_spin = sqrt(Pow2(s[0])+Pow2(s[1])+Pow2(s[2]));
  m = sqrt(Pow2(irr_mass)+Pow2(net_spin)/(4*Pow2(irr_mass)));
  Psetd("BH_Christodoulou_mass_current",m);
  
  Psetd("BH_chi_x_current",s[0]/Pow2(m));
  Psetd("BH_chi_y_current",s[1]/Pow2(m));
  Psetd("BH_chi_z_current",s[2]/Pow2(m));
  
  /* BHNS: */
  observe(phys,"ADM(M)",Pgets(P_"Observe_ADM_M"),&m);
  Psetd(P_"adm_mass",m);
  
  observe(phys,"Komar(M)",Pgets(P_"Observe_Komar_M"),&m);
  Psetd(P_"Komar_mass",m);
  
  observe(phys,"ADM(P)",Pgets(P_"Observe_ADM_P"),p);
  Psetd(P_"Px_ADM",p[0]);
  Psetd(P_"Py_ADM",p[1]);
  Psetd(P_"Pz_ADM",p[2]);

  observe(phys,"ADM(J)",Pgets(P_"Observe_ADM_J"),j);
  Psetd(P_"Jx_ADM",j[0]);
  Psetd(P_"Jy_ADM",j[1]);
  Psetd(P_"Jz_ADM",j[2]);

  /* mass ratio */
  double q = 
    Pgetd("BH_Christodoulou_mass_current")/Pgetd("NS_TOV_ADM_mass");
  Psetd(P_"mass_ratio",q);

  /* binding energy */
  double bin_e = Pgetd(P_"adm_mass")-
    (Pgetd("BH_Christodoulou_mass_current")+Pgetd("NS_TOV_ADM_mass"));
  Psetd(P_"binding_energy",bin_e);

  /* virial error */
  double v_e = fabs(1.-Pgetd(P_"adm_mass")/Pgetd(P_"komar_mass"));
  Psetd(P_"virial_error",v_e);
  
  /* number of orbits (lowest order PN) */
  if (0)/* too much inaccurate */
  {
    double m1    = Pgetd("NS_adm_mass");
    double m2    = Pgetd("BH_Christodoulou_mass_current");
    double nu    = m1*m2/Pow2(m1+m2);
    double omega = Pgetd(P_"angular_velocity");
    double m_tot = Pgetd(P_"ADM_mass");
    double N_orb = pow(m_tot*omega,-5./3.)/(32.*nu)/(2.*M_PI);
    /* initially some vars might be off or negative like total adm_mass */
    N_orb = (isfinite(N_orb) && N_orb > 0. ? N_orb : 0.);
    Psetd(P_"number_of_orbits_1PN",N_orb);
  }
  
  free_physics(ns);
  free_physics(bh);
}
