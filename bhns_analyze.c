/*
// Alireza Rashti
// January 2021
*/

/* analyzing initial data such as mass, momentum, constraints etc.  */

#include "bhns_analyze.h"

/* analyzing physics properties, constraints etc */
void bhns_analyze(Physics_T *const phys,const int iteration)
{
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
  bhns_print_physical_system_properties(phys,file,iteration,1);
  Fclose(file);

  /* open properties file in "Diagnostics" and save */
  sprintf(str,"%s/%s",Pgets(P_"Diagnostics"),properties_file_name);
  file = Fopen(str,"a");
  bhns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);
  
  /* prints */
  print_fields_3D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_1D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  
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
// NOTE: order of parameters matter. */
static void compute_properties(Physics_T *const phys/* bhns */)
{
  Physics_T *const ns = init_physics(phys,NS);
  double m     = 0.;
  double p[3]  = {0};
  double j[3]  = {0};
  double s[3]  = {0};
  double cm[3] = {0};
  
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

  
  /* NS: */
  observe(ns,"ADM(M)",Pgets("NS_Observe_ADM_M"),&m);
  Psetd("NS_ADM_mass",m);
  
  observe(ns,"Komar(M)",Pgets("NS_Observe_Komar_M"),&m);
  Psetd("NS_Komar_mass",m);
  
  observe(ns,"Baryonic(M)",Pgets("NS_Observe_baryonic_M"),&m);
  Psetd("NS_baryonic_mass_current",m);
  
  TOV_T *tov = TOV_init();
  tov->phys  = ns;
  tov->bar_m = Pgetd("NS_baryonic_mass_current");
  tov = TOV_solution(tov);
  Psetd("NS_TOV_ADM_mass",tov->ADM_m);
  Psetd("NS_TOV_compactness",tov->ADM_m/tov->rbar[tov->N-1]);
  TOV_free(tov);
  
  Psetd("NS_shedding_indicator",star_NS_mass_shedding_indicator(ns));
  
  observe(ns,"CM",Pgets("NS_Observe_CM"),cm);
  Psetd("NS_x_CM",cm[0]);
  Psetd("NS_y_CM",cm[1]);
  Psetd("NS_z_CM",cm[2]);

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
  
  free_physics(ns);
}
