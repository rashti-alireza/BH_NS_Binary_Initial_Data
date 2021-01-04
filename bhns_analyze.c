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
  FUNC_TIC
  
  Physics_T *const bh = init_physics(phys,BH);

  if (pr_screen)
  {
    printf(Pretty0"iteration = %d:\n",iteration);
  }
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %d\n",iteration);
  fprintf(file,"\n");
  
  bh_print_properties(bh,Pgets(P_"BH_properties"),file,pr_screen);
  sys_print_properties(phys,Pgets(P_"BHNS_properties"),file,pr_screen);
  
  free_physics(bh);
  
  FUNC_TOC
}
