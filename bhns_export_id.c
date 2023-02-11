/*
// Alireza Rashti
// January 2021
*/

/* exporting initial data for evolution codes */

#include "bhns_export_id.h"


/* exporting initial data for bam.
// it writes the required fields into a file to be read by bam. */
void bhns_export_id_bam_generic(void *vp)
{
  FUNC_TIC
  
  Physics_T *bhns    = 0;
  ID_Reader_T *points = idr_init();
  FILE *file          = 0;
  char fields_name[STR_LEN_MAX] = {'\0'};
  char **sfield = 0;
  Uint f;
  
  /* don't stop */
  Pseti(CHECKPOINT_SET_PARAM_ P_"STOP",0);
  
  /* go from Omega x r to inertial coords sys asymptotically.
  // BAM needs this!? (it will be set in the reader) */
  // Psets(CHECKPOINT_SET_PARAM_ "ADM_B1I_form","zero");
  
  /* read physics from checkpoint */
  Psets("checkpoint_file_path",Pgets(P_ BAM_"checkpoint_file_path"));
  bhns = bhns_read_physics_from_checkpoint();
  points->grid = bhns->grid;
  
  physics(bhns,ADM_UPDATE_Kij);/* before filling */
  /* fill BH */
  Physics_T *const bh  = init_physics(bhns,BH);
  Psets("BH_filler_method",Pgets(P_ BAM_"filler_method"));
  Pseti("BH_filler_verbose",1);/* make it verbose anyway. */
  /* fill these fields */
  Psets("BH_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                           "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                           "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                           "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                           "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  physics(bh,BH_FILL);
  free_physics(bh);
    
  /* set bam fields based on initial data to be usable for bam */
  bhns_set_bam_fields_generic(bhns->grid);
 
  /* read (x,y,z) points from bam file to be interpolated on them */
  idr_load_Cartesian_coordinates_from_file
    (Pgets(P_ BAM_"coords_file_path"),points);
  
  /* open a binary file to write fields in it. */
  file = idr_new_binary_file_to_write
    (Pgets(P_ BAM_"fields_file_path"),Pgets(P_ BAM_"fields_name"));
  
  /* adapt fields_notations for Elliptica */
  assert(sprintf(fields_name,"%s",Pgets(P_ BAM_"fields_name")));
  
  /* metric fields */
  regex_replace(fields_name,"\\balpha\\b",BAM_"alpha",fields_name);
  
  regex_replace(fields_name,"\\bbetax\\b",BAM_"beta_U0",fields_name);
  regex_replace(fields_name,"\\bbetay\\b",BAM_"beta_U1",fields_name);
  regex_replace(fields_name,"\\bbetaz\\b",BAM_"beta_U2",fields_name);
  
  regex_replace(fields_name,"\\badm_gxx\\b",BAM_"adm_g_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_gxy\\b",BAM_"adm_g_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_gxz\\b",BAM_"adm_g_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_gyy\\b",BAM_"adm_g_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_gyz\\b",BAM_"adm_g_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_gzz\\b",BAM_"adm_g_D2D2",fields_name);
  
  regex_replace(fields_name,"\\badm_Kxx\\b",BAM_"adm_Kij_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_Kxy\\b",BAM_"adm_Kij_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_Kxz\\b",BAM_"adm_Kij_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_Kyy\\b",BAM_"adm_Kij_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_Kyz\\b",BAM_"adm_Kij_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_Kzz\\b",BAM_"adm_Kij_D2D2",fields_name);
  
  /* matter fields */
  regex_replace(fields_name,"\\bgrhd_vx\\b",BAM_"grhd_v_U0",fields_name);
  regex_replace(fields_name,"\\bgrhd_vy\\b",BAM_"grhd_v_U1",fields_name);
  regex_replace(fields_name,"\\bgrhd_vz\\b",BAM_"grhd_v_U2",fields_name);
  regex_replace(fields_name,"\\bgrhd_rho\\b",BAM_"grhd_rho",fields_name);
  regex_replace(fields_name,"\\bgrhd_p\\b",BAM_"grhd_p",fields_name);
  regex_replace(fields_name,"\\bgrhd_epsl\\b",BAM_"grhd_epsl",fields_name);
  
  /* check if all fields are expected */
  sfield = read_separated_items_in_string(fields_name,',');
  f = 0;
  while(sfield[f])
  {
    if (!regex_search("^"BAM_,sfield[f]))
    {
      printf("field '%s' is Unexpected!\n",sfield[f]);
      Error1("Unexpected field! Add me!");
    }
    f++;
  }
  free_2d(sfield);
  
  /* write into file */
  idr_interpolate_fields_and_write_to_file
    (file,points,fields_name,Pgets(P_ BAM_"fields_name"));
  
  /* finishing up */
  idr_close_file(file);
  idr_free(points);
  free_physics(bhns);
  
  UNUSED(vp);
  FUNC_TOC  
}

/* export ID for a general evolution code */
void bhns_export_id_generic(Elliptica_ID_Reader_T *const idr)
{
  FUNC_TIC
  
  Physics_T *bhns = 0;
  ID_Reader_T *points = idr_init();
  double CM[3] = {0.};
  char fields_name[STR_LEN_MAX] = {'\0'};
  char **sfield = 0;
  Uint f;
  
  /* don't stop */
  Pseti(CHECKPOINT_SET_PARAM_ P_"STOP",0);
  
  /* set all parameters set in evo code. */
  // NOTE: adding both plain and modify version to make sure things kick in 
  // after loading from checkpoint file.
  for (Uint n = 0; n < idr->nparams; ++n)
  {
    char lv[STR_LEN_MAX];
    sprintf(lv,CHECKPOINT_SET_PARAM_"%s",idr->params_lv[n]);
    
    Psets(idr->params_lv[n],idr->params_rv[n]);
    Psets(lv,idr->params_rv[n]);
  }
  
  /* read physics from checkpoint */
  Psets("checkpoint_file_path",idr->checkpoint_path);
  bhns = bhns_read_physics_from_checkpoint();
  points->grid = bhns->grid;
  
  /* go from Omega x r to inertial coords sys asymptotically.
  // EVO needs this!? (it will be set in the reader) */
  // Psets(CHECKPOINT_SET_PARAM_ "ADM_B1I_form","zero");
  physics(bhns,ADM_UPDATE_Kij);/* before filling */
  
  /* fill BH */
  Physics_T *const bh  = init_physics(bhns,BH);
  Pseti("BH_filler_verbose",1);/* make it verbose anyway. */
  /* fill these fields */
  Psets("BH_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                           "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                           "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                           "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                           "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  physics(bh,BH_FILL);
  free_physics(bh);
    
  /* set bam fields based on initial data to be usable for evo */
  bhns_set_evo_fields_generic(bhns->grid);
 
  /* get (x,y,z) points from evo. NOTE: no allocation done for (x,y,z) */
  CM[0] = Pgetd(P_"x_CM");
  CM[1] = Pgetd(P_"y_CM");
  CM[2] = Pgetd(P_"z_CM");
  idr_find_XYZ_from_xyz(idr,points,CM);
  
  /* adapt fields_notations for Elliptica */
  assert(sprintf(fields_name,"%s",idr->ifields));
  
  /* metric fields */
  regex_replace(fields_name,"\\balpha\\b",EVO_"alpha",fields_name);
  
  regex_replace(fields_name,"\\bbetax\\b",EVO_"beta_U0",fields_name);
  regex_replace(fields_name,"\\bbetay\\b",EVO_"beta_U1",fields_name);
  regex_replace(fields_name,"\\bbetaz\\b",EVO_"beta_U2",fields_name);
  
  regex_replace(fields_name,"\\badm_gxx\\b",EVO_"adm_g_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_gxy\\b",EVO_"adm_g_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_gxz\\b",EVO_"adm_g_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_gyy\\b",EVO_"adm_g_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_gyz\\b",EVO_"adm_g_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_gzz\\b",EVO_"adm_g_D2D2",fields_name);
  
  regex_replace(fields_name,"\\badm_Kxx\\b",EVO_"adm_Kij_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_Kxy\\b",EVO_"adm_Kij_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_Kxz\\b",EVO_"adm_Kij_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_Kyy\\b",EVO_"adm_Kij_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_Kyz\\b",EVO_"adm_Kij_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_Kzz\\b",EVO_"adm_Kij_D2D2",fields_name);
  
  /* matter fields */
  regex_replace(fields_name,"\\bgrhd_vx\\b",EVO_"grhd_v_U0",fields_name);
  regex_replace(fields_name,"\\bgrhd_vy\\b",EVO_"grhd_v_U1",fields_name);
  regex_replace(fields_name,"\\bgrhd_vz\\b",EVO_"grhd_v_U2",fields_name);
  regex_replace(fields_name,"\\bgrhd_rho\\b",EVO_"grhd_rho",fields_name);
  regex_replace(fields_name,"\\bgrhd_p\\b",EVO_"grhd_p",fields_name);
  regex_replace(fields_name,"\\bgrhd_epsl\\b",EVO_"grhd_epsl",fields_name);
  
  /* check if all fields are expected */
  sfield = read_separated_items_in_string(fields_name,',');
  f = 0;
  while(sfield[f])
  {
    if (!regex_search("^"EVO_,sfield[f]))
    {
      printf("field '%s' is Unexpected!\n",sfield[f]);
      Error1("Unexpected field! Add me!");
    }
    f++;
  }
  free_2d(sfield);
  
  /* write into array */
  idr_interpolate_fields_and_save_in_array(idr,points,fields_name,idr->ifields);
  
  /* finishing up */
  // since no alocation done for (x,y,z):
  points->x = 0;
  points->y = 0;
  points->z = 0;
  idr_free(points);
  free_physics(bhns);
  
  FUNC_TOC  
}

