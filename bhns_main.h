#include "bhns_header.h"

#define STR_LEN_MAX (400)

int BH_NS_Binary_Initial_Data(void *vp);
void bhns_export_id_bam_generic(void *vp);
void bhns_export_id_generic(void *vp);
void bhns_export_id_generic_mt_safe(void *vp);
static void construct_initial_data(void *vp);
static void set_default_parameters(void);
static void fill_blackhole(Physics_T *const phys);

