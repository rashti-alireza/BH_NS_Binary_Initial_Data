#include "bhns_header.h"
#include "id_reader_lib.h"
#include "elliptica_id_reader_lib.h"

/* prefix parameters came from evo codes, should be lower case  */
#define BAM_ "bam_"
#define EVO_ "evo_"

#define STR_LEN_MAX (999)

void bhns_export_id_bam_generic(void *vp);
void bhns_export_id_generic(Elliptica_ID_Reader_T *const idr);
void bhns_set_bam_fields_generic(Grid_T *const grid);
void bhns_set_evo_fields_generic(Grid_T *const grid);

