#include "bhns_header.h"
#include "elliptica_id_reader_lib.h"

/* prefix parameters came from evo codes, should be lower case  */
#define BAM_ "bam_"
#define EVO_ "evo_"

#define STR_LEN_MAX (999)

void bhns_bam_exporting_initial_data(void *vp);
void bhns_evo_exporting_initial_data(Elliptica_ID_Reader_T *const idr);
void bhns_set_bam_fields(Grid_T *const grid);
void bhns_set_evo_fields(Grid_T *const grid);

