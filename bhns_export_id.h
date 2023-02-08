#include "bhns_header.h"
#include "elliptica_id_reader_lib.h"

/* prefix parameters came from evo codes, should be lower case  */
#define BAM_ "bam_"
#define EVO_ "evo_"

#define STR_LEN_MAX (999)

void bhns_bam_read_id_asymptotically_inertial(void *vp);
void bhns_read_id_asymptotically_inertial(Elliptica_ID_Reader_T *const idr);
void bhns_set_bam_fields(Grid_T *const grid);
void bhns_set_evo_fields_asymptotically_inertial(Grid_T *const grid);

