#include "bhns_header.h"

/* useful macro */
#define MAX_STR_LEN (400)
#define LINE_STR    "-------------------------------------------------------------------------"


void bhns_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen);

void bhns_analyze(Physics_T *const phys,const int iteration);


                                          
