#include "sbh_header.h"

/* useful macro */
#define MAX_STR_LEN (400)
#define LINE_STR    "-------------------------------------------------------------------------"


void sbh_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen);

void sbh_analyze(Physics_T *const phys,const int iteration);


                                          
