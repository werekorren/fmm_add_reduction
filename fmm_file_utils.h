#ifndef FMM_FILE_UTILS_H
#define FMM_FILE_UTILS_H
#include <stdio.h>
#include "fmm_algorithm_specification.h"

/* reads algorithm and meta data from file
 * Grey format
 * A has k * l rows and q columns
 * B has l * m rows and q columns
 * M has q rows and k * m columns
 * */
int read_from_file(fmm_alg *alg, const char *file_name, int t_capacity, int verbose);

#endif
