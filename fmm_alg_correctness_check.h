#ifndef FMM_ALG_CORRECTNESS_CHECK_H
#define FMM_ALG_CORRECTNESS_CHECK_H
#include "fmm_algorithm_specification.h"

/* Checks if algorithm does in fact produce the expected output matrix.
 * This implementation supports arbitrary (signed) int entries in (A, B, M).
 * Assumes only that partial sums do not overflow with platform-size int type,
 * which is a non-issue for the coefficient sizes under consideration.
 */
int fmm_alg_is_correct(fmm_alg *alg, int verbose);

#endif
