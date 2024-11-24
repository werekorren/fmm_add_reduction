#ifndef FMM_ADDITION_REDUCTION
#define FMM_ADDITION_REDUCTION
#include "fmm_algorithm_specification.h"

/* reduction methods */
typedef enum {
  reduction_method_brute_force,
  reduction_method_greedy_vanilla,
  reduction_method_greedy_potential,
} reduction_method;

typedef union {
  struct { /* no additional parameters */ } _brute_force;
  struct { /* no additional parameters */ } _vanilla;
  struct { int k1; int k2; } _potential;
} reduction_parameters;

void fmm_addition_reduction(fmm_matrix *A, reduction_method red, reduction_parameters *red_par, int verbose);

void find_best_greedy_potential_parameters(fmm_matrix *A, double alpha_start, double alpha_end, int num_steps, int *k1, int *k2, int *num_add_after_reduction, int verbose);

void fmm_reduction_move(fmm_matrix *A, int row1, int row2, int move_type);
void fmm_reduction_move_undo(fmm_matrix *A);

#endif
