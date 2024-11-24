#include "fmm_brute_force.h"
#include "fmm_addition_reduction.h"
#include <stdio.h>
#include <assert.h>
#include <malloc.h>

static fmm_matrix best_m;
static int best_num_ops;

/* returns the number of times that a substitution of type t = A_row1 +/- A_row2 can be applied */
/* (if it can be applied t times, it can save t-1 addition operations) */
static void overlap_counter(fmm_matrix *A, int row1, int row2, int *p, int *n) {
  int pp = 0; // count number of matches of type t = A_row1 + A_row2
  int nn = 0; // count number of matches of type t = A_row1 - A_row2
  for (int i=0; i<(A->cols + A->t); i++) { // for all columns
    int e1 = fmm_matrix_entry(A, row1, i);
    if (e1 == 0) {
      continue;
    }
    int e2 = fmm_matrix_entry(A, row2, i);
    if (e2 == 0) {
      continue;
    }
    if (e1 == e2) {
      pp++;
    } else if (e1 == -e2) {
      nn++;
    }
  }
  *n = nn;
  *p = pp;
}

static void potential(fmm_matrix *m, int *s, int *num_non_zero) {
  int ss = 0;
  int nz = 0;
  int t = m->rows + m->t;
  for (int i=0; i<t; i++) {
    for (int j=i+1; j<t; j++) {
      int p, n;
      overlap_counter(m, i, j, &p, &n);
      if (p > 1) {
        ss += p - 1;
        nz++;
      }
      if (n > 1) {
        ss += n - 1;
        nz++;
      }
    }
  }
  *s = ss;
  *num_non_zero = nz;
}

static void printfn(int n, const char *s) {
  for (int i=0; i<n; i++) {
    printf(s);
  }
}

static void fmm_brute_force_matrix_internal(fmm_matrix *m, int depth, int verbose) {
  int s, tot;
  potential(m, &s, &tot);
  int a = fmm_matrix_num_additions(m);
  if (a - s >= best_num_ops) { // the number of ops we can save with all actions combined is less than we need to find a better solution
    return; // prune search tree here
  }

  int t = m->rows + m->t;
//  int tot = t * (t - 1) / 2;
  int iteration = 0;
  for (int i=0; i<t; i++) {
    for (int j=i+1; j<t; j++) {
      int p, n;
      overlap_counter(m, i, j, &p, &n);

      if (p > 1) { // substitution u = row i + row j saves more than one addition
        fmm_reduction_move(m, i, j, 1);
        int a = fmm_matrix_num_additions(m);
        if (a < best_num_ops) {
          fmm_matrix_copy(&best_m, m);
          best_num_ops = a;
          if (verbose) {
            printf("current best uses %d additions, found at depth %d\n", best_num_ops, depth);
          }
        }
        fmm_brute_force_matrix_internal(m, depth + 1, verbose);
        fmm_reduction_move_undo(m);
        iteration++;
        if (verbose && depth < verbose) {
          printfn(depth, "------------------------------- ");
          printf("[depth %d, %3d of %3d = %6.2f%%]\n", depth, iteration, tot, iteration*100/(double)tot);
        }
      }

      if (n > 1) { // substitution u = row i - row j saves more than one addition
        fmm_reduction_move(m, i, j, -1);
        int a = fmm_matrix_num_additions(m);
        if (a < best_num_ops) {
          fmm_matrix_copy(&best_m, m);
          best_num_ops = a;
          if (verbose) {
            printf("current best uses %d additions, found at depth %d\n", best_num_ops, depth);
          }
        }
        fmm_brute_force_matrix_internal(m, depth + 1, verbose);
        fmm_reduction_move_undo(m);
        iteration++;
        if (verbose && depth < verbose) {
          printfn(depth, "------------------------------- ");
          printf("[depth %d, %3d of %3d = %6.2f%%]\n", depth, iteration, tot, iteration*100/(double)tot);
        }
      }

    }
  }
}

void fmm_brute_force_matrix(fmm_matrix *m, int verbose) {
  fmm_matrix_init(&best_m, m->rows, m->cols, m->t_capacity);
  if (verbose) {
    printf("original uses %d additions\n", fmm_matrix_num_additions(m));
  }

  /* (over) estimate optimal number of additions to improve pruning */
  /* we can use greedy methods to provide a very good approximation */
  //fmm_matrix_copy(&best_m, m);
  //best_num_ops = fmm_matrix_num_additions(m); // initial approximation can be set manually like this
  fmm_matrix temp;
  fmm_matrix_init(&temp, m->rows, m->cols, m->t_capacity);
  fmm_matrix_copy(&temp, m);
#if 0 // using greedy potential to determine pruning level
  const char *greedy_method = "potential";
  reduction_parameters par;
  par._potential.k1 = 10;
  par._potential.k2 = 1;
  fmm_addition_reduction(&temp, reduction_method_greedy_potential, &par, 0 /* silent mode */);
#else // using greedy vanilla to determine pruning level
  const char *greedy_method = "vanilla";
  fmm_addition_reduction(&temp, reduction_method_greedy_vanilla, NULL, 0 /* silent mode */);
#endif
  fmm_matrix_copy(&best_m, &temp);
  best_num_ops = fmm_matrix_num_additions(&temp);
  fmm_matrix_destroy(&temp);
  if (verbose) {
    printf("initial best greedy %s solution uses %d additions (improved pruning speeds up search)\n", greedy_method, best_num_ops);
  }

  /* perform the brute-force search */
  fmm_brute_force_matrix_internal(m, 0, verbose);

  /* restore optimal result to m */
  fmm_matrix_copy(m, &best_m);
  if (verbose) {
    printf("final best uses %d additions\n", best_num_ops);
  }
  fmm_matrix_destroy(&best_m);
}

void fmm_brute_force_alg(fmm_alg *alg, int verbose) {
  fmm_brute_force_matrix(&alg->A, verbose);
  fmm_brute_force_matrix(&alg->B, verbose);
  fmm_brute_force_matrix(&alg->C, verbose);
}
