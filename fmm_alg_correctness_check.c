#include "fmm_alg_correctness_check.h"
#include <memory.h>
#include <stdio.h>
#include <assert.h>
#include "fmm_matrix_storage_utils.h"
#include "fmm_addition_reduction.h"

static void add_m_to_pattern(fmm_matrix *pattern, fmm_alg *alg, int mi, int factor) {
  for (int ai=0; ai<alg->na; ai++) {
    int aa = fmm_matrix_entry(&(alg->A), ai, mi);
    if (aa) {
      for (int bi=0; bi<alg->nb; bi++) {
        int bb = fmm_matrix_entry(&(alg->B), bi, mi);
        int v = fmm_matrix_entry(pattern, ai, bi) + factor * aa * bb;
        fmm_matrix_entry_set(pattern, ai, bi, v);
      }
    }
  }
}

static int e_mod_p(int e, int p) {
  if (p == 0) {
    return e; // computations in Z, so we do not reduce
  }
  assert(p > 0 && "unexpected parameter p");
  while (e < 0) {
    e += p;
  }
  return e % p;
}

/* note: modifies pattern */
static int pattern_matches_ci(fmm_matrix *pattern, fmm_alg *alg, int ci, int verbose) {
  /* figure out which row in A and col in B that correspond to entry ci in C */
  int row = ci / alg->m;
  int col = ci % alg->m;

  /* check each of the l terms in the ab-sum corresponding to ci */
  for (int i=0; i<alg->l; i++) {
    int ai = alg->l * row + i;
    int bi = col + i * alg->m;
    int e = fmm_matrix_entry(pattern, ai, bi);
    int emp = e_mod_p(e, alg->p); // reduce mod p if Z_p (does not reduce for Z)
    if (emp != 1) {
      if (verbose) {
        printf("*** checking p[%d][%d] = %d [INCORRECT - should be 1] ***\n", ai, bi, e);
      }
      return 0;
    }
    if (verbose) {
      printf("checking p[%d][%d] = 1 [correct]\n", ai, bi);
    }

    /* set to zero (modifies pattern) */
    fmm_matrix_entry_set(pattern, ai, bi, 0);
  }

  /* check that all remaining entries are zero now */
  for (int ai=0; ai<alg->na; ai++) {
    for (int bi=0; bi<alg->nb; bi++) {
      int e = fmm_matrix_entry(pattern, ai, bi);
      int emp = e_mod_p(e, alg->p); // reduce mod p if Z_p (does not reduce for Z)
      if (emp) {
        if (verbose) {
          printf("*** checking p[%d][%d] = %d [INCORRECT - should be 0] ***\n", ai, bi, e);
        }
        return 0;
      }
      if (verbose) {
        printf("checking p[%d][%d] = 0 [correct]\n", ai, bi);
      }
    }
  }

  /* all checks ok (and pattern zeroized as a side-effect) */
  return 1;
}

int fmm_alg_is_correct(fmm_alg *alg, int verbose) {
  fmm_matrix pattern;

  /* pattern for ab terms */
  fmm_matrix_init(&pattern, alg->na, alg->nb, 0 /* t_capacity */);
  fmm_matrix_set_all_entries(&pattern, 0); // zero-initialize pattern

  /* copy alg */
  fmm_alg alg_copy;
  fmm_alg_deep_copy(&alg_copy, alg);

  /* back-substitute t-space */
  while (alg_copy.A.t) {
    fmm_reduction_move_undo(&alg_copy.A);
  }
  while (alg_copy.B.t) {
    fmm_reduction_move_undo(&alg_copy.B);
  }
  while (alg_copy.C.t) {
    fmm_reduction_move_undo(&alg_copy.C);
  }

  /* verify ab-patterns for each c_i separately */
  for (int ci=0; ci<alg_copy.nc; ci++) {

    /* overlay the M-patterns that correspond to entry c_i in the C-matrix */
    for (int mi=0; mi<alg_copy.q; mi++) {
      if (fmm_matrix_entry(&alg_copy.C, mi, ci)) {
        int factor = fmm_matrix_entry(&alg_copy.C, mi, ci);
        add_m_to_pattern(&pattern, &alg_copy, mi, factor);
      }
    }

    if (verbose) {
      /* aux test printing */
      printf("pattern at ci = %d\n", ci);
      fmm_matrix_print(&pattern, " %d");
    }

    /* check pattern for ci */
    if (!pattern_matches_ci(&pattern, &alg_copy, ci, verbose)) {
      fmm_matrix_destroy(&pattern);
      return 0;
    }

    /* note: pattern is implicitly zero-initialized in pattern_matches_ci */
  }

  /* cleanup */
  fmm_alg_destroy(&alg_copy);
  fmm_matrix_destroy(&pattern);
  return 1;
}
















#if 0
static void add_m_to_pattern(fmm_matrix *pattern, fmm_alg *alg, int mi, int factor) {
  for (int ai=0; ai<alg->na; ai++) {
    int aa = fmm_matrix_entry(&(alg->A), ai, mi);
    if (aa) {
      for (int bi=0; bi<alg->nb; bi++) {
        int bb = fmm_matrix_entry(&(alg->B), bi, mi);
        int v = fmm_matrix_entry(pattern, ai, bi) + factor * aa * bb;
        fmm_matrix_entry_set(pattern, ai, bi, v);
      }
    }
  }
}

static int e_mod_p(int e, int p) {
  if (p == 0) {
    return e; // computations in Z, so we do not reduce
  }
  assert(p > 0 && "unexpected parameter p");
  while (e < 0) {
    e += p;
  }
  return e % p;
}

/* note: modifies pattern */
static int pattern_matches_ci(fmm_matrix *pattern, fmm_alg *alg, int ci, int verbose) {
  /* figure out which row in A and col in B that correspond to entry ci in C */
  int row = ci / alg->m;
  int col = ci % alg->m;

  /* check each of the l terms in the ab-sum corresponding to ci */
  for (int i=0; i<alg->l; i++) {
    int ai = alg->l * row + i;
    int bi = col + i * alg->m;
    int e = fmm_matrix_entry(pattern, ai, bi);
    int emp = e_mod_p(e, alg->p); // reduce mod p if Z_p (does not reduce for Z)
    if (emp != 1) {
      if (verbose) {
        printf("*** checking p[%d][%d] = %d [INCORRECT - should be 1] ***\n", ai, bi, e);
      }
      return 0;
    }
    if (verbose) {
      printf("checking p[%d][%d] = 1 [correct]\n", ai, bi);
    }

    /* set to zero (modifies pattern) */
    fmm_matrix_entry_set(pattern, ai, bi, 0);
  }

  /* check that all remaining entries are zero now */
  for (int ai=0; ai<alg->na; ai++) {
    for (int bi=0; bi<alg->nb; bi++) {
      int e = fmm_matrix_entry(pattern, ai, bi);
      int emp = e_mod_p(e, alg->p); // reduce mod p if Z_p (does not reduce for Z)
      if (emp) {
        if (verbose) {
          printf("*** checking p[%d][%d] = %d [INCORRECT - should be 0] ***\n", ai, bi, e);
        }
        return 0;
      }
      if (verbose) {
        printf("checking p[%d][%d] = 0 [correct]\n", ai, bi);
      }
    }
  }

  /* all checks ok (and pattern zeroized as a side-effect) */
  return 1;
}

/* implementation assumes that  */
int fmm_alg_is_correct(fmm_alg *alg, int verbose) {
  fmm_matrix pattern;

  fmm_matrix_init(&pattern, alg->na, alg->nb, 0 /* t_capacity */);
  fmm_matrix_set_all_entries(&pattern, 0); // zero-initialize pattern

  for (int ci=0; ci<alg->nc; ci++) {

    /* overlay the M-patterns that correspond to entry c_i in the C-matrix */
    for (int mi=0; mi<alg->q; mi++) {
      if (fmm_matrix_entry(&(alg->C), mi, ci)) {
        int factor = fmm_matrix_entry(&(alg->C), mi, ci);
        add_m_to_pattern(&pattern, alg, mi, factor);
      }
    }

    if (verbose) {
      /* aux test printing */
      printf("pattern at ci = %d\n", ci);
      fmm_matrix_print(&pattern, " %d");
    }

    /* check pattern for ci */
    if (!pattern_matches_ci(&pattern, alg, ci, verbose)) {
      fmm_matrix_destroy(&pattern);
      return 0;
    }

    /* note: pattern is implicitly zero-initialized in pattern_matches_ci */
  }

  /* cleanup */
  fmm_matrix_destroy(&pattern);
  return 1;
}

#endif
