#ifndef FMM_ALGORITHM_SPECIFICATION_H
#define FMM_ALGORITHM_SPECIFICATION_H
#include "fmm_matrix_storage_utils.h"

typedef struct {
  int p; // p=0 for Z, p>0 for Z_p
  int k, l, m;
  int q;
  int na, nb, nc;
  char name[256];
  fmm_matrix A;
  fmm_matrix B;
  fmm_matrix C;
} fmm_alg;

void fmm_alg_init(fmm_alg *alg, int p, int k, int l, int m, int q, int t_capacity);
void fmm_alg_destroy(fmm_alg *alg);
void fmm_alg_deep_copy(fmm_alg *dst, fmm_alg *src);

void fmm_alg_print(fmm_alg *alg);
void fmm_alg_print_t(fmm_alg *alg);
void fmm_alg_print_full_capacity(fmm_alg *alg);
void fmm_alg_print_latex(fmm_alg *alg, const char *a, const char *b, const char *m, const char *c, const char *t, const char *u, const char *v);

#endif
