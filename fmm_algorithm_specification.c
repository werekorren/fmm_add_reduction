#include "fmm_algorithm_specification.h"
#include "fmm_matrix_storage_utils.h"
#include <stdio.h>
#include <assert.h>

void fmm_alg_init(fmm_alg *alg, int p, int k, int l, int m, int q, int t_capacity) {
  alg->p = p;
  alg->k = k;
  alg->l = l;
  alg->m = m;
  alg->q = q;
  alg->na = k * l; // superflouous utility value
  alg->nb = l * m; // superflouous utility value
  alg->nc = k * m; // superflouous utility value
  fmm_matrix_init(&(alg->A), alg->na, q, t_capacity); // matrix of dimension (k*l) x q
  fmm_matrix_init(&(alg->B), alg->nb, q, t_capacity); // matrix of dimension (l*m) x q
  fmm_matrix_init(&(alg->C), q, alg->nc, t_capacity); // matrix of dimension q x (k*m)
  /* note: matrix entries not initialized */
}

void fmm_alg_destroy(fmm_alg *alg) {
  fmm_matrix_destroy(&(alg->A));
  fmm_matrix_destroy(&(alg->B));
  fmm_matrix_destroy(&(alg->C));
}

void fmm_alg_deep_copy(fmm_alg *dst, fmm_alg *src) {
  int max_t_capacity = src->A.t_capacity;
  if (src->B.t_capacity > max_t_capacity) {
    max_t_capacity = src->B.t_capacity;
  }
  if (src->C.t_capacity > max_t_capacity) {
    max_t_capacity = src->C.t_capacity;
  }
  fmm_alg_init(dst, src->p, src->k, src->l, src->m, src->q, max_t_capacity);
  fmm_matrix_copy(&dst->A, &src->A);
  fmm_matrix_copy(&dst->B, &src->B);
  fmm_matrix_copy(&dst->C, &src->C);
}

void fmm_alg_print(fmm_alg *alg) {
  printf("(k, l, m, q) = (%d, %d, %d, %d)\n\n", alg->k, alg->l, alg->m, alg->q);
  printf("A =\n");
  fmm_matrix_print(&(alg->A), " %2d");
  printf("\nB =\n");
  fmm_matrix_print(&(alg->B), " %2d");
  printf("\nM =\n");
  fmm_matrix_print(&(alg->C), " %2d");
}

void fmm_alg_print_t(fmm_alg *alg) {
  printf("(k, l, m, q) = (%d, %d, %d, %d)\n\n", alg->k, alg->l, alg->m, alg->q);
  printf("A =\n");
  fmm_matrix_print_t(&(alg->A), " %2d");
  printf("\nB =\n");
  fmm_matrix_print_t(&(alg->B), " %2d");
  printf("\nM =\n");
  fmm_matrix_print_t(&(alg->C), " %2d");
}

void fmm_alg_print_full_capacity(fmm_alg *alg) {
  printf("(k, l, m, q) = (%d, %d, %d, %d)\n\n", alg->k, alg->l, alg->m, alg->q);
  printf("A =\n");
  fmm_matrix_print_full_capacity(&(alg->A), " %2d");
  printf("\nB =\n");
  fmm_matrix_print_full_capacity(&(alg->B), " %2d");
  printf("\nM =\n");
  fmm_matrix_print_full_capacity(&(alg->C), " %2d");
}

static int first_nonzero_col_entry_index(fmm_matrix *A, int col) {
  for (int i=0; i<(A->rows + A->t); i++) { // for every row
    int e = fmm_matrix_entry(A, i, col);
    if (e) {
      return i;
    }
  }
  assert(0 && "first_nonzero_col_entry_index: unexpected input, function assumes column contains precisely two nonzero entries");
  return -1;
}

static void print_term(int coeff, const char *term_symbol, int index, int suppress_plus) {
  if (!coeff) {
    return;
  }
  if (!(suppress_plus && coeff > 0)) {
    printf("%s", coeff < 0 ? "-" : "+");
  }// else { printf(" "); }
  int c = coeff < 0 ? -coeff : coeff;
  if (c > 1) {
    printf("%d", c);
  }
  printf("%s_{%d}", term_symbol, index);
}

static void print_t_col(fmm_matrix *matrix, int t_col, const char *a, const char *t) {
  int ts = matrix->t; // size of t-space
  int mc = matrix->cols; // size of matrix-space col
  int mr = matrix->rows; // size of matrix-space row
  int a1 = first_nonzero_col_entry_index(matrix, mc + t_col);
  int e1 = fmm_matrix_entry(matrix, a1, mc + t_col);
  printf("%s_{%d} &=& ", t, t_col);
  print_term(e1, a1 < mr ? a : t, a1 < mr ? a1 : a1 - mr, 1);
  for (int j=a1+1; j<(matrix->rows + ts); j++) { // for all remaining rows
    e1 = fmm_matrix_entry(matrix, j, mc + t_col);
    if (e1) {
      print_term(e1, j < mr ? a : t, j < mr ? j : j - mr, 0);
    }
  }
  printf("\\\\\n");
}

static void print_m_col(fmm_matrix *matrix, int col, const char *a, const char *t, int use_parenthesis) {
  int w = fmm_matrix_column_weight(matrix, col);
  if (use_parenthesis && w > 1) {
    printf("(");
  }
  int ts = matrix->t; // size of t-space
  int mr = matrix->rows; // size of matrix-space row
  int a1 = first_nonzero_col_entry_index(matrix, col);
  int e1 = fmm_matrix_entry(matrix, a1, col);
  print_term(e1, a1 < mr ? a : t, a1 < mr ? a1 : a1 - mr, 1);
  for (int j=a1+1; j<(matrix->rows + ts); j++) { // for all remaining rows
    e1 = fmm_matrix_entry(matrix, j, col);
    if (e1) {
      print_term(e1, j < mr ? a : t, j < mr ? j : j - mr, 0);
    }
  }
  if (use_parenthesis && w > 1) {
    printf(")");
  }
}

static void matrix_print_explicit_t_substitutions(fmm_matrix *matrix, const char *a, const char *t) {
  int ts = matrix->t; // size of t-space
  for (int i=0; i<ts; i++) { // for all columns in t-space
    print_t_col(matrix, i, a, t);
  }
}

void fmm_alg_print_latex(fmm_alg *alg, const char *a, const char *b, const char *m, const char *c, const char *t, const char *u, const char *v) {
  printf("\\begin{eqnarray*}\n");
  /* print substitutions */
  matrix_print_explicit_t_substitutions(&alg->A, a, t);
  matrix_print_explicit_t_substitutions(&alg->B, b, u);
  matrix_print_explicit_t_substitutions(&alg->C, m, v);
  /* print M expressions */
  assert(alg->A.cols == alg->C.rows && "Dimension mismatch between A and C when printing");
  assert(alg->B.cols == alg->C.rows && "Dimension mismatch between A and C when printing");
  for (int i=0; i<alg->A.cols; i++) {
    printf("%s_{%d} &=& ", m, i);
    print_m_col(&alg->A, i, a, t, 1);
    printf("\\times{}");
    print_m_col(&alg->B, i, b, u, 1);
    printf("\\\\\n");
  }
  /* print C expressions */
  for (int i=0; i<alg->C.cols; i++) {
    printf("%s_{%d} &=& ", c, i);
    print_m_col(&alg->C, i, m, v, 0);
    printf("\\\\\n");
  }
  printf("\\end{eqnarray*}\n");
}
