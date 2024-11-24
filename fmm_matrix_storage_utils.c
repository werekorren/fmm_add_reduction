#include "fmm_matrix_storage_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <memory.h>

void fmm_matrix_init(fmm_matrix *matrix, int rows, int cols, int t_capacity) {
  matrix->rows = rows;
  matrix->cols = cols;
  matrix->t = 0;
  matrix->t_capacity = t_capacity;
  matrix->e = malloc((rows + t_capacity) * sizeof(int*));
  assert(matrix->e && "Failed to allocate matrix row storage!");
  for (int i=0; i<(matrix->rows + matrix->t_capacity); i++) {
    int *r = malloc((matrix->cols + matrix->t_capacity) * sizeof(int));
    assert(r && "Failed to allocate matrix row!");
    memset(r, 0, (matrix->cols + matrix->t_capacity) * sizeof(int)); // zeroize all entries
    matrix->e[i] = r;
  }
}

void fmm_matrix_destroy(fmm_matrix *matrix) {
  for (int i=0; i<(matrix->rows + matrix->t_capacity); i++) {
    free(matrix->e[i]);
  }
  free(matrix->e);
}

void fmm_matrix_increase_capacity(fmm_matrix *matrix, int t_capacity_delta) {
  /* increase row storage capacity (old row pointers copied) */
  int old_row_dim = matrix->rows + matrix->t_capacity;
  int new_row_dim = old_row_dim + t_capacity_delta;
  int old_col_dim = matrix->cols + matrix->t_capacity;
  int new_col_dim = old_col_dim + t_capacity_delta;
  int **ee = realloc(matrix->e, new_row_dim * sizeof(int*));
  assert(ee && "fmm_matrix_increase_capacity: could not realloc row container");
  matrix->e = ee;
  for (int i=0; i<old_row_dim; i++) { // realloc old rows
    int *r = realloc(matrix->e[i], new_col_dim * sizeof(int));
    assert(r && "fmm_matrix_increase_capacity: could not realloc row");
    memset(r + old_col_dim, 0, (new_col_dim - old_col_dim) * sizeof(int)); // new entries zero-initialized
    matrix->e[i] = r;
  }
  for (int i=old_row_dim; i<new_row_dim; i++) { // malloc new rows
    int *r = malloc(new_col_dim * sizeof(int));
    assert(r && "fmm_matrix_increase_capacity: could not malloc row");
    memset(r, 0, new_col_dim * sizeof(int)); // new entries zero-initialized
    matrix->e[i] = r;
  }
  matrix->t_capacity += t_capacity_delta;
}

void fmm_matrix_copy(fmm_matrix *dst, fmm_matrix *src) {
  /* dst assumed to be initialized */
  /* resize if necessary */
  if (dst->cols != src->cols || dst->rows != src->rows) {
    fmm_matrix_destroy(dst);
    fmm_matrix_init(dst, src->rows, src->cols, src->t_capacity);
  }
  assert(dst->cols == src->cols && dst->rows == src->rows && "dimension mismatch");
  if (dst->t_capacity < src->t_capacity) {
    fmm_matrix_increase_capacity(dst, src->t_capacity - dst->t_capacity);
  }
  assert(dst->t_capacity >= src->t_capacity && "capacity mismatch");
  fmm_matrix_set_all_entries(dst, 0);
  for (int i=0; i<(src->rows + src->t_capacity); i++) {
    for (int j=0; j<(src->cols + src->t_capacity); j++) {
      fmm_matrix_entry_set(dst, i, j, fmm_matrix_entry(src, i, j));
    }
  }
  dst->t = src->t;
}

inline int fmm_matrix_entry(fmm_matrix *matrix, int row, int col) { return matrix->e[row][col]; }

inline void fmm_matrix_entry_set(fmm_matrix *matrix, int row, int col, int v) { matrix->e[row][col] = v; }

void fmm_matrix_set_all_entries(fmm_matrix *matrix, int v) {
  for (int i=0; i<(matrix->rows + matrix->t_capacity); i++) {
    for (int j=0; j<(matrix->cols + matrix->t_capacity); j++) {
      fmm_matrix_entry_set(matrix, i, j, v);
    }
  }
}

static void multiprint(const char *s, int n) {
  for (int i=0; i<n; i++) {
    printf(s);
  }
}

static void fmm_matrix_print_including_partial_t_space(fmm_matrix *matrix, const char *cell_format_string, int num_t) {
  for (int i=0; i<(matrix->rows + num_t); i++) {
    if (i == matrix->rows) {
      multiprint("-", (matrix->cols + num_t + 1) * 3);
      printf("\n");
    }
    for (int j=0; j<(matrix->cols + num_t); j++) {
      if (j == matrix->cols) {
        printf("  |");
      }
      printf(cell_format_string, fmm_matrix_entry(matrix, i, j));
    }
    printf("\n");
  }
}

void fmm_matrix_print(fmm_matrix *matrix, const char *cell_format_string) {
  fmm_matrix_print_including_partial_t_space(matrix, cell_format_string, 0);
}

void fmm_matrix_print_t(fmm_matrix *matrix, const char *cell_format_string) {
  fmm_matrix_print_including_partial_t_space(matrix, cell_format_string, matrix->t);
}

void fmm_matrix_print_full_capacity(fmm_matrix *matrix, const char *cell_format_string) {
  fmm_matrix_print_including_partial_t_space(matrix, cell_format_string, matrix->t_capacity);
}

int fmm_matrix_column_weight(fmm_matrix *matrix, int col) {
  int num_nonzero_entries = 0;
  for (int i=0; i<(matrix->rows + matrix->t); i++) {
    if (fmm_matrix_entry(matrix, i, col)) {
      num_nonzero_entries++;
    }
  }
  return num_nonzero_entries;
}

int fmm_matrix_num_additions_at_col(fmm_matrix *matrix, int col) {
  int num_nonzero_entries = fmm_matrix_column_weight(matrix, col);
  return num_nonzero_entries > 0 ? num_nonzero_entries - 1 : 0;
}

int fmm_matrix_num_additions(fmm_matrix *matrix) {
  int n = 0;
  for (int i=0; i<(matrix->cols + matrix->t); i++) {
    n += fmm_matrix_num_additions_at_col(matrix, i);
  }
  return n;
}

void fmm_matrix_column_matches_nonzero_equality(fmm_matrix *matrix, int col1, int col2, int *num_same_sign, int *num_opposite_sign) {
  int s = 0; // same sign counter
  int o = 0; // opposite sign counter
  for (int i=0; i<(matrix->rows + matrix->t); i++) {
    int e1 = fmm_matrix_entry(matrix, i, col1);
    int e2 = fmm_matrix_entry(matrix, i, col2);
    if (e1 && e2) {
      if (e1 == e2) {
        s++;
      }
      if (e1 == -e2) {
        o++;
      }
    }
  }
  *num_same_sign = s;
  *num_opposite_sign = o;
}
