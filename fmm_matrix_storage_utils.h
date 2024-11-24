#ifndef FMM_MATRIX_STORAGE_UTILS_H
#define FMM_MATRIX_STORAGE_UTILS_H

typedef struct {
  int rows;
  int cols;
  int t;
  int t_capacity;
  int **e; // entries
} fmm_matrix;

void fmm_matrix_init(fmm_matrix *matrix, int rows, int cols, int t_capacity);
void fmm_matrix_destroy(fmm_matrix *matrix);
void fmm_matrix_increase_capacity(fmm_matrix *matrix, int capacity_delta);
void fmm_matrix_copy(fmm_matrix *dst, fmm_matrix *src);

int fmm_matrix_entry(fmm_matrix *matrix, int row, int col);
void fmm_matrix_entry_set(fmm_matrix *matrix, int row, int col, int value);
void fmm_matrix_set_all_entries(fmm_matrix *matrix, int value);

void fmm_matrix_print(fmm_matrix *matrix, const char *cell_format_string);
void fmm_matrix_print_t(fmm_matrix *matrix, const char *cell_format_string);
void fmm_matrix_print_full_capacity(fmm_matrix *matrix, const char *cell_format_string);

int fmm_matrix_column_weight(fmm_matrix *matrix, int col);
int fmm_matrix_num_additions_at_col(fmm_matrix *matrix, int col);
int fmm_matrix_num_additions(fmm_matrix *matrix);

void fmm_matrix_column_matches_nonzero_equality(fmm_matrix *matrix, int col1, int col2, int *num_same_sign, int *num_opposite_sign);

#endif
