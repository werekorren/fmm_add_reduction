#include "fmm_file_utils.h"
#include "fmm_matrix_storage_utils.h"
#include "fmm_algorithm_specification.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* read parameters from file name
 * Assumptions (lots): using Grey txt file format and naming conventions, klm parameters are each single digit, dash used as delimiter in file name and no dashes in path.
 */
static int read_params_from_file_name(const char *file_name, int *zp, int *k, int *l, int *m, int *q, char * name) {
  int len = strlen(file_name);
  if (len >= 256) {
    return 1; // file name longer than expected (assuming max 256 characters in code below)
  }

  /* get zp */
  *zp = 0;
  char *p = strstr(file_name, "mod");
  if (p) {
    *zp = (int)(p[3] - '0'); // could not find any dash delimiter in file name
    assert(0 <= *zp && *zp <= 9 && "could not detect parameter p");
    if (*zp < 0 || 9 < *zp) {
      return 3; // could not detect parameter p
    }
  }

  /* get name */
  p = strchr(file_name, '-'); // find first dash (if any)
  if (!p) {
    return 2; // could not find any dash delimiter in file name
  }
  int first_dash_index = (int)(p - file_name);
  char *s = strrchr(file_name, '/'); // find last slash (if any)
  int last_slash_index = s ? (int)(s - file_name) : -1;
  int name_len = first_dash_index - last_slash_index - 1;
  strncpy(name, file_name + last_slash_index + 1, name_len);
  name[name_len] = 0;

  /* get k, l, m */
  *k = (int)(p[1] - '0');
  *l = (int)(p[2] - '0');
  *m = (int)(p[3] - '0');

  /* get q */
  char *pp = strchr(p + 1, '-');
  if (!pp) {
    return 3; // could not find second dash delimiter in file name
  }
  char *ppp = strchr(pp + 1, '-');
  if (!ppp) {
    return 4; // could not find third dash delimiter in file name
  }
  char t[256];
  int q_len = (int)(ppp - pp - 1);
  strncpy(t, pp+1, q_len);
  *q = atoi(t);
  return 0;
}

static int looks_like_a_number(char *buf) {
  return *buf == '-' || (*buf >= '0' && *buf <= '9');
}

static char *skip_to_next_integer(char *buf) {
  while (*buf != '\0' && looks_like_a_number(buf)) {
    buf++;
  }
  while (*buf != '\0' && !looks_like_a_number(buf)) {
    buf++;
  }
  return *buf == '\0' ? NULL : buf;
}

static int read_alg_from_grey_file(const char *file_name, fmm_alg *alg, int verbose) {
  FILE *f = fopen(file_name, "r");
  if (!f) {
    if (verbose) {
      printf("Could not open Grey txt file %s\n", file_name);
    }
    return 1;
  }

  int state1 = 0; // counts which row of actual matrix data we are at (global line counter for file)
  int state2 = 0; // counts which matrix we are at (0 for A, 1 for B, 2 for M)
  int state3 = 0; // counts which matrix row we are at (zero-indexed)
  int max_line_len = 2048;
  char line[max_line_len];
  while(fgets(line, max_line_len, f)) {
    int len = (int)strlen(line);
    if (len == 0 || line[0] == '#') { // skip empty lines or lines beginning with a # comment character
      if (verbose) {
        printf("[line length %3d,             comment] %s", len, line);
      }
      continue; // no state updates when skipping comment rows
    }
    if (state1 == alg->k * alg->l ||
        state1 == alg->k * alg->l + alg->l * alg->m) { // reached data row for matrix B or M
      state2++;
      state3 = 0;
    }
    if (verbose) {
      const char *matrix_name[] = {"A", "B", "M"};
      const char *rc[] = {"row", "row", "col"};
      printf("[line length %3d, %s %3d of matrix %s] %s", len, rc[state2], state3, matrix_name[state2], line);
    }

    /* read numbers from line here, put into A, B and M storage */
    char *buf = line;
    for (int i=0; i<alg->q; i++) { // every line is expected to contain q integers
      if (!looks_like_a_number(buf)) {
        return 1;
      }
      int n = atoi(buf);
      switch (state2) {
      case 0: fmm_matrix_entry_set(&(alg->A), state3, i, n); break;
      case 1: fmm_matrix_entry_set(&(alg->B), state3, i, n); break;
      case 2: fmm_matrix_entry_set(&(alg->C), i, state3, n); break; // transposing on purpose
      default: exit(1); /* unexpected error when reading file */
      }
      buf = skip_to_next_integer(buf);
    }

    /* update state counters */
    state1++;
    state3++;
  }

  fclose(f);
  return 0;
}

/* moosebauer file format transposes c entries, so this adjusts the element placement enumeration */
static int transposed_running_index(int rows, int cols, int i) {
  int r = i / rows;
  int c = i % rows;
  return cols * c + r;
}

static int read_alg_from_moosbauer_m_file(const char *file_name, fmm_alg *alg, int verbose) {
  FILE *f = fopen(file_name, "r");
  if (!f) {
    if (verbose) {
      printf("Could not open Moosbauer m file %s\n", file_name);
    }
    return 1; // could not open file
  }

  int matrix_row = 0; // counts which row of actual matrix data we are at (global line counter for file)
  int max_line_len = 2048;
  char line[max_line_len];
  char *buf = line;
  int which_matrix = 0; // counts which matrix we are at (0 for A, 1 for B, 2 for M)
  int num_read = 0; // counts how many integers we have read from the current line
  int num_read_a = 0, num_read_b = 0, num_read_m = 0; // counts how many integers we have read from the current line for each matrix
  int kl = alg->k * alg->l;
  int lm = alg->l * alg->m;
  int km = alg->k * alg->m;
  const int tot_to_read = (alg->k * alg->l + alg->l * alg->m + alg->k * alg->m) * alg->q;
  if (!fgets(line, max_line_len, f)) {
    if (verbose) {
      printf("Could not read from Moosbauer m file %s\n", file_name);
    }
    return 2; // could not read from file
  }
  if (verbose) {
    int len = (int)strlen(line);
    printf("[line length %3d, matrix data] %s", len, line);
  }

  /* read numbers from line here, put into A, B and M storage */
  for (;;) { // every line is expected to contain q integers
    buf = skip_to_next_integer(buf);
    if (!buf) { // reached end prematurely, no more numbers to read on this line, so we need to read the next line and continue
      if (!fgets(line, max_line_len, f)) {
        break; // end of file reached
      }
      buf = skip_to_next_integer(line);
    }
    int n = atoi(buf);
    switch (which_matrix) {
    case 0:
      fmm_matrix_entry_set(&(alg->A), num_read_a, matrix_row, n);
      num_read_a++;
      if (num_read_a % kl == 0) {
        num_read_a = 0;
        which_matrix++;
      }
      break;
    case 1:
      fmm_matrix_entry_set(&(alg->B), num_read_b, matrix_row, n);
      num_read_b++;
      if (num_read_b % lm == 0) {
        num_read_b = 0;
        which_matrix++;
      }
      break;
    case 2:
      fmm_matrix_entry_set(&(alg->C), matrix_row, transposed_running_index(alg->k, alg->m, num_read_m), n);
      num_read_m++;
      if (num_read_m % km == 0) {
        num_read_m = 0;
        which_matrix = 0;
        matrix_row++;
      }
      break; // transposing on purpose
    }
    if (++num_read == tot_to_read) {
      break; // all numbers read from this line, so we can proceed to next line
    }
  }

  fclose(f);
  return 0;
}

static int is_abc(char ch) { return ch == 'a' || ch == 'b' || ch == 'c'; }

static int abc_2_matrix_index(char ch) {
  switch (ch) {
  case 'a': return 0;
  case 'b': return 1;
  case 'c': return 2;
  }
  return -1;
}

static const char *next_abc(const char *buf) {
  for (;; buf++) {
    if (*buf == '\0') {
      return NULL;
    }
    if (is_abc(*buf)) {
      return buf; // found a, b or c
    }
  }
}

static const char * next_exp_symbol(const char *buf, int *coeff, int *which_matrix, int *row, int *col) {
  int len = strlen(buf);

  const char *abc = next_abc(buf);
  if (!abc) {
    return NULL;
  }

  int m = abc_2_matrix_index(*abc);
  assert(0 <= m && m <= 2 && "unexpected matrix classification");

  if (buf + len - abc < 2) {
    return NULL;
  }
  int r = (int)(abc[1] - '0');
  int c = (int)(abc[2] - '0');
  assert(0 <= r && r <= 9 && "unexpected row index");
  assert(0 <= r && r <= 9 && "unexpected col index");

  int k = 0;
  int neg = 0;
  for (const char *p=buf; p<abc; p++) {
    if (*p == '-') { // minus sign
      neg = 1;
    }
    if ('0' <= *p && *p <= '9') { // digit
      k = 10 * k + (int)(*p - '0');
    }
  }
  if (k == 0) { // no coefficient found, so we go with 1
    k = 1;
  }

  *coeff = neg ? -k : k;
  *which_matrix = m;
  *row = r;
  *col = c;

  return abc + 3;
}

static int read_alg_from_moosbauer_exp_file(const char *file_name, fmm_alg *alg, int verbose) {
  FILE *f = fopen(file_name, "r");
  if (!f) {
    if (verbose) {
      printf("Could not open Moosbauer exp file %s\n", file_name);
    }
    return 1; // could not open file
  }

  int matrix_row = 0; // counts which row of actual matrix data we are at (global line counter for file)
  int max_line_len = 2048;
  char line[max_line_len];

  /* read numbers from line here, put into A, B and M storage */
  while (fgets(line, max_line_len, f)) {
    int len = (int)strlen(line);
    if (len == 0 || line[0] == '#') { // skip empty lines or lines beginning with a # comment character
      if (verbose) {
        printf("[line length %3d,             comment] %s", len, line);
      }
      continue; // no state updates when skipping comment rows
    }

    const char *buf = line;
    int which_matrix = 0; // counts which matrix we are at (0 for A, 1 for B, 2 for M)
    /* extract, in turn, all a's, b's and c's on this line */
    int coeff, row, col;
    while ((buf = next_exp_symbol(buf, &coeff, &which_matrix, &row, &col))) {
      /* store in matrix */
      switch (which_matrix) {
      case 0:
        fmm_matrix_entry_set(&(alg->A), alg->l*(row-1) + (col-1), matrix_row, coeff);
        break;
      case 1:
        fmm_matrix_entry_set(&(alg->B), alg->m*(row-1) + (col-1), matrix_row, coeff);
        break;
      case 2:
        fmm_matrix_entry_set(&(alg->C), matrix_row, transposed_running_index(alg->k, alg->m, alg->k*(row-1) + (col-1)), coeff);
        break; // transposing on purpose
      }
    }
    matrix_row++;
  }
  assert(matrix_row == alg->q && "unexpected number of lines in input file");
  if (matrix_row != alg->q) {
    if (verbose) {
      printf("Unexpected number of lines in input file %s\n", file_name);
    }
    return 2; // unexpected number of lines in input file
  }

  fclose(f);
  return 0;
}

static int read_alg_from_file(const char *file_name, fmm_alg *alg, int verbose) {
  int len = strlen(file_name);
  if (len > 4 && !strcmp(file_name + len - 4, ".txt")) {
    return read_alg_from_grey_file(file_name, alg, verbose);
  }
  if (len > 2 && !strcmp(file_name + len - 2, ".m")) {
    return read_alg_from_moosbauer_m_file(file_name, alg, verbose);
  }
  if (len > 4 && !strcmp(file_name + len - 4, ".exp")) {
    return read_alg_from_moosbauer_exp_file(file_name, alg, verbose);
  }
  return -1; // unsupported file format
}

/* alg is assumed to not be initialized */
int read_from_file(fmm_alg *alg, const char *file_name, int t_capacity, int verbose) {
  int p, k, l, m, q;
  if (read_params_from_file_name(file_name, &p, &k, &l, &m, &q, alg->name)) {
    return 1; // failed to extract parameters from file name
  }
  fmm_alg_init(alg, p, k, l, m, q, t_capacity);
  /* populate matrices */
  return read_alg_from_file(file_name, alg, verbose);
}
