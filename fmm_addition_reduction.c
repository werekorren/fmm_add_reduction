#include "fmm_addition_reduction.h"
#include <assert.h>
#include <stdio.h>
#include <memory.h>

typedef void (*move_set_fitness_evaluation_function)(fmm_matrix *, int, int, reduction_parameters *, fmm_matrix *, int *, int *);

typedef struct {
  int step;
  reduction_method red;
  reduction_parameters par; /* optional parameters for some reduction methods */
  fmm_matrix move_overlap_value; // storage for move overlap evaluation (counts the number of times a variable substitution is applicable)
} reduction_state;

void fmm_reduction_move(fmm_matrix *A, int row1, int row2, int move_type) {
  /* make sure that there is enough storage room in the matrices */
  if (A->t == A->t_capacity) {
    fmm_matrix_increase_capacity(A, 10);
  }

  /* update new col in A */
  assert((move_type == 1 || move_type == -1) && "update_matrix_with_move: unexpected parameter move type");
  fmm_matrix_entry_set(A, row1, A->cols + A->t, 1);
  fmm_matrix_entry_set(A, row2, A->cols + A->t, move_type == 1 ? 1 : -1);

  /* update rows in A - scan for overlaps */
  for (int i=0; i<(A->cols + A->t); i++) { // for all columns
    int e1 = fmm_matrix_entry(A, row1, i);
    if (e1 == 0) {
      continue;
    }
    int e2 = fmm_matrix_entry(A, row2, i);
    if (e2 == 0) {
      continue;
    }
    if ((move_type == 1 && e1 == e2) ||
        (move_type == -1 && e1 == -e2)) {
      fmm_matrix_entry_set(A, row1, i, 0);
      fmm_matrix_entry_set(A, row2, i, 0);
      fmm_matrix_entry_set(A, A->rows + A->t, i, e1);
    }
  }

  /* we added a row and col to A, so update t-space dimension of A */
  A->t += 1;
}

void fmm_reduction_move_undo(fmm_matrix *A) {
  int row1 = -1, row2 = -1, move_type = 0;

  /* remove row and col of A (decrease t-space dimension) */
  A->t -= 1;
  int rowt = A->rows + A->t;
  int colt = A->cols + A->t;
  for (int i=0; i<rowt; i++) { // for all columns except the one that was added last
    int e = fmm_matrix_entry(A, i, colt);
    if (e == 0) { continue; }
    if (e == 1) { row1 = i; break; }
  }
  for (int i=rowt-1; i>row1; i--) { // for all columns except the one that was added last
    int e = fmm_matrix_entry(A, i, colt);
    if (e == 0) { continue; }
    row2 = i;
    move_type = e;
    break;
  }
  if (row1 >= row2) { printf("ERROR: row1 = %d, row2 = %d\n", row1, row2); }
  if (move_type != 1 && move_type != -1) { printf("ERROR: move_type = %d\n", move_type); }
  assert(row1 < row2 && "update_matrix_move_undo: unexpected row indices");
  assert(fmm_matrix_entry(A, row1, colt) == 1 && "update_matrix_move_undo: unexpected move entriy in row1 (not +1)");
  assert((move_type == 1 || move_type == -1) && "update_matrix_move_undo: unexpected move type");
  assert(((move_type == 1 && fmm_matrix_entry(A, row2, colt) == 1) ||
          (move_type == -1 && fmm_matrix_entry(A, row2, colt) == -1)) &&
            "update_matrix_move_undo: unexpected move entry in row2 (not +/-1)");

  /* nullify last col in A */
  fmm_matrix_entry_set(A, row1, colt, 0);
  fmm_matrix_entry_set(A, row2, colt, 0);

  /* revert rows in A - scan for overlaps */
  for (int i=0; i<colt; i++) { // for all columns
    int e1 = fmm_matrix_entry(A, rowt, i);
    if (e1 == 0) { continue; }
    fmm_matrix_entry_set(A, rowt, i, 0);
    fmm_matrix_entry_set(A, row1, i, e1);
    fmm_matrix_entry_set(A, row2, i, move_type == 1 ? e1 : -e1);
  }
}

/* returns the number of times that a substitution of type t = A_row1 +/- A_row2 can be applied */
/* (if it can be applied t times, it can save t-1 addition operations) */
void compute_move_value_greedy_vanilla(fmm_matrix *A, int row1, int row2, int *p, int *n) {
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

static void update_move_overlap_values_for_one_move_set(reduction_state *state, fmm_matrix *A, int i, int j) {
  int p, n;
  //get_fitness_values_for_one_move_set(state, A, i, j, &p, &n);

  assert(i < j && "get_fitness_values_for_one_move_set: unecpected parameters i and j");
  compute_move_value_greedy_vanilla(A, i, j, &p, &n);

  fmm_matrix_entry_set(&state->move_overlap_value, i, j, p); // fitness value for addition in top right half
  fmm_matrix_entry_set(&state->move_overlap_value, j, i, n); // fitness value for subtraction in bottom left half
}

static void update_all_move_overlap_values(reduction_state *state, fmm_matrix *A) {
  int t = state->move_overlap_value.rows + state->move_overlap_value.t;
  for (int i=0; i<t; i++) {
    for (int j=i+1; j<t; j++) { // for all pairs or rows
      update_move_overlap_values_for_one_move_set(state, A, i, j);
    }
  }
}

static void init_move_overlap_values(reduction_state *state, fmm_matrix *A) {
  fmm_matrix_init(&state->move_overlap_value, A->rows, A->rows, A->t_capacity); // allocate storage, square matrix since we consider pairs of A_i's (row index corresponds to an A_i)
  update_all_move_overlap_values(state, A);
}

static void fmm_addition_reduction_destroy_state(reduction_state *state) {
  fmm_matrix_destroy(&state->move_overlap_value);
}

static void print_move_overlap_values(reduction_state *state) {
  printf("move values =\n");
  fmm_matrix_print_t(&state->move_overlap_value, " %2d");
}

static void update_all_affected_move_overlap_values(reduction_state *state, fmm_matrix *A, int row1, int row2) {
  int t = state->move_overlap_value.cols + state->move_overlap_value.t;
  int tv[t-1];
  memset(tv, 0, sizeof(int) * (t-1));
  for (int i=0; i<t-1; i++) {
    for (int j=i+1; j<t-1; j++) { // for all pairs or rows
      if (i != row1 && j != row1 && i != row2 && j != row2) {
        continue; // no updates needed for non-overlapping
      }
      /* update (i, j, op) only if it has potential */
      int p = fmm_matrix_entry(&state->move_overlap_value, i, j);
      int n = fmm_matrix_entry(&state->move_overlap_value, j, i);
      if (p > 1 || n > 1) {
        update_move_overlap_values_for_one_move_set(state, A, i, j);
        /* if i OR j is in {row1, row2}, then (i, t-1, op) and (j, t-1, op) may both need to be updated */
        if (i != row1 && i != row2) { // only if i does not coincide with move indices
          tv[i] = 1;
        }
        if (j != row1 && j != row2) { // only if j does not coincide with move indices
          tv[j] = 1;
        }
      }
    }
  }
  for (int i=0; i<t-1; i++) {
    if (tv[i]) {
      update_move_overlap_values_for_one_move_set(state, A, i, t-1); // update (i, t-1, op)
    }
  }
}

static int best_move_value_greedy_vanilla(reduction_state *state, int *row1, int *row2, int *move_type) {
  /* overlap counters assumed to be updated */
  /* lower left half stores fitness values for subtraction moves */
  /* top right half stores fitness values for addition moves */
  /* all values are non-negative, higher fitness vale is better */
  int best = 0;
  int t = state->move_overlap_value.rows + state->move_overlap_value.t;
  for (int i=0; i<t; i++) {
    for (int j=i+1; j<t; j++) { // for all pairs or rows in A
      int move_value = fmm_matrix_entry(&state->move_overlap_value, i, j);
      if (move_value > best) {
        *row1 = i;
        *row2 = j;
        best = move_value;
        *move_type = 1; // move corresponds to addition of variables
      }
      move_value = fmm_matrix_entry(&state->move_overlap_value, j, i);
      if (move_value > best) {
        *row1 = i;
        *row2 = j;
        best = move_value;
        *move_type = -1; // move corresponds to subtraction of variables
      }
    }
  }
  return best; // return best move value
}

/* potential = total num saved ops across all remaining operations */
static int total_potential_after_new_move(fmm_matrix *A, reduction_parameters *par, int row1, int row2, fmm_matrix *old_move_value) {
  int s = 0;
  int p, n;
  int t = A->rows + A->t;
  int tv[t-1];
  memset(tv, 0, sizeof(int) * (t-1));
  for (int i=0; i<t - 1; i++) {
    for (int j=i+1; j<t - 1; j++) {

      p = fmm_matrix_entry(old_move_value, i, j);
      n = fmm_matrix_entry(old_move_value, j, i);
      if (i != row1 && j != row1 && i != row2 && j != row2) {
        s += p < 2 ? 0 : p - 1; // num saved ops for positive move
        s += n < 2 ? 0 : n - 1; // num saved ops for negative move
        continue;
      }

      if (p > 1 || n > 1) {
        compute_move_value_greedy_vanilla(A, i, j, &p, &n);
        if (i != row1 && i != row2) {
          tv[i] = 1;
        }
        if (j != row1 && j != row2) {
          tv[j] = 1;
        }
      }
      s += p < 2 ? 0 : p - 1; // num saved ops for positive move
      s += n < 2 ? 0 : n - 1; // num saved ops for negative move

      /* if i & j are unaffected by the move, use previously computed fitness values */
      /* if i & j are affected by the move, recompute fitness value */
    }
  }
  for (int i=0; i<t-1; i++) {
    if (tv[i]) {
      compute_move_value_greedy_vanilla(A, i, t-1, &p, &n);
      s += p < 2 ? 0 : p - 1; // num saved ops for positive move
      s += n < 2 ? 0 : n - 1; // num saved ops for negative move
    }
  }
  return s;
}

static void compute_move_value_greedy_potential(fmm_matrix *A, int row1, int row2, reduction_parameters *par, fmm_matrix *old_move_value, int *p, int *n) {
  /* compute positive and negative overlaps */
  int pp, nn, potential = 0;
  //fitness_function_greedy_vanilla(A, row1, row2, par, old_move_value, &pp, &nn);
  assert(row1 < row2 && "unexpected paramaters row1 and row2");
  pp = fmm_matrix_entry(old_move_value, row1, row2);
  nn = fmm_matrix_entry(old_move_value, row2, row1);
  *p = pp; // overlap (not counting potential yet)
  *n = nn; // overlap (not counting potential yet)

  /* abort if no additions can be reduced */
  if (pp < 2 && nn < 2) {
    return; // return best overlap, but without zero potential
  }

  /* compute potential for positive move; t = A_row1 + A_row2 */
  int k1 = par->_potential.k1;
  int k2 = par->_potential.k2;
  if (pp >= 2) {
    if (k2) {
      fmm_reduction_move(A, row1, row2, 1); // update matrix A with move (row1, row2, +)
      potential = total_potential_after_new_move(A, par, row1, row2, old_move_value); // compute potential after move
      fmm_reduction_move_undo(A); // restore matrix A
    }
    *p = k1*pp + k2*potential;
  }

  /* compute potential for negative move; t = A_row1 - A_row2 */
  if (nn >= 2) {
    if (k2) {
      fmm_reduction_move(A, row1, row2, -1); // update matrix A with move (row1, row2, -)
      potential = total_potential_after_new_move(A, par, row1, row2, old_move_value); // compute potential after move
      fmm_reduction_move_undo(A); // restore matrix A
    }
    *n = k1*nn + k2*potential;
  }
}

static int best_move_value_greedy_potential(fmm_matrix *A, reduction_state *state, int *row1, int *row2, int *move_type) {
  /* overlap counters assumed to be updated */
  /* lower left half stores fitness values for subtraction moves */
  /* top right half stores fitness values for addition moves */
  /* all values are non-negative, higher fitness vale is better */

  int best = 0;
  int t = state->move_overlap_value.rows + state->move_overlap_value.t;
  for (int i=0; i<t; i++) {
    for (int j=i+1; j<t; j++) { // for all pairs or rows in A
      int move_value_p, move_value_n;
      compute_move_value_greedy_potential(A, i, j, &state->par, &state->move_overlap_value, &move_value_p, &move_value_n);
      if (move_value_p > best) {
        *row1 = i;
        *row2 = j;
        best = move_value_p;
        *move_type = 1; // move corresponds to addition of variables
      }
      if (move_value_n > best) {
        *row1 = i;
        *row2 = j;
        best = move_value_n;
        *move_type = -1; // move corresponds to subtraction of variables
      }
    }
  }
  return best; // return best move value
}


#define DEFAULT_K1 5
#define DEFAULT_K2 1
static void fmm_addition_reduction_init_state(reduction_state *state, fmm_matrix *A, reduction_method red, reduction_parameters *red_par) {
  state->step = 0;
  state->red = red;
  switch (red) {
  case reduction_method_brute_force:
    /* nothing to do here */
    break;
  case reduction_method_greedy_vanilla:
    /* nothing to do here */
    break;
  case reduction_method_greedy_potential:
    state->par._potential.k1 = red_par ? red_par->_potential.k1 : DEFAULT_K1; /* set optional parameters (if any) for greedy potential */
    state->par._potential.k2 = red_par ? red_par->_potential.k2 : DEFAULT_K2; /* set optional parameters (if any) for greedy potential */
    break;
  default: /* intentional fall-through */
    assert(0 && "reduction method not implemented");
    break;
  }
  init_move_overlap_values(state, A); // note: move_set_fitness_eval needs to be set
}

/* returns 1 if reduction was applied, or zero if no additions could be reduced */
static int greedy_step(reduction_state *state, fmm_matrix *A) {
  /* find best pair of variables to combine (assumes that overlap counters with fitness values are updated) */
  int row1 = 0, row2 = 0, move_type = 0; // assignments to avoid compiler warnings only
  int move_value = 0; // assignments to avoid compiler warnings only

  switch (state->red) {
  case reduction_method_greedy_vanilla:
    move_value = best_move_value_greedy_vanilla(state, &row1, &row2, &move_type); // number of pattern overlaps is used as move value for greedy vanilla
    break;
  case reduction_method_greedy_potential:
    move_value = best_move_value_greedy_potential(A, state, &row1, &row2, &move_type); // weighted sum of number of pattern overlaps and remaining potential is used as move value for greedy potential
    break;
  default:
    assert(0 && "unsupported reduction method");
  }
  if (move_value < 2) {
    return 0; // no improvement could be found
  }

  /* update affected rows and columns in A */
  /* stores corresponding variable substitution as col t+1 (appends new row/col) */
  fmm_reduction_move(A, row1, row2, move_type);

  /* update move values correspondingly */
  if (state->move_overlap_value.t == state->move_overlap_value.t_capacity) { // make sure that there is enough storage room in the matrices
    fmm_matrix_increase_capacity(&state->move_overlap_value, 10);
  }
  state->move_overlap_value.t += 1;

  /* update all move overlap evaluations that need to be updated */
  update_all_affected_move_overlap_values(state, A, row1, row2);

  /* update step counter */
  state->step += 1;
  return 1; // indicate that non-trivial reduction was performed
}

void fmm_addition_reduction(fmm_matrix *A, reduction_method red, reduction_parameters *red_par, int verbose) {
  reduction_state state;
  int a = fmm_matrix_num_additions(A);
  fmm_addition_reduction_init_state(&state, A, red, red_par);
  if (verbose) {
    printf("Applying addition reduction\nHeuristic: ");
    switch (red) {
    case reduction_method_greedy_vanilla: printf("greedy vanilla"); break;
    case reduction_method_greedy_potential: printf("greedy potential"); break;
    default: printf("unknown");
    }
    printf("\n\nBefore reduction\nInitial use of %d additions\n", a);
    if (verbose >= 2){
      switch (verbose) {
      case 2: fmm_matrix_print_t(A, " %2d"); break;
      case 3: fmm_matrix_print_full_capacity(A, " %2d"); break;
      }
      print_move_overlap_values(&state);
      printf("********************************************************************************************************************************************************");
    }
    printf("\n");
  }
  for (;;) {
    int r = greedy_step(&state, A);
    if (r == 0) { // if no reduced additions
      if (verbose) {
        printf("Addition reduction iteration %d\n", state.step + 1);
        printf("No further improvement, exiting\n");
        int aa = fmm_matrix_num_additions(A);
        printf("Overall reduction from %d to %d additions (reduced %d)\n\n", a, aa, a - aa);
      }
      fmm_addition_reduction_destroy_state(&state);
      return;
    }
    if (verbose) {
      printf("Addition reduction iteration %d\n", state.step);
      printf("Reduced %d addition%s\n", r, r > 1 ? "s" : "");
      printf("Now using %d additions\n", fmm_matrix_num_additions(A));
      if (verbose >= 2){
        switch (verbose) {
        case 2: fmm_matrix_print_t(A, " %2d"); break;
        case 3: fmm_matrix_print_full_capacity(A, " %2d"); break;
        }
        print_move_overlap_values(&state);
        printf("********************************************************************************************************************************************************");
      }
      printf("\n");
    }
  }
}

static int myround(double n) { return (int)(n + 0.5); } // for positive numbers only
//static int myround(double n) { return (int)(n + 0.5 - (n<0)); } // for positive or negative numbers

#define MAX_DENOMINATOR 200000
static void find_approx_fraction(double f, int *nom, int *denom) {
  int best_nom = myround(f);
  int best_denom = 1;
  double best_fractional_part = f > best_nom ? f - best_nom : best_nom - f;
  /* try all denominators up to 1000 to see if they provide a better approximation */
  for (int i=2; i<(MAX_DENOMINATOR+1); i++) {
    double ff = f * i;
    int a = myround(ff);
    double fractional_part = a > ff ? a - ff : ff - a;
    if (fractional_part < best_fractional_part) {
      best_fractional_part = fractional_part;
      best_nom = a;
      best_denom = i;
    }
  }
  /* best_nom and best_denom can optionally be checked for common factors afterwards (here) if it is important to obtain a fully reduced fraction */
  /* performance not needed here, so keeping it really simple */
  for (int i=2; i<MAX_DENOMINATOR/2+1; i++) {
    while ((best_nom % i == 0) && (best_denom % i == 0)) {
      best_nom /= i;
      best_denom /= i;
    }
  }
  *nom = best_nom;
  *denom = best_denom;
}

void find_best_greedy_potential_parameters(fmm_matrix *A, double alpha_start, double alpha_end, int num_steps, int *k1, int *k2, int *num_add_after_reduction, int verbose) {
  fmm_matrix M; // perform reduction on local copy of A
  fmm_matrix_init(&M, A->rows, A->cols, A->t_capacity);
  fmm_matrix_copy(&M, A);
  fmm_matrix bestM; // best reduction of A stored and returned (avoids later recomputation)
  fmm_matrix_init(&bestM, A->rows, A->cols, A->t_capacity);
  fmm_matrix_copy(&bestM, A);
  int a = fmm_matrix_num_additions(&M);
  if (verbose == 1) {
    printf("Best number of additions found: %d", a);
  }
  int best_k1 = 1;
  int best_k2 = 0;
  int best_a = a;
  double alpha_step = (alpha_end - alpha_start) / num_steps;
  double f = alpha_start;
  do {
    fmm_matrix_copy(&M, A);
    reduction_parameters rp;
    find_approx_fraction(f, &rp._potential.k2, &rp._potential.k1);
    fmm_addition_reduction(&M, reduction_method_greedy_potential, &rp, 0);
    int aa = fmm_matrix_num_additions(&M);
    if (aa < best_a) {
      best_a = aa;
      best_k1 = rp._potential.k1;
      best_k2 = rp._potential.k2;
      fmm_matrix_copy(&bestM, &M); // save best reduction
      if (verbose == 1) {
        printf(" --> %d", aa);
      }
    }
    if (verbose > 1) {
      printf("For alpha = %f = k2/k1 = [approx] = %4d/%5d = %f: %4d -> %4d\n", f, rp._potential.k2, rp._potential.k1, (double)rp._potential.k2/rp._potential.k1, a, aa);
    }
    f += alpha_step;
  } while (f < alpha_end + alpha_step/2);
  if (verbose == 1) {
    printf("\n");
  }
//  printf("finished with reduction loop\n");
  fmm_matrix_destroy(&M);
  fmm_matrix_copy(A, &bestM); // copy best reduction into A
  fmm_matrix_destroy(&bestM);
  *k1 = best_k1;
  *k2 = best_k2;
  *num_add_after_reduction = best_a;
}
