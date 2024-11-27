#include "fmm_file_utils.h"
#include "fmm_matrix_storage_utils.h"
#include "fmm_alg_correctness_check.h"
#include "fmm_algorithm_specification.h"
#include "fmm_addition_reduction.h"
#include "fmm_brute_force.h"
#include <time.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#define LATEX_PRINT_DEFAULT 0
#define VERBOSE_FILE_READ   0
#define VERBOSE_CORRECTNESS 0
#define VERBOSE_REDUCTION   2

//#define REDUCTION_METHOD reduction_method_brute_force
//#define REDUCTION_METHOD reduction_method_greedy_vanilla
#define REDUCTION_METHOD reduction_method_greedy_potential

#if 1 /* same alpha paramater domain for A, B and C */
#define ALPHA_START    0.0
#define ALPHA_END      0.5
#define ALPHA_NUM_STEPS 5
#else
/* Individual alpha parameter domains for A, B and C */
/* These sample parameters will run separate alpha values for A, B and C, (a single alpha value per A, B or C matrix)
 * and these particular alphas are near-optimal for the Moosbauer beast 666 algoritm.
 * Near-optimal should be interpreted as the best that the potential heuristic can do given the current code base.
 * The order of moves can matter, since ties are resolved by taking the first move with the best heuristic score.
 * This is the reason why the Python implementation and the C implementations sometimes output slightly different,
 * results (an op or two at most according to our simulations).
 * For comparison, and for your convenience, here is the final output that I got when running greedy potential with
 * the given alphas on Moosbauers beast on my somewhat aged but still decent ASUS ZenBook laptop:
 *
 *   Algorithm uses 850 + 844 + 1254 = 2948 additions [naive]
 *   Algorithm uses 251 + 222 + 399 = 872 additions after reduction with alpha = ( 0.130000, 0.101000, 0.033400) [greedy potential]
 *                  ---------------------
 *   Total savings: 599 + 622 + 855 = 2076 additions (70.42%)
 *
 *   Reduction runtime: 46.16 sec
 */
/* A */
#define ALPHA_START_A     0.13
#define ALPHA_END_A       0.13
#define ALPHA_NUM_STEPS_A 1
/* B */
#define ALPHA_START_B     0.101
#define ALPHA_END_B       0.101
#define ALPHA_NUM_STEPS_B 1
/* C */
#define ALPHA_START_C     0.0334
#define ALPHA_END_C       0.0334
#define ALPHA_NUM_STEPS_C 1
#endif

/*
 * X -> Y/Z, (a, b, c)
 *
 * X = number of naive additions
 * Y = number of additions after greedy vanilla reduction
 * Z = number of additions after greedy potential reduction
 * a = best alpha_A found for greedy potential reduction
 * b = best alpha_B found for greedy potential reduction
 * c = best alpha_C found for greedy potential reduction
 */
#define ALG_FOLDER_OTHER "algorithms/other/"
//#define DEFAULT_FILE ALG_FOLDER_OTHER "does_not_exist.txt";
/* Z */
//#define DEFAULT_FILE ALG_FOLDER_OTHER "Strassen-222-7-18.txt"; // 18 -> 18/18
//#define DEFAULT_FILE ALG_FOLDER_OTHER "Strassen-222-7-24.txt"; // 24 -> 15/15
//#define DEFAULT_FILE ALG_FOLDER_OTHER "Smirnov-333-23-139.txt"; // 84 -> 68/68
#define DEFAULT_FILE ALG_FOLDER_OTHER "Laderman-333-23-98.txt"; // 98 -> 70/62, (0.1, 0.1, 0.2)
/* Z_2 */
//#define DEFAULT_FILE ALG_FOLDER_OTHER "Adaptive-455-73-mod2.m"; // 1,766 -> 516/480, ( 0.01, 0.072, 0.036)
//#define DEFAULT_FILE ALG_FOLDER_OTHER "Adaptive-555-94-mod2.m"; // 1,563 -> 507/492, (0.11, 0.132, 0.0255)

#define ALG_FOLDER_GREY "algorithms/grey/"
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-234-20-144.txt"; // 96 -> 63/62, (0.1, 0, 0)
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-322-11-50.txt"; // 22 -> 21/21
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-323-15-103.txt"; // 64 -> 44/44
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-422-14-84.txt"; // 48 -> 37/37
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-423-20-144.txt"; // 92 -> 59/58, (0, 0.1, 0)
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-522-18-99.txt"; // 53 -> 40/40

//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-333-23-142.txt"; // 87 -> 65/65
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-333-23-143.txt"; // 88 -> 66/64, (0, 0.1, 0)
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-333-23-144.txt"; // 89 -> 68/68
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-333-23-152.txt"; // 97 -> 64/63, (0.1, 0, 0)
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-333-23-221.txt"; // 166 -> 85/82, (0, 0, 0.1)
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-433-29-234.txt"; // 164 -> 99/98, (0, 0, 0.3)
//#define DEFAULT_FILE ALG_FOLDER_GREY "Grey-343-29-234.txt"; // 167 -> 103/101, (0, 0.4, 0)

#define ALG_FOLDER_MOOSBAUER "algorithms/moosbauer/"
/* Z */
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-223-11-mod0.m"; // 42 -> 31/31, (0,0,0)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-224-14-mod0.m"; // 99 -> 66/65, (0,0,0.16)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-225-18-mod0.m"; // 130 -> 79/77, (0,0.2,0.2)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-233-15-mod0.m"; // 76 -> 48/48, (0,0,0)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-234-20-mod0.m"; // 161 -> 100/99, (0.34,0,0)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-235-25-mod0.m"; // 199 -> 106/103, (0.3,0.3,0)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-244-26-mod0.m"; // 319 -> 175/173, (0,0.22,0)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-245-33-mod0.exp"; // 359 -> 177/172, (0.2, 0, 0.1)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-255-40-mod0.m"; // 682 -> 323/313, (0.10,0.12,0.06)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-266-56-mod0.m"; // 1133 -> 442/419, (0.2,0.1,0.08)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-334-29-mod0.m"; // 202 -> 107/105, (0.01,0.43,0)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-335-36-mod0.m"; // 370 -> 181/176, (0.01,0.24,0.1)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-336-42-mod0.m"; // 728 -> 274/264, (0,0.136,0.065)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-344-38-mod0.m"; // 401 -> 205/198, (0.2,0.145,0.17)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-345-47-mod0.m"; // 714 -> 284/276, (0,0.105,0.075)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-346-56-mod0.m"; // 795 -> 325/319, (0.235,0,0.085)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-355-58-mod0.m"; // 689 -> 331/330, (0,0,0.1)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-356-71-mod0.m"; // 950 -> 399/386, (0.255,0.11,0.14)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-444-49-mod0.m"; // 629 -> 335/325, (0.215,0.24,0.115)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-445-62-mod0.m"; // 1,175 -> 501/476, (0.18,0.133,0.08)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-446-74-mod0.m"; // 1,187 -> 488/467, (0,0.104,0.06)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-455-76-mod0.m"; // 1,107 -> 469/451, (0.17,0.1,0.105)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-456-93-mod0.m"; // 1,152 -> 567/544, (0.2025,0.13145,0.04)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-555-97-mod0.m"; // 1,710 -> 729/701, (0.002,0.162,0.061)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-555-97-mod0a.m"; // 2,208 -> 886/852, (0.002,0.087,0.056)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-556-116-mod0.m"; // 1,990 -> 747/712, (0.002,0.119,0.07)

/* Z_2 */
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-226-21-mod2.exp"; // 162 -> 74/71, (0, 0, 0.101)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-236-30-mod2.exp"; // 401 -> 159/149, ( 0.37, 0.125, 0.08)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-246-39-mod2.exp"; // 709 -> 253/242, (0, 0.01, 0.07)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-256-48-mod2.exp"; // 588 -> 240/228, (0.21, 0.101, 0.08)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-266-56-mod2.exp"; // 948 -> 332/227, (0, 0, 0.05)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-336-42-mod2.exp"; // 650 -> 231/222, (0.26, 0.09, 0.06)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-346-56-mod2.exp"; // 747 -> 293/278, (0.18, 0.101, 0.059)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-356-71-mod2.exp"; // 819 -> 335/322, (0.1,0.11,0.055)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-366-86-mod2.exp"; // 1,020 -> 415/385, (0.22, 0.06, 0.03)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-444-47-mod2.exp"; // 567 -> 234/223, (0.28, 0.13, 0.01)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-445-60-mod2.exp"; // 817 -> 293/282, (0.067, 0.132, 0.046)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-446-74-mod2.exp"; // 1,077 -> 377/360, (0.077, 0.063, 0.052)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-455-76-mod2.exp"; // 1,015 -> 388/372, (0.15, 0.11, 0.01)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-456-93-mod2.exp"; // 1,409 -> 508/465, (0.23, 0.123, 0.042)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-466-116-mod2.exp"; // 1,745 -> 621/581, (0.0001, 0.0931, 0.031)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-555-95-mod2.exp"; // 808 -> 418/403, (0.167, 0.209, 0.001)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-556-116-mod2.exp"; // 1,777 -> 601/580, (0.1269, 0.001, 0.059)
//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-566-144-mod2.exp"; // 2,325 -> 811/766, (0.05, 0.095, 0.035)

//#define DEFAULT_FILE ALG_FOLDER_MOOSBAUER "Moosbauer-666-164-mod2.exp"; // 2,948 -> 927/872 = 251 + 222 + 399, (0.13, 0.101, 0.0334)


static void get_commandline_args(int argc, char **argv, char *fn, reduction_method *red, int *verbosity_level_reduction, double *alpha_start, double *alpha_end, int *alpha_steps, int *print_latex) {
  if (argc > 2) { // reduction method
    const char *reduction_method_text = argv[2];
    assert(reduction_method_text && "unexpected reduction method parameter");
    if (!strcmp(reduction_method_text, "brute-force") || !strcmp(reduction_method_text, "bruteforce") || !strcmp(reduction_method_text, "bf")) {
      *red = reduction_method_brute_force;
      if (argc > 3) { // verbosity level
        *verbosity_level_reduction = atoi(argv[3]);
      }
      if (argc > 4) { // latex printing
        *print_latex = !strcmp("latex", argv[4]);
      }
    }
    if (!strcmp(reduction_method_text, "vanilla") || !strcmp(reduction_method_text, "greedy vanilla") || !strcmp(reduction_method_text, "gv") || !strcmp(reduction_method_text, "v")) {
      *red = reduction_method_greedy_vanilla;
      if (argc > 3) { // verbosity level
        *verbosity_level_reduction = atoi(argv[3]);
      }
      if (argc > 4) { // latex printing
        *print_latex = !strcmp("latex", argv[4]);
      }
    }
    if (!strcmp(reduction_method_text, "potential") || !strcmp(reduction_method_text, "greedy potential") || !strcmp(reduction_method_text, "gp") || !strcmp(reduction_method_text, "p")) {
      *red = reduction_method_greedy_potential;
      switch (argc) {
      case 4:
      case 5:
        *verbosity_level_reduction = atoi(argv[3]);
        if (argc > 4) { // latex printing
          *print_latex = !strcmp("latex", argv[4]);
        }
        break;
      case 6:
      case 7:
      case 8:
        {
          double a = atof(argv[3]);
          double b = atof(argv[4]);
          int c = atoi(argv[5]);
          for (int i=0; i<3; i++) {
            alpha_start[i] = a;
            alpha_end[i] = b;
            alpha_steps[i] = c;
          }
        }
        if (argc > 6) { // verbosity level
          *verbosity_level_reduction = atoi(argv[6]);
        }
        if (argc > 7) { // latex printing
          *print_latex = !strcmp("latex", argv[7]);
        }
        break;
      case 12:
      case 13:
      case 14:
        for (int i=0; i<3; i++) {
          alpha_start[i] = atof(argv[3+3*i]);
          alpha_end[i] = atof(argv[4+3*i]);
          alpha_steps[i] = atoi(argv[5+3*i]);
        }
        if (argc > 12) { // verbosity level
          *verbosity_level_reduction = atoi(argv[12]);
        }
        if (argc > 13) { // latex printing
          *print_latex = !strcmp("latex", argv[13]);
        }
        break;
      default:
        /* intentionally left blank */
        /* use default parameters in all other cases */
        break;
      }
    }
  }
}


int main(int argc, char **argv) {
  /* set default parameters */
  const char *file_name = DEFAULT_FILE;
  char fn[256];
  strncpy(fn, file_name, 256);
  if (argc > 1) { // file name
    strcpy(fn, argv[1]);
  }
  reduction_method red = REDUCTION_METHOD;
  int verbosity_level_reduction = VERBOSE_REDUCTION;
#if defined(ALPHA_START)
  double alpha_start[3] = {ALPHA_START, ALPHA_START, ALPHA_START};
  double alpha_end[3] = {ALPHA_END, ALPHA_END, ALPHA_END};
  int alpha_steps[3] = {ALPHA_NUM_STEPS, ALPHA_NUM_STEPS, ALPHA_NUM_STEPS};
#else
  double alpha_start[3] = {ALPHA_START_A, ALPHA_START_B, ALPHA_START_C};
  double alpha_end[3] = {ALPHA_END_A, ALPHA_END_B, ALPHA_END_C};
  int alpha_steps[3] = {ALPHA_NUM_STEPS_A, ALPHA_NUM_STEPS_B, ALPHA_NUM_STEPS_C};
#endif
  int print_latex = LATEX_PRINT_DEFAULT;

  /* replace default parameters with command-line arguments, if any */
  get_commandline_args(argc, argv, fn, &red, &verbosity_level_reduction, alpha_start, alpha_end, alpha_steps, &print_latex);

  /* read algorithm from file */
  fmm_alg alg;
  int t_capacity = 10;
  printf("Using file %s\n\n", fn);
  int read_error;
  if ((read_error = read_from_file(&alg, fn, t_capacity, VERBOSE_FILE_READ))) {
    printf("Could not read file %s, error %d\n", fn, read_error);
    return 1; // could not read algorithm specification data
  }

  /* print naive algorithm in compact matrix form */
  printf("Naive algorithm in compact form:\n");
  fmm_alg_print(&alg);
  printf("\n");

  /* verify correctness before reduction */
  int is_correct = fmm_alg_is_correct(&alg, VERBOSE_CORRECTNESS);
  printf("Algorithm correctness verification before reduction: %s\n\n", is_correct ? "ok" : "NOT OK");
  if (!is_correct) {
    return 1;
  }

  /* print latex code for naive algorithm (before reduction) */
  if (print_latex) {
    printf("Complete naive algorithm in LaTeX format:\n");
    fmm_alg_print_latex(&alg, "A", "B", "M", "C", "t", "u", "v");
    printf("\n");
  }

  /* count naive additions before reduction */
  const char *ABC[] = {"A", "B", "C"};
  fmm_matrix *mm[] = {&alg.A, &alg.B, &alg.C};
  int a[3], ra[3]; // number of additions before and after reduction
  for (int i=0; i<3; i++) {
    a[i] = fmm_matrix_num_additions(mm[i]);
  }

  /* reduce algorithm */
  clock_t start = clock();
  int k1[3], k2[3];
  switch (red) {
  case reduction_method_brute_force:
    for (int i=0; i<3; i++) {
      fmm_brute_force_matrix(mm[i], verbosity_level_reduction);
      ra[i] = fmm_matrix_num_additions(mm[i]);
    }
    break;
  case reduction_method_greedy_vanilla:
    for (int i=0; i<3; i++) {
      fmm_addition_reduction(mm[i], red, NULL, verbosity_level_reduction);
      ra[i] = fmm_matrix_num_additions(mm[i]);
    }
    break;
  case reduction_method_greedy_potential:
    for (int i=0; i<3; i++) {
      if (verbosity_level_reduction) { printf("%s:\n", ABC[i]); }
      find_best_greedy_potential_parameters(mm[i], alpha_start[i], alpha_end[i], alpha_steps[i], &k1[i], &k2[i], &ra[i], alpha_steps[i] > 2000 ? 1 : verbosity_level_reduction);
      printf("Best greedy potential params found for %s are alpha = %f, or ( k1, k2) = (%4d,%4d), reducing from %3d to %3d additions\n%s", ABC[i], (double)k2[i]/k1[i], k1[i], k2[i], a[i], ra[i], verbosity_level_reduction ? "\n" : "");
    }
    if (!verbosity_level_reduction) { printf("\n"); }
  }
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;

  /* print reduced algorithm in compact matrix form */
  printf("Reduced algorithm in compact form:\n");
  fmm_alg_print_t(&alg);
  printf("\n");

  /* verify correctness after reduction */
  is_correct = fmm_alg_is_correct(&alg, VERBOSE_CORRECTNESS);
  printf("Algorithm correctness verification after reduction: %s\n\n", is_correct ? "ok" : "NOT OK");
  if (!is_correct) {
    return 1;
  }

  /* print latex code for reduced algorithm */
  if (print_latex) {
    printf("Complete reduced algorithm in LaTeX format:\n");
    fmm_alg_print_latex(&alg, "A", "B", "M", "C", "t", "u", "v");
    printf("\n\n");
  }

  /* print stats */
  int t = a[0] + a[1] + a[2];
  int tt = ra[0] + ra[1] + ra[2];
  printf("Algorithm uses %3d + %3d + %3d = %3d additions [naive]\n", a[0], a[1], a[2], t);
  printf("Algorithm uses %3d + %3d + %3d = %3d additions after reduction", ra[0], ra[1], ra[2], tt);
  switch (red) {
  case reduction_method_brute_force:
    printf(" [brute force]");
    break;
  case reduction_method_greedy_vanilla:
    printf(" [greedy vanilla]");
    break;
  case reduction_method_greedy_potential:
    printf(" with alpha = ( %f, %f, %f) [greedy potential]", (double)k2[0]/k1[0], (double)k2[1]/k1[1], (double)k2[2]/k1[2]);
    break;
  }
  printf("\n");
  printf("               ---------------------\n");
  printf("Total savings: %3d + %3d + %3d = %3d additions (%5.2f%%)\n\n", a[0] - ra[0], a[1] - ra[1], a[2] - ra[2], t - tt, 100*((double)t - tt)/t);

  printf(seconds < 1 ? "Reduction runtime: %f sec\n\n" : "Reduction runtime: %.2f sec\n\n", seconds);

  fflush(stdout);
  fmm_alg_destroy(&alg);
  return 0;
}
