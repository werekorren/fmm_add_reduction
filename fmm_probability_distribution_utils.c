#include "fmm_probability_distribution_utils.h"
#include <malloc.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>//system
#include <string.h>

int probability_distribution_init(probability_distribution_int *pd, int start_value, int num_items) {
  pd->start_value = start_value;
  pd->num_items = num_items;
  pd->v = calloc(num_items, sizeof(int));
  return 0;
}

void probability_distribution_free(probability_distribution_int *pd) {
  if (pd) {
    free(pd->v);
  }
}

int probability_distribution_add(probability_distribution_int *pd, int value, int count) {
  assert(pd && "unexpected parameter, usage error");
  if (value < pd->start_value || (pd->start_value + pd->num_items) <= value) {
    return 1;
  }
  pd->v[value - pd->start_value] += count;
  return 0;
}

int probability_distribution_remove(probability_distribution_int *pd, int value, int count) {
  return probability_distribution_add(pd, value, -count);
}

static void get_probability_distribution_nonzero_min_max(probability_distribution_int *pd, int *min, int *max) {
  assert(pd && "unexpected parameter, usage error");
  assert(min && "unexpected parameter, usage error");
  assert(max && "unexpected parameter, usage error");
  *min = -1;
  *max = -1;
  /* find first non-zero entry */
  for (int i=0; i<pd->num_items; i++) {
    if (pd->v[i]) {
      *min = i;
      break;
    }
  }
  if (*min == -1) {
    return; /* zero distribution, nothing to print */
  }
  /* find last non-zero entry */
  for (int i=pd->num_items-1; i>=0; i--) {
    if (pd->v[i]) {
      *max = i;
      break;
    }
  }
}

static void get_probability_distribution_sum(probability_distribution_int *pd, int *sum) {
  assert(pd && "unexpected parameter, usage error");
  assert(sum && "unexpected parameter, usage error");
  int s = 0;
  for (int i=0; i<pd->num_items; i++) {
    s += pd->v[i];
  }
  *sum = s;
}

void probability_distribution_print(probability_distribution_int *pd, const char *delimiter) {
  assert(pd && "unexpected parameter, usage error");
  int min_non_zero, max_non_zero;
  get_probability_distribution_nonzero_min_max(pd, &min_non_zero, &max_non_zero);
  if (min_non_zero == -1) {
    return; /* zero distribution, nothing to print */
  }
  assert(min_non_zero != -1 && max_non_zero != -1 && "unexpected computation for min and max indices in probability distribution");
  for (int i=min_non_zero; i<max_non_zero; i++) {
    printf("%4d: %d%s", i, pd->v[i], delimiter);
  }
  printf("%4d: %d", max_non_zero, pd->v[max_non_zero]);
}

static char *unique_file_name(char *file_name, char *prefix) {
  time_t t = time(NULL);
  struct tm *currentTime = localtime(&t);
  sprintf(file_name, "%s_%04d%02d%02d_%02d%02d%02d.py",
    prefix,
    currentTime->tm_year + 1900, currentTime->tm_mon + 1, currentTime->tm_mday,
    currentTime->tm_hour, currentTime->tm_min, currentTime->tm_sec);
  return file_name;
}

int probability_distribution_python_plot(char *file_name, probability_distribution_int *pd, char *file_name_prefix, const char *python_commmand) {
  assert(pd && "unexpected parameter, usage error");
  assert(file_name_prefix && "unexpected parameter, usage error");
  int min_non_zero, max_non_zero;
  get_probability_distribution_nonzero_min_max(pd, &min_non_zero, &max_non_zero);
  if (min_non_zero == -1) {
    return 1; /* zero distribution, nothing to print */
  }
  assert(min_non_zero != -1 && max_non_zero != -1 && "unexpected computation for min and max indices in probability distribution");

  /* create output file */
  unique_file_name(file_name, file_name_prefix);
  int file_name_len = strlen(file_name);
  char file_name_png[256];
  strcpy(file_name_png, file_name);// including .py extension
  strcpy(file_name_png + file_name_len - 3, ".png");// replace .py extension with .png
  FILE *f = fopen(file_name, "w");
  if (!f) {
    return 2;
  }

  /* write python plot */
  int pd_sum;
  get_probability_distribution_sum(pd, &pd_sum);
  fprintf(f,
    "import matplotlib.pyplot as plt\n"
    "import seaborn as sns # pip include seaborn\n"
    "import numpy as np\n"
    "\n"
    "raw_data = { # %d entries\n", pd_sum);
  for (int i=min_non_zero; i<max_non_zero; i++) {
    fprintf(f, "%4d: %d,\n", i, pd->v[i]);
  }
  fprintf(f, "%4d: %d\n", max_non_zero, pd->v[max_non_zero]);
  fprintf(f,
    "}\n"
    "\n"
    "def main():\n"
    "    data = []\n"
    "    for num_additions, count in raw_data.items():\n"
    "        data.extend([num_additions] * count)\n"
    "\n"
    "    # create histogram\n"
    "    plt.figure(figsize=(9, 6)) # Optional: Adjust figure size\n"
    "    #ax = sns.histplot(data, bins=sorted(raw_data.keys()), kde=True, color='red', edgecolor='black', facecolor='skyblue', discrete=True)\n"
    "    ax = sns.histplot(data, bins=sorted(raw_data.keys()), kde=False, color='red', edgecolor='black', facecolor='skyblue', discrete=True)\n"
    "\n"
    "    # add labels and title\n"
    "    ax.set_xlabel(\"Num Additions After Reduction\")\n"
    "    ax.set_ylabel(\"Num FMM Algorithms\")\n"
    "    #ax.set_title(\"Addition Count Distribution\")\n"
    "    ax.set_xticks(sorted(raw_data.keys())) #  x-axis ticks to match addition count\n"
    "\n"
    "    # save plot to png file\n"
    "    plt.savefig('%s')\n"
    "\n"
    "    # display plot\n"
    "    #plt.show()\n"
    "\n"
    "if __name__ == \"__main__\":\n"
    "    main()\n", file_name_png);
  fclose(f);
  if (python_commmand) {
    char command[512];
    sprintf(command, "%s %s", python_commmand, file_name);
    system(command);
  }
  return 0;
}
