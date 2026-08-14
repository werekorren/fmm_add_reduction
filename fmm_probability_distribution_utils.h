#ifndef FMM_PROBABILITY_DISTRIBUTION_UTILS
#define FMM_PROBABILITY_DISTRIBUTION_UTILS

typedef struct {
  int start_value;
  int num_items;
  int *v;
} probability_distribution_int;

int probability_distribution_init(probability_distribution_int *pd, int start_value, int num_items);
void probability_distribution_free(probability_distribution_int *pd);

int probability_distribution_add(probability_distribution_int *pd, int value, int count);
int probability_distribution_remove(probability_distribution_int *pd, int value, int count);
void probability_distribution_print(probability_distribution_int *pd, const char *delimiter);
int probability_distribution_python_plot(char *file_name, probability_distribution_int *pd, char *file_name_prefix, const char *python_commmand);

#endif // FMM_PROBABILITY_DISTRIBUTION_UTILS
