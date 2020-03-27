#ifndef MYMATH_H
#define MYMATH_H

int find_index(int *ar, int len, int val);
int unique_count(int *ar, int len, int val);

void create_adjacency_matrix(int **mat, int *ar, int len);
void union_adjacency_matrix(int **um, int **am1, int **am2, int len);
void adjacency_matrix_remove_element(int **mat, int len, int val);
void print_adjacency_matrix(int **mat, int nrow, int ncol);
void calc_matrix_row_lengths(int **mat, int nrow, int ncol, int *len);

_Bool rand_bool();
int rand_int(int min, int max);
double rand_real(double min, double max);
void rand_perm(int *array, int length);
int rand_entry(int *ar, int len);
void shuffle_perm(int *array, int len);

void int_to_binary(unsigned int val, _Bool *ar, int len);
int binary_to_int(_Bool *ar, int len);
void binary_to_gray(_Bool *b, int len);
void gray_to_binary(_Bool *g, int len);
int find_binary_length(int val);

double choose(double n, double k);

#endif
