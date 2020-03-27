#include "mymath.h"

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <stdio.h>

int find_index(int *ar, int len, int val)
{
	int ind = 0;
	while(ind < len && ar[ind] != val)
		++ind;

	return(ind == len ? -1 : ind);
}

int unique_count(int *ar, int len, int val)
{
	int i;
	int cnt = 0;

	for(i = 0; i < len; ++i)
		if(ar[i] == val)
			cnt++;

	return cnt;
}

void create_adjacency_matrix(int **mat, int *ar, int len)
{
	int i;
	int loc, ind1, ind2;

	for(i = 0; i < len; ++i) {
		loc = find_index(ar, len, i+1);
		if(0 == loc) {
			ind1 = len-1;
			ind2 = loc+1;
		} else if(loc == len-1) {
			ind1 = loc-1;
			ind2 = 0;
		} else {
			ind1 = loc-1;
			ind2 = loc+1;
		}

		if(ar[ind1] > ar[ind2]) {
			mat[i][0] = ar[ind2];
			mat[i][1] = ar[ind1];
		} else {
			mat[i][0] = ar[ind1];
			mat[i][1] = ar[ind2];
		}
	}

}

void union_adjacency_matrix(int **um, int **am1, int **am2, int len)
{
	int i, j;
	int tmp;
	int *array;

	array = malloc(4*sizeof(int));

	for(i = 0; i < len; ++i) {
		// sort values
		array[0] = (am1[i][0] < am2[i][0]) ? am1[i][0] : am2[i][0];
		array[1] = (am1[i][0]< am2[i][0]) ? am2[i][0] : am1[i][0];
		array[2] = (am1[i][1] > am2[i][1]) ? am2[i][1] : am1[i][1];
		array[3] = (am1[i][1] > am2[i][1]) ? am1[i][1] : am2[i][1];

		if(array[1] > array[2]) {
			tmp = array[1];
			array[1] = array[2];
			array[2] = tmp;
		}

		// update array
		tmp = 0;
		um[i][tmp++] = array[0];
		for(j = 1; j < 4; ++j)
			if(array[j] != array[j-1])
				um[i][tmp++] = array[j];
	}

	free(array);
}

void adjacency_matrix_remove_element(int **mat, int len, int val)
{
	int j, k;
	int found;

	for(j = 0; j < len; ++j) {
		found = 0;
		for(k = 0; k < 4; ++k) {
			if(mat[j][k] == val)
				found = 1;
			if(1 == found)
					mat[j][k] = mat[j][k+1];
		}
		if(1 == found)
			mat[j][3] = -1;
	}
}

void print_adjacency_matrix(int **mat, int nrow, int ncol)
{
	int i, j;

	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol; ++j) {
			printf("%i  ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void calc_matrix_row_lengths(int **mat, int nrow, int ncol, int *len)
{
	int i, j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol; ++j)
			if(-1 == mat[i][j])
				break;
		len[i] = j;
	}
}

_Bool rand_bool()
{
	return rand() % 2;
}

int rand_int(int min, int max) {
	return rand() % (max-min+1) + min;
}

double rand_real(double min, double max) {
	double range = max - min;
	double drand = (double)rand()/(double)RAND_MAX;

	return drand*range + min; // old school
}

void rand_perm(int *array, int len)
{
	int i;

	for(i = 0; i < len; ++i)
		array[i] = i+1;

	shuffle_perm(array, len);
}

int rand_entry(int *ar, int len)
{
	return ar[rand_int(0,len-1)];
}

void shuffle_perm(int *array, int len) // knuth shuffle
{
	int i, j;
	int tmp;

	for(i = len-1; i > 0; --i) {
		j = rand() % i;

		tmp = array[j];
		array[j] = array[i];
		array[i] = tmp;
	}
}

void int_to_binary(unsigned int val, _Bool *ar, int len)
{
	int k;

	// val = fabs(val);

	for(k = len - 1; k >= 0; --k) {
		if(val > 0) {
			ar[k] = (_Bool)(val % 2);
			val /= 2;
		}
		else
			ar[k] = 0;
	}
}

int binary_to_int(_Bool *ar, int len)
{
	int i, bval;
	int retval = 0;

	for(i = 0; i < len; ++i) {
		bval = (0 == ar[len - 1 - i]) ? 0 : 1;
		retval += (int)pow(2, i)*bval;
	}

	return retval;
}

void binary_to_gray(_Bool *b, int len)
{
	int i;
	_Bool g[len];

	g[0] = b[0];
	for(i = 1; i < len; ++i) {
		if(1 == b[i-1])
			g[i] = !b[i];
		else
			g[i] = b[i];
	}

	for(i = 0; i < len; ++i)
		b[i] = g[i];
}

void gray_to_binary(_Bool *g, int len)
{
	int i;
	_Bool b[len];

	b[0] = g[0];

	for(i = 1; i < len; ++i) {
		b[i] = g[i] ^ b[i-1];
	}

	for(i = 0; i < len; ++i)
		g[i] = b[i];
}

int find_binary_length(int val)
{
	int result;

	if(0 == val)
		result = 0;
	else
		result = (int)floor(log(fabs(val))/log(2.0))+1;

	return result;
}

// binomial coefficient, combination
double choose(double n, double k)
{
	int i;
	double r = 1;

	if(k == 0)
		return 1;

	k = fmin(k, n-k);
	for(i = 0; i < k; ++i)
		r *= (n-i)/(i+1.0);

	return r;
}


























