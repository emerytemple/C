#include "ga.h"
#include "mymath.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

void check_gaparam(struct GAParam *p)
{
	int i;

	assert(0 <= p->max_iter);
	assert(0 <= p->pop_size);

	assert(0 <= p->num_bool);
	assert(0 <= p->num_int);
	for(i = 0; i < p->num_int; ++i)
		assert(p->imin[i] <= p->imax[i]);

	assert(0 <= p->num_real);
	for(i = 0; i < p->num_real; ++i)
		assert(p->rmin[i] <= p->rmax[i]);

	// TODO(emery): make sure imin and imax arrays
	// have length of num_int

	assert(0 <= p->num_perm);
	for(i = 0; i < p->num_perm; ++i)
		assert(0 <= p->perm_len);

	assert(0.0 <= p->pm && p->pm <= 1.0);
	assert(0 <= p->npts);
	assert(0 <= p->tsk);
	assert(0.0 <= p->alpha && p->alpha <= 1.0);
	assert(0.0 <= p->c);
}

void print_gaparam(struct GAParam *p)
{
	int i;

	printf("max iter = %i\n", p->max_iter);
	printf("pop size = %i\n", p->pop_size);
	printf("\n");

	printf("num bool = %u\n\n", p->num_bool);
	printf("num int = %i\n", p->num_int);
	for(i = 0; i < p->num_int; ++i)
		printf("%i	[%i, %i]\n", i+1, p->imin[i], p->imax[i]);
	printf("\n");
	printf("num real = %i\n", p->num_real);
	for(i = 0; i < p->num_real; ++i)
		printf("%i	[%f, %f]\n", i+1, p->rmin[i], p->rmax[i]);
	printf("\n");
	printf("num perm = %i\n", p->num_perm);
	for(i = 0; i < p->num_perm; ++i)
		printf("%i	%i\n", i+1, p->perm_len[i]);
	printf("\n");

	printf("pm = %f\n", p->pm);
	printf("npts = %i\n", p->npts);
	printf("tsk = %i\n", p->tsk);
	printf("alpha = %f\n", p->alpha);
	printf("\n");
}


void evolve(struct GAParam *p, double (*funcp)(struct Chromosome*, struct GAParam*))
{
	int i, j;

	srand(time(NULL));

	double sum, avg, stdev;

	struct Chromosome *max_chrome;
	struct Chromosome **chrome;

	check_gaparam(p);

	chrome = malloc(p->pop_size*sizeof(struct Chromosome*));
	for(i = 0; i < p->pop_size; ++i) {
		chrome[i] = new_chromosome(p);
		populate_chrome(chrome[i], p);
		chrome[i]->fitness = (*funcp)(chrome[i], p);
	}

	max_chrome = new_chromosome(p);
	populate_chrome(max_chrome, p);
	max_chrome->fitness = (*funcp)(max_chrome, p);

	elitism(max_chrome, chrome, p);

	calc_stats(&sum, &avg, &stdev, chrome, p);
	print_stats(0, sum, avg, stdev, max_chrome->fitness);

	for(i = 1; i <= p->max_iter; ++i)
	{
		parent_selection(max_chrome, chrome, p, avg, stdev);

		crossover(chrome, p);
		mutation(chrome, p, i);

		for(j = 0; j < p->pop_size; ++j)
			chrome[j]->fitness = (*funcp)(chrome[j], p);

		elitism(max_chrome, chrome, p);

		calc_stats(&sum, &avg, &stdev, chrome, p);
		print_stats(i, sum, avg, stdev, max_chrome->fitness);
	}
	printf("\n");
	print_chrome(max_chrome, p);

	for(i = 0; i < p->pop_size; ++i)
		delete_chrome(chrome[i], p);
	free(chrome);
	delete_chrome(max_chrome, p);
}

void parent_selection(struct Chromosome *max_chrome, struct Chromosome **chrome, struct GAParam *p, double avg, double stdev)
{
	int i;

	struct Chromosome **mating_pool;

	mating_pool = malloc(p->pop_size*sizeof(struct Chromosome*));
	for(i = 0; i < p->pop_size; ++i)
		mating_pool[i] = new_chromosome(p);

	copy_chrome(mating_pool[0], max_chrome, p);

	switch(p->parent_select) {
		case PS_FPS:
				fitness_proportional_selection(chrome, p, avg, stdev);
				selection_probability(mating_pool, chrome, p);
			break;
		case PS_RANKING:
				ranking_selection(chrome, p);
				selection_probability(mating_pool, chrome, p);
			break;
		case PS_TOURNAMENT:
			tournament_selection(mating_pool, chrome, p);
			break;
		case PS_TRUNCATION:
			truncation_selection(mating_pool, chrome, p);
			break;
		case PS_UNIFORM_CROSSOVER:
			hboa_probabilistic_uniform_crossover_selection(mating_pool, chrome, p);
			break;
	}

	for(i = 0; i < p->pop_size; ++i) {
		copy_chrome(chrome[i], mating_pool[i], p);
		delete_chrome(mating_pool[i], p);
	}
	free(mating_pool);
}

void crossover(struct Chromosome **chrome, struct GAParam *p)
{
	int i;

	for(i = 0; i < p->pop_size; i += 2) {
		switch(p->bool_cross) {
			case BC_UNIFORM:
				bool_uniform_crossover(chrome[i]->bool_chrome, chrome[i+1]->bool_chrome, p->num_bool, p->pm);
				break;
			case BC_NPOINT:
				bool_npoint_crossover(chrome[i]->bool_chrome, chrome[i+1]->bool_chrome, p->num_bool, p->npts);
				break;
		}

		// TODO(emery): clean up int crossover functions

		switch(p->int_cross) {
			case IC_UNIFORM:
				int_uniform_crossover(chrome[i]->int_chrome, chrome[i+1]->int_chrome, p->num_int, p->imin, p->imax, p->pm);
				break;
			case IC_NPOINT:
				int_npoint_crossover(chrome[i]->int_chrome, chrome[i+1]->int_chrome, p->num_int, p->imin, p->imax, p->npts);
				break;
		}

		switch(p->real_cross) {
			case RC_DISCRETE:
				real_discrete_crossover(chrome[i]->real_chrome, chrome[i+1]->real_chrome, p->num_real);
				break;
			case RC_SIMPLE:
				real_simple_crossover(chrome[i]->real_chrome, chrome[i+1]->real_chrome, p->num_real, p->alpha);
				break;
			case RC_SINGLE:
				real_single_crossover(chrome[i]->real_chrome, chrome[i+1]->real_chrome, p->num_real, p->alpha);
				break;
			case RC_WHOLE:
				real_whole_crossover(chrome[i]->real_chrome, chrome[i+1]->real_chrome, p->num_real, p->alpha);
				break;
		}

		switch(p->perm_cross) {
			case PC_CYCLE:
				perm_cycle_crossover(chrome[i]->perm_chrome, chrome[i+1]->perm_chrome, p->num_perm, p->perm_len);
				break;
			case PC_EDGE:
				perm_edge_crossover(chrome[i]->perm_chrome, chrome[i+1]->perm_chrome, p->num_perm, p->perm_len);
				break;
			case PC_ORDER:
				perm_order_crossover(chrome[i]->perm_chrome, chrome[i+1]->perm_chrome, p->num_perm, p->perm_len);
				break;
			case PC_PMX:
				perm_pmx_crossover(chrome[i]->perm_chrome, chrome[i+1]->perm_chrome, p->num_perm, p->perm_len);
				break;
		}
	}
}

void mutation(struct Chromosome **chrome, struct GAParam *p, int curr_iter)
{
	int i;

	for(i = 0; i < p->pop_size; ++i) {
		switch(p->bool_mutate) {
			case BM_BITWISE:
				bool_bitwise_mutation(chrome[i]->bool_chrome, p->num_bool, p->pm);
				break;
		}

		switch(p->int_mutate) {
			case IM_UNIFORM:
				int_uniform_mutation(chrome[i]->int_chrome, p->num_int, p->imin, p->imax, p->pm);
				break;
			case IM_NONUNIFORM:
				int_nonuniform_mutation(chrome[i]->int_chrome, p->num_int, p->imin, p->imax, p->pm, p->b, curr_iter, p->max_iter);
				break;
			case IM_BOUNDARY:
				int_boundary_mutation(chrome[i]->int_chrome, p->num_int, p->imin, p->imax, p->pm);
				break;
			case IM_CREEP:
				int_creep_mutation(chrome[i]->int_chrome, p->num_int, p->imin, p->imax, p->pm, p->s);
				break;
		}

		switch(p->real_mutate) {
			case RM_UNIFORM:
				real_uniform_mutation(chrome[i]->real_chrome, p->num_real, p->rmin, p->rmax, p->pm);
				break;
			case RM_NONUNIFORM:
				// real_nonuniform_mutation(chrome[i]->real_chrome, p->num_real, p->rmin, p->rmax, p->pm, p->b, curr_iter, p->max_iter);
				break;
			case RM_BOUNDARY:
				real_boundary_mutation(chrome[i]->real_chrome, p->num_real, p->rmin, p->rmax, p->pm);
				break;
			case RM_CREEP:
				// real_creep_mutation(chrome[i]->real_chrome, p->num_real, p->rmin, p->rmax, p->pm, p->s);
				break;
		}

		switch(p->perm_mutate) {
			case PM_INSERT:
				perm_insert_mutation(chrome[i]->perm_chrome, p->num_perm, p->perm_len, p->pm);
				break;
			case PM_INVERSION:
				perm_inversion_mutation(chrome[i]->perm_chrome, p->num_perm, p->perm_len, p->pm);
				break;
			case PM_SCRAMBLE:
				perm_scramble_mutation(chrome[i]->perm_chrome, p->num_perm, p->perm_len, p->pm);
				break;
			case PM_SWAP:
				perm_swap_mutation(chrome[i]->perm_chrome, p->num_perm, p->perm_len, p->pm);
				break;
		}
	}
}

struct Chromosome *new_chromosome(struct GAParam *p)
{
	int i;

	struct Chromosome *chrome;

	chrome = malloc(sizeof(struct Chromosome));

	chrome->bool_chrome = malloc(p->num_bool*sizeof(_Bool));
	chrome->int_chrome = malloc(p->num_int*sizeof(int));
	chrome->real_chrome = malloc(p->num_real*sizeof(double));
	chrome->perm_chrome = malloc(p->num_perm*sizeof(int*));
	for(i = 0; i < p->num_perm; ++i)
		chrome->perm_chrome[i] = malloc(p->perm_len[i]*sizeof(int));

	return chrome;
}

void populate_chrome(struct Chromosome *self, struct GAParam *p)
{
	int i;

	for(i = 0; i < p->num_bool; ++i)
		self->bool_chrome[i] = rand_bool();

	for(i = 0; i < p->num_int; ++i)
		self->int_chrome[i] = rand_int(p->imin[i], p->imax[i]);

	for(i = 0; i < p->num_real; ++i)
		self->real_chrome[i] = rand_real(p->rmin[i], p->rmax[i]);

	for(i = 0; i < p->num_perm; ++i)
		rand_perm(self->perm_chrome[i], p->perm_len[i]);
}


void print_chrome(struct Chromosome *self, struct GAParam *p)
{
	int i, j;

	for(i = 0; i < p->num_bool; ++i)
		printf("%i", self->bool_chrome[i]);
	printf("%s","	");

	for(i = 0; i < p->num_int; ++i)
		printf("%3i ", self->int_chrome[i]);
	printf("%s","	");

	for(i = 0; i < p->num_real; ++i)
		printf("%10.6f ", self->real_chrome[i]);
	printf("%s","	");

	for(i = 0; i < p->num_perm; ++i) {
		for(j = 0; j < p->perm_len[i]; ++j)
			printf("%i ", self->perm_chrome[i][j]);
		if(i != p->num_perm-1)
			printf(", ");
	}

	printf("%f\n", self->fitness);
}

void copy_chrome(struct Chromosome *to, struct Chromosome *from, struct GAParam *p)
{
	int i, j;

	for(i = 0; i < p->num_bool; ++i)
		to->bool_chrome[i] = from->bool_chrome[i];

	for(i = 0; i < p->num_int; ++i)
		to->int_chrome[i] = from->int_chrome[i];

	for(i = 0; i < p->num_real; ++i)
		to->real_chrome[i] = from->real_chrome[i];

	for(i = 0; i < p->num_perm; ++i)
		for(j = 0; j < p->perm_len[i]; ++j)
			to->perm_chrome[i][j] = from->perm_chrome[i][j];

	to->fitness = from->fitness;
}

void delete_chrome(struct Chromosome *self, struct GAParam *p)
{
	int i;

	free(self->bool_chrome);
	free(self->int_chrome);
	free(self->real_chrome);
	for(i = 0; i < p->num_perm; ++i)
		free(self->perm_chrome[i]);
	free(self->perm_chrome);

	free(self);
}

void calc_stats(double *sum, double *avg, double *stdev, struct Chromosome **chrome, struct GAParam *p)
{
	int i;

	*sum = 0.0;
	for(i = 0; i < p->pop_size; ++i)
		*sum += chrome[i]->fitness;
	*avg = *sum/p->pop_size;

	*stdev = 0.0;
	for(i = 0; i < p->pop_size; ++i)
		*stdev += pow(chrome[i]->fitness - *avg, 2.0);
	*stdev = sqrt(*stdev/p->pop_size);
}

void print_stats(int curr_iter, double sum, double avg, double stdev, double max)
{
	printf("%i	", curr_iter);
	printf("%10.6f	", sum);
	printf("%10.6f	", avg);
	printf("%10.6f	", stdev);
	printf("%10.6f\n", max);
}

void print_population(struct Chromosome *max_chrome, struct Chromosome **chrome, struct GAParam *p)
{
	int i;

	for(i = 0; i < p->pop_size; ++i)
		print_chrome(chrome[i], p);
	printf("\n");
	print_chrome(max_chrome, p);
	printf("\n");
}

void elitism(struct Chromosome *max_chrome, struct Chromosome **chrome, struct GAParam *p)
{
	int i;
	int flag = 0;

	for(i = 0; i < p->pop_size; ++i) {
		if(max_chrome->fitness <= chrome[i]->fitness) {
			copy_chrome(max_chrome, chrome[i], p);
			flag = 1;
		}
	}

	// NOTE(emery): if max_chrome has highest fitness
	// and is not in population

	if(flag == 0)
		copy_chrome(chrome[0], max_chrome, p);
}

void fitness_proportional_selection(struct Chromosome **chrome, struct GAParam *p, double avg, double stdev)
{
	int i;
	double scale;

	// Goldberg's sigma scaling
	scale = avg - (p->c*stdev);
	for(i = 0; i < p->pop_size; ++i)
		chrome[i]->fitness = fmax(chrome[i]->fitness - scale, 0.0);

	for(i = 1; i < p->pop_size; ++i)
		chrome[i]->fitness += chrome[i-1]->fitness;

	for(i = 0; i < p->pop_size; ++i)
		chrome[i]->fitness = chrome[i]->fitness/chrome[p->pop_size-1]->fitness;
}

void ranking_selection(struct Chromosome **chrome, struct GAParam *p)
{
	int i;
	int *rank;

	rank = (int*)malloc(p->pop_size*sizeof(int));

	for(i = 0; i < p->pop_size; ++i)
		rank[i] = i;

	insertion_sort(chrome, p, rank, 1); // sort based on fitness

	for(i = 0; i < p->pop_size; ++i)
		chrome[i]->fitness = ((2-p->s)/p->pop_size) + ((2*i*(p->s-1))/(p->pop_size*(p->pop_size-1)));

	insertion_sort(chrome, p, rank, 2); // sort based on chrome number

	for(i = 1; i < p->pop_size; ++i)
		chrome[i]->fitness += chrome[i-1]->fitness;

	free(rank);
}

void insertion_sort(struct Chromosome **chrome, struct GAParam *p, int *rank, int choice)
{
	int i;
	int hole_pos;
	int int_temp;
	double dbl_temp;

	// choice = 1 for sorting by fitness, choice = 2 for sorting by index

	for(i = 1; i < p->pop_size; ++i) {
		int_temp = rank[i];
		dbl_temp = chrome[i]->fitness;

		for(hole_pos = i; hole_pos > 0; --hole_pos) {
			if(1 == choice) {
				if(dbl_temp >= chrome[hole_pos-1]->fitness)
					break;
			} else if(2 == choice) {
				if(int_temp >= rank[hole_pos-1])
					break;
			}

			chrome[hole_pos]->fitness = chrome[hole_pos-1]->fitness;
			rank[hole_pos] = rank[hole_pos-1];
		}

		chrome[hole_pos]->fitness = dbl_temp;
		rank[hole_pos] = int_temp;
	}
}

void selection_probability(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p)
{
	int i, j;
	int ind;
	double random;

	switch(p->selection_prob) {
		case SP_ROULETTE:
			for(i = 1; i < p->pop_size; ++i) {

				ind = 0;
				random = rand_real(0.0, 1.0);
				for(j = 0; j < p->pop_size; ++j) {
					if(random <= chrome[j]->fitness) {
						ind = j;
						break;
					}
				}

				copy_chrome(mating_pool[i], chrome[ind], p);
			}
			break;
		case SP_SUS:
			random = rand_real(0.0, 1.0);
			for(i = 1; i < p->pop_size; ++i) {
				random /= (double)p->pop_size;
				random += 1.0/(double)p->pop_size;

				ind = 0;
				for(j = 0; j < p->pop_size; ++j) {
					if(random <= chrome[j]->fitness) {
						ind = j;
						break;
					}
				}

				copy_chrome(mating_pool[i], chrome[ind], p);
			}
			break;
	}
}

void tournament_selection(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p)
{
	int i, j;
	int random, tmp;

	// TODO(emery): maybe make max number in tournament
	// less than the population size

	for(i = 1; i < p->pop_size; ++i) {
		tmp = rand_int(0, p->pop_size-1);
		for(j = 1; j < p->tsk; ++j) {
			random = rand_int(0, p->pop_size-1);
			if(chrome[random]->fitness > chrome[tmp]->fitness)
				tmp = random;
		}

		copy_chrome(mating_pool[i], chrome[tmp], p);
	}
}

void truncation_selection(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p)
{
	int i, ind;
	int *rank;

	rank = (int*)malloc(p->pop_size*sizeof(int));

	for(i = 0; i < p->pop_size; ++i)
		rank[i] = i;

	insertion_sort(chrome, p, rank, 1); // sort based on fitness

	for(i = 0; i < p->pop_size; ++i)
		copy_chrome(mating_pool[i], chrome[p->pop_size-1-(i%p->tsk)], p);

	free(rank);
}

void hboa_probabilistic_uniform_crossover_selection(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p)
{
	int i, j;
	double *pm;
	double rand;

	pm = malloc(p->num_bool*sizeof(double));

	for(i = 0; i < p->pop_size; ++i)
		print_chrome(chrome[i], p);
	printf("\n");

	for(j = 0; j < p->num_bool; ++j) {
		pm[j] = 0.0;
		for(i = 0; i < p->pop_size; ++i) {
			if(1 == chrome[i]->bool_chrome[j]) {
				pm[j]++;
			}
		}
		pm[j] /= p->pop_size;
	}

	for(i = 0; i < p->pop_size; ++i) {
		copy_chrome(mating_pool[i], chrome[i], p);
		for(j = 0; j < p->num_bool; ++j) {
			if(rand_real(0.0, 1.0) <= pm[j]) {
				mating_pool[i]->bool_chrome[j] = 1;
			} else {
				mating_pool[i]->bool_chrome[j] = 0;
			}
		}
	}

	for(i = 0; i < p->pop_size; ++i)
		print_chrome(mating_pool[i], p);
	printf("\n");
}

void bool_bitwise_mutation(_Bool *chrome, int len, double pm)
{
	int i;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm)
			chrome[i] = !chrome[i];
	}
}


void int_uniform_mutation(int *chrome, int len, int *min, int *max, double pm)
{
	int i;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm)
			chrome[i] = rand_int(min[i], max[i]);
	}
}

void int_nonuniform_mutation(int *chrome, int len, int *min, int *max, double pm, double b, int curr_iter, int max_iter)
{
	int i;
	double x, frac, del;

	frac = (double)curr_iter/ (double)max_iter;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm) {
			x = chrome[i];
			del = (1.0-pow(rand_real(0.0, 1.0),pow(1.0-frac, b)));

			if(rand_real(0.0, 1.0) > 0.5)
				chrome[i] = (int)floor(x+(x-min[i])*del);
			else
				chrome[i] = (int)floor(x-(max[i]-x)*del);
		}
	}
}

void int_boundary_mutation(int *chrome, int len, int *min, int *max, double pm)
{
	int i;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm) {
			if(rand_real(0.0, 1.0) > 0.5)
				chrome[i] = max[i];
			else
				chrome[i] = min[i];
		}
	}
}

void int_creep_mutation(int *chrome, int len, int *min, int *max, double pm, double s)
{
	int i;
	double x, r;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm) {
			x = chrome[i];
			r = rand_real(-1.0,1.0);

			x = fmin(max[i], fmax(min[i], x+(r*s*(max[i]-min[i]))));
			chrome[i] = (int)floor(x);
		}
	}
}

void real_uniform_mutation(double *chrome, int len, double *min, double *max, double pm)
{
	int i;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm)
			chrome[i] = rand_real(min[i], max[i]);
	}
}

void dbl_nonuniform_mutation(double *chrome, int len, double *min, double *max, double pm, double b, int curr_iter, int max_iter)
{
	int i;
	double x, frac, del;

	frac = (double)curr_iter/ (double)max_iter;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm) {
			x = chrome[i];
			del = (1.0-pow(rand_real(0.0, 1.0),pow(1.0-frac, b)));

			if(rand_real(0.0, 1.0) > 0.5)
				chrome[i] = x+(x-min[i])*del;
			else
				chrome[i] = x-(max[i]-x)*del;
		}
	}
}

void real_boundary_mutation(double *chrome, int len, double *min, double *max, double pm)
{
	int i;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm) {
			if(rand_real(0.0, 1.0) > 0.5)
				chrome[i] = max[i];
			else
				chrome[i] = min[i];
		}
	}
}

void dbl_creep_mutation(double *chrome, int len, double *min, double *max, double pm, double s)
{
	int i;
	double x, r;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0, 1.0) < pm) {
			x = chrome[i];
			r = rand_real(-1.0,1.0);

			x = fmin(max[i], fmax(min[i], x+(r*s*(max[i]-min[i]))));
			chrome[i] = x;
		}
	}
}

void perm_insert_mutation(int **chrome, int num, int *len, double pm)
{
	int i, j;
	int p1, p2;
	int tmp;

	for(i = 0; i < num; ++i) {
		if((rand_real(0.0, 1.0) < pm) && (len[i] > 1)) {
			p1 = rand_int(0, len[i]-1);
			p2 = rand_int(1, len[i]-1);

			if(p1 == p2)
				p2 = 0;

			tmp = chrome[i][p2];
			if(p1 < p2) {
				for(j = p2; j > p1; --j)
					chrome[i][j] = chrome[i][j-1];
				chrome[i][p1+1] = tmp;
			} else { // p2 < p1
				for(j = p2; j < p1; ++j)
					chrome[i][j] = chrome[i][j+1];
				chrome[i][p1] = tmp;
			}
		}
	}
}

void perm_inversion_mutation(int **chrome, int num, int *len, double pm)
{
	int i, j;
	int p1, p2;
	int tmp, ind;

	for(i = 0; i < num; ++i) {
		if((rand_real(0.0, 1.0) < pm) && (len[i] > 1)) {
			// one iteration of knuth shuffle
			p1 = rand_int(0, len[i]-1);
			p2 = rand_int(1, len[i]-1);

			if(p1 == p2)
				p2 = 0;

			// find middle
			if(0 == ((p2 - p1) % 2))
				ind = p1 + ((p2-p1)/2);
			else
				ind = p1 + ((p2-p1)/2) + 1;

			// reverse the order
			for(j = p1; j < ind; ++j) {
				tmp = chrome[i][j];
				chrome[i][j] = chrome[i][p2+p1-j];
				chrome[i][p2+p1-j] = tmp;
			}
		}
	}
}

void perm_scramble_mutation(int **chrome, int num, int *len, double pm)
{
	int i, j;
	int p1, p2;
	int *array, arlen;

	for(i = 0; i < num; ++i) {
		if((rand_real(0.0, 1.0) < pm) && (len[i] > 1)) {
			p1 = rand_int(0, len[i]-1);
			p2 = rand_int(1, len[i]-1);

			if(p1 == p2)
				p2 = 0;

			arlen = p2-p1+1;
			array  = malloc(arlen*sizeof(int));
			for(j = 0; j < arlen; ++j)
				array[j] = chrome[i][p1+j];

			shuffle_perm(array, arlen);

			for(j = 0; j < arlen; ++j)
				chrome[i][p1+j] = array[j];

			free(array);
		}
	}
}

void perm_swap_mutation(int **chrome, int num, int *len, double pm)
{
	int i;
	int p1, p2;
	int tmp;

	for(i = 0; i < num; ++i) {
		if((rand_real(0.0, 1.0) < pm) && (len[i] > 1)) {
			p1 = rand_int(0, len[i]-1);
			p2 = rand_int(1, len[i]-1);

			if(p1 == p2)
				p2 = 0;

			tmp = chrome[i][p1];
			chrome[i][p1] = chrome[i][p2];
			chrome[i][p2] = tmp;
		}
	}
}

void bool_uniform_crossover(_Bool *chrome1, _Bool *chrome2, int len, double pm)
{
	_Bool tmp;
	int i;

	for(i = 0; i < len; ++i) {
		if(rand_real(0.0,1.0) < pm) {
			tmp = chrome1[i];
			chrome1[i] = chrome2[i];
			chrome2[i] = tmp;
		}
	}
}

void bool_npoint_crossover(_Bool *chrome1, _Bool *chrome2, int len, int npts)
{
	_Bool tmp;
	int i, j;

	// TODO(emery): maybe make crossover point unique
	// and limit the number to less than len

	for(i = 0; i < npts; ++i) {
		for(j = rand_int(0, len-1); j < len; ++j) {
			if(chrome1[j] != chrome2[j]) {
				tmp = chrome1[j];
				chrome1[j] = chrome2[j];
				chrome2[j] = tmp;
			}
		}
	}
}

void int_uniform_crossover(int *chrome1, int *chrome2, int num, int *min, int *max, double pm)
{
	_Bool *ar1, *ar2;
	int i;
	int int1, int2;
	int sgn1, sgn2;
	int len1, len2, len;

	for(i = 0; i < num; ++i) {
		len1 = find_binary_length(min[i]);
		len2 = find_binary_length(max[i]);
		len = (len1 > len2) ? len1: len2;

		ar1 = malloc(len*sizeof(_Bool));
		ar2 = malloc(len*sizeof(_Bool));

		int1 = chrome1[i];
		int2 = chrome2[i];

		sgn1 = (int1 < 0) ? -1: 1;
		sgn2 = (int2 < 0) ? -1: 1;

		int_to_binary(fabs(int1), ar1, len);
		int_to_binary(fabs(int2), ar2, len);

		binary_to_gray(ar1, len);
		binary_to_gray(ar2, len);

		bool_uniform_crossover(ar1, ar2, len, pm);

		gray_to_binary(ar1, len);
		gray_to_binary(ar2, len);

		int1 = binary_to_int(ar1, len);
		int2 = binary_to_int(ar2, len);

		int1 = copysign(int1, sgn1);
		int2 = copysign(int2, sgn2);

		int1 = (int1 > max[i]) ? max[i]: int1;
		int1 = (int1 < min[i]) ? min[i]: int1;

		int2 = (int2 > max[i]) ? max[i]: int2;
		int2 = (int2 < min[i]) ? min[i]: int2;

		chrome1[i] = int1;
		chrome2[i] = int2;

		free(ar1);
		free(ar2);
	}
}

void int_npoint_crossover(int *chrome1, int *chrome2, int num, int *min, int *max, int npts)
{
	_Bool *ar1, *ar2;
	int i;
	int int1, int2;
	int sgn1, sgn2;
	int len1, len2, len;

	for(i = 0; i < num; ++i) {
		len1 = find_binary_length(min[i]);
		len2 = find_binary_length(max[i]);
		len = (len1 > len2) ? len1: len2;

		ar1 = malloc(len*sizeof(_Bool));
		ar2 = malloc(len*sizeof(_Bool));

		int1 = chrome1[i];
		int2 = chrome2[i];

		sgn1 = (int1 < 0) ? -1: 1;
		sgn2 = (int2 < 0) ? -1: 1;

		int_to_binary(fabs(int1), ar1, len);
		int_to_binary(fabs(int2), ar2, len);

		binary_to_gray(ar1, len);
		binary_to_gray(ar2, len);

		bool_npoint_crossover(ar1, ar2, len, npts);

		gray_to_binary(ar1, len);
		gray_to_binary(ar2, len);

		int1 = binary_to_int(ar1, len);
		int2 = binary_to_int(ar2, len);

		int1 = copysign(int1, sgn1);
		int2 = copysign(int2, sgn2);

		int1 = (int1 > max[i]) ? max[i]: int1;
		int1 = (int1 < min[i]) ? min[i]: int1;

		int2 = (int2 > max[i]) ? max[i]: int2;
		int2 = (int2 < min[i]) ? min[i]: int2;

		chrome1[i] = int1;
		chrome2[i] = int2;

		free(ar1);
		free(ar2);
	}
}

void real_discrete_crossover(double *chrome1, double *chrome2, int len)
{
	int i;
	double tmp;

	for(i = rand_int(1, len-1); i < len; ++i) {
		tmp = chrome1[i];
		chrome1[i] = chrome2[i];
		chrome2[i] = tmp;
	}
}

void real_simple_crossover(double *chrome1, double *chrome2, int len, double alpha)
{
	int i;
	double tmp;

	for(i = rand_int(1, len-1); i < len; ++i) {
		tmp = chrome1[i];
		chrome1[i] = ((1.0-alpha)*tmp) + (alpha*chrome2[i]);
		chrome2[i] = (alpha*tmp) + ((1.0-alpha)*chrome2[i]);
	}
}

void real_single_crossover(double *chrome1, double *chrome2, int len, double alpha)
{
	int i;
	double tmp;

	i = rand_int(1, len-1);
	tmp = chrome1[i];
	chrome1[i] = ((1.0-alpha)*tmp) + (alpha*chrome2[i]);
	chrome2[i] = (alpha*tmp) + ((1.0-alpha)*chrome2[i]);
}

void real_whole_crossover(double *chrome1, double *chrome2, int len, double alpha)
{
	int i;
	double tmp;

	for(i = 0; i < len; ++i) {
		tmp = chrome1[i];
		chrome1[i] = ((1.0-alpha)*tmp) + (alpha*chrome2[i]);
		chrome2[i] = (alpha*tmp) + ((1.0-alpha)*chrome2[i]);
	}
}

void perm_cycle_crossover(int **chrome1, int **chrome2, int num, int *len)
{
	int i, j, k;
	int *c1;
	int inst, ind;
	int tmp;

	for(k = 0; k < num; ++k) {
		c1 = malloc(len[k]*sizeof(int));

		for(i = 0; i < len[k]; ++i)
			c1[i] = -1;

		inst = 0;
		for(i = 0; i < len[k]; ++i) {
			if(-1 == c1[i]) {
				ind = i;
				c1[i] = inst;
				for(j = 0; j < len[k]; ++j) {
					ind = find_index(chrome1[k], len[k], chrome2[k][ind]);

					if(ind == i) {
						inst++;
						break;
					}
					else
						c1[ind] = inst;
				}
			}
		}

		for(i = 0; i < len[k]; ++i) {
			if(1 == c1[i] % 2) {
				tmp = chrome1[k][i];
				chrome1[k][i] = chrome2[k][i];
				chrome2[k][i] = tmp;
			}
		}

		free(c1);
	}
}

void perm_edge_crossover(int **chrome1, int **chrome2, int num, int *len)
{
	int i, j, k;
	int curr, val, imin;
	int *offspring, *lenarr;
	int **al1, **al2;
	int **tab;
	int tmp[4] = {0,0,0,0};
	int tmplen;

	offspring = malloc(len[1]*sizeof(int));
	lenarr = malloc(len[1]*sizeof(int));

	al1 = malloc(len[1]*sizeof(int*));
	for(i = 0; i < len[1]; ++i)
		al1[i] = malloc(2*sizeof(int));

	al2 = malloc(len[1]*sizeof(int*));
	for(i = 0; i < len[1]; ++i)
		al2[i] = malloc(2*sizeof(int));

	tab = malloc(len[1]*sizeof(int*));
	for(i = 0; i < len[1]; ++i) {
		tab[i] = malloc(4*sizeof(int));
		for(j = 0; j < 4; ++j)
			tab[i][j] = -1;
	}

	// construct adjacency matrix
	create_adjacency_matrix(al1, chrome1[1], len[1]);
	create_adjacency_matrix(al2, chrome2[1], len[1]);
	union_adjacency_matrix(tab, al1, al2, len[1]);

	calc_matrix_row_lengths(tab, len[1], 4, lenarr);

	for(i = 0; i < len[1]; ++i) {
		if(i == 0)
			curr = rand_int(1, len[1]);
		else {
			// if common edge
			tmplen = 0;
			for(j = 0; j < 2; ++j) {
				for(k= 0; k < 2; ++k) {
					if(al1[curr-1][j] == al2[curr-1][k] && unique_count(offspring, i, al1[curr-1][j]) == 0) {
						tmp[tmplen++] = al1[curr-1][j];
					}
				}
			}
			if(tmplen > 0)
				curr = rand_entry(tmp, tmplen);

			// else if shortest length
			if(0 == tmplen) {
				imin = 0;
				for(j = 0; j < 4; ++j) {
					val = lenarr[tab[curr-1][j]-1];
					if(0 == val)
						break;
					else if(0 == j)
						imin = val;
					else if(val < imin)
						imin = val;
				}
				for(j = 0; j < 4; ++j) {
					val = lenarr[tab[curr-1][j]-1];
					if(0 == val)
						break;
					else if(val == imin)
						tmp[tmplen++] = tab[curr-1][j];
				}
				if(0 == imin)
					// else randomly chosen
					curr = rand_entry(tab[curr-1], lenarr[curr-1]);
				else
					curr = rand_entry(tmp, tmplen);		
			}
		}

		// remove curr from table
		adjacency_matrix_remove_element(tab, len[1], curr);
		calc_matrix_row_lengths(tab, len[1], 4, lenarr);

		// place in offspring
		offspring[i] = curr;
	}

	for(i = 0; i < len[1]; ++i) {
		chrome1[1][i] = offspring[i];
		chrome2[1][i] = offspring[i];
	}

	free(offspring);
	free(lenarr);

	for(i = 0; i < len[1]; ++i) {
		free(al1[i]);
		free(al2[i]);
		free(tab[i]);
	}
	free(al1);
	free(tab);
	free(al2);
}

void perm_order_crossover(int **chrome1, int **chrome2, int num, int *len)
{
	int i;
	int r1, r2;
	int *c1, *c2;
	int start, end;
	int ind, p1, p2;
	int val1, val2;
	int res1, res2;

	c1 = (int*)malloc(len[1]*sizeof(int));
	c2 = (int*)malloc(len[1]*sizeof(int));

	r1 = rand_int(0, len[1] - 2);
	r2 = rand_int(0, len[1] - 1);

	start = fmin(r1, r2);
	end = fmax(r1, r2);

	for(i = 0; i < len[1]; ++i) {
		if(start <= i && i <= end) {
			c1[i] = chrome1[1][i];
			c2[i] = chrome2[1][i];
		} else {
	 		c1[i] = 0;
			c2[i] = 0;
		}
	}

	p1 = end;
	p2 = end;
	for(i = 0; i < len[1]; ++i) {
		ind = (end+i+1) % len[1];

		val1 = chrome1[1][ind];
		val2 = chrome2[1][ind];

		res1 = unique_count(c1, len[1], val2);
		res2 = unique_count(c2, len[1], val1);

		if(0 == res1) {
			p1++;
			p1 %= len[1];
			c1[p1] = val2;
		}
		if(0 == res2) {
			p2++;
			p2 %= len[1];
			c2[p2] = val1;
		}
	}

	for(i = 0; i < len[1]; ++i) {
		chrome1[1][i] = c1[i];
		chrome2[1][i] = c2[i];
	}

	free(c1);
	free(c2);
}

void perm_pmx_crossover(int **chrome1, int **chrome2, int num, int *len)
{
	int i, j;
	int temp, index;
	int r1, r2;
	int start, end;
	int *c1, *c2;

	c1 = (int*)malloc(len[1]*sizeof(int));
	c2 = (int*)malloc(len[1]*sizeof(int));

	r1 = rand_int(0, len[1] - 2);
	r2 = rand_int(0, len[1] - 1);

	start = fmin(r1, r2);
	end = fmax(r1, r2);

	for(i = 0; i < len[1]; ++i) {
		if(start <= i && i <= end) {
			c1[i] = chrome1[1][i];
			c2[i] = chrome2[1][i];
		} else {
	 		c1[i] = 0;
			c2[i] = 0;
		}
	}

	for(j = start; j <= end; ++j) {
		temp = 1;
		for(i = start; i <= end; ++i)
			if(chrome2[1][j] == chrome1[1][i])
				temp = 0;
		if(0 != temp) {
			index = recursive(j, chrome1[1], chrome2[1], len[1], start, end);
			c1[index] = chrome2[1][j];
		}

		temp = 1;
		for(i = start; i <= end; ++i) {
			if(chrome1[1][j] == chrome2[1][i])
				temp = 0;
		}
		if(0 != temp) {
			index = recursive(j, chrome2[1], chrome1[1], len[1], start, end);
			c2[index] = chrome1[1][j];
		}
	}

	for(j = 0; j < len[1]; ++j) {
		if(0 == c1[j])
			c1[j] = chrome2[1][j];
		if(0 == c2[j])
			c2[j] = chrome1[1][j];
	}

	for(i = 0; i < len[1]; ++i) {
		chrome1[1][i] = c1[i];
		chrome2[1][i] = c2[i];
	}

	free(c1);
	free(c2);
}

int recursive(int curr, int *p1, int *p2, int len, int r1, int r2)
{
	int k;
	int retval = 1;
	int temp, temp2;
	
	if((curr >= r1) && (curr <= r2))
		temp = 0;
	else
		temp = 1;

	if(1 == temp)
		retval = curr;
	else {
		temp2 = -1;
		for(k = 0; k < len; ++k)
			if(p2[k] == p1[curr])
				temp2 = k;
		retval = recursive(temp2, p1, p2, len, r1, r2);
	}
	return retval;
}

