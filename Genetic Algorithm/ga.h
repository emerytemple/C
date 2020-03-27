#ifndef GA_H
#define GA_H

struct Chromosome {
	_Bool *bool_chrome;
	int *int_chrome;
	double *real_chrome;
	int **perm_chrome;

	double fitness;
};

struct GAParam {
	int max_iter;
	int pop_size;

	int num_bool, num_int, num_real, num_perm;

	int *imin, *imax;
	double *rmin, *rmax;
	int *perm_len;

	enum {PS_FPS, PS_RANKING, PS_TOURNAMENT, PS_TRUNCATION, PS_UNIFORM_CROSSOVER} parent_select;
	enum {SP_ROULETTE, SP_SUS} selection_prob;

	enum {BC_UNIFORM, BC_NPOINT} bool_cross;
	enum {IC_UNIFORM, IC_NPOINT} int_cross;
	enum {RC_DISCRETE, RC_SIMPLE, RC_SINGLE, RC_WHOLE} real_cross;
	enum {PC_CYCLE, PC_EDGE, PC_ORDER, PC_PMX} perm_cross;

	enum {BM_BITWISE} bool_mutate;
	enum {IM_UNIFORM, IM_NONUNIFORM, IM_BOUNDARY, IM_CREEP} int_mutate;
	enum {RM_UNIFORM, RM_NONUNIFORM, RM_BOUNDARY, RM_CREEP} real_mutate;
	enum {PM_INSERT, PM_INVERSION, PM_SCRAMBLE, PM_SWAP} perm_mutate;

	double pm; // mutation rate (0.1 - 0.3)
	double pc; // crossover rate - add this functionality (0.7 - 0.9)
	int npts;
	int tsk;
	double alpha;
	double c;
	double s;
	double b;
};

// NOTE(emery): main functions
void evolve(struct GAParam *p, double (*funcp)(struct Chromosome*, struct GAParam*));
void elitism(struct Chromosome *max_chrome, struct Chromosome **chrome, struct GAParam *p);
void parent_selection(struct Chromosome *max_chrome, struct Chromosome **chrome, struct GAParam *p, double avg, double stdev);
void crossover(struct Chromosome **chrome, struct GAParam *p);
void mutation(struct Chromosome **chrome, struct GAParam *p, int curr_iter);
void calc_stats(double *sum, double *avg, double *stdev, struct Chromosome **chrome, struct GAParam *p);
void print_stats(int curr_iter, double sum, double avg, double max, double stdev);
void check_gaparam(struct GAParam *p);
void print_gaparam(struct GAParam *p);

// NOTE(emery): chromosome functions
struct Chromosome *new_chromosome(struct GAParam *p);
void populate_chrome(struct Chromosome *self, struct GAParam *p);
void print_chrome(struct Chromosome *self, struct GAParam *p);
void copy_chrome(struct Chromosome *to, struct Chromosome *from, struct GAParam *p);
void delete_chrome(struct Chromosome *self, struct GAParam *p);

// NOTE(emery): parent selection functions
void fitness_proportional_selection(struct Chromosome **chrome, struct GAParam *p, double avg, double stdev);
void ranking_selection(struct Chromosome **chrome, struct GAParam *p);
void insertion_sort(struct Chromosome **chrome, struct GAParam *p, int *rank, int choice); // helper for ranking selection
void selection_probability(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p);
void tournament_selection(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p);
void truncation_selection(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p);
void hboa_probabilistic_uniform_crossover_selection(struct Chromosome **mating_pool, struct Chromosome **chrome, struct GAParam *p);

// NOTE(emery): mutation functions
void bool_bitwise_mutation(_Bool *chrome, int len, double pm);

void int_uniform_mutation(int *chrome, int len, int *min, int *max, double pm);
void int_nonuniform_mutation(int *chrome, int len, int *min, int *max, double pm, double b, int curr_iter, int max_iter);
void int_boundary_mutation(int *chrome, int len, int *min, int *max, double pm);
void int_creep_mutation(int *chrome, int len, int *min, int *max, double pm, double s);

void real_uniform_mutation(double *chrome, int len, double *min, double *max, double pm);
void real_nonuniform_mutation(double *chrome, int len, double *min, double *max, double pm, double b, int curr_iter, int max_iter);
void real_boundary_mutation(double *chrome, int len, double *min, double *max, double pm);
void real_creep_mutation(double *chrome, int len, double *min, double *max, double pm, double s);

void perm_insert_mutation(int **chrome, int num, int *len, double pm);
void perm_inversion_mutation(int **chrome, int num, int *len, double pm);
void perm_scramble_mutation(int **chrome, int num, int *len, double pm);
void perm_swap_mutation(int **chrome, int num, int *len, double pm);

// NOTE(emery): crossover functions
void bool_uniform_crossover(_Bool *chrome1, _Bool *chrome2, int len, double pm);
void bool_npoint_crossover(_Bool *chrome1, _Bool *chrome2, int len, int npts);

void int_uniform_crossover(int *chrome1, int *chrome2, int num, int *min, int *max, double pm);
void int_npoint_crossover(int *chrome1, int *chrome2, int num, int *min, int *max, int npts);

void real_discrete_crossover(double *chrome1, double *chrome2, int len);
void real_simple_crossover(double *chrome1, double *chrome2, int len, double alpha);
void real_single_crossover(double *chrome1, double *chrome2, int len, double alpha);
void real_whole_crossover(double *chrome1, double *chrome2, int len, double alpha);

void perm_cycle_crossover(int **chrome1, int **chrome2, int num, int *len);

// TODO(emery): make these last three perm crossover functions
// act on all the perms, not just one

void perm_edge_crossover(int **chrome1, int **chrome2, int num, int *len);
void perm_order_crossover(int **chrome1, int **chrome2, int num, int *len);
void perm_pmx_crossover(int **chrome1, int **chrome2, int num, int *len);
int recursive(int curr, int *p1, int *p2, int len, int r1, int r2); // helper for pmx

#endif
