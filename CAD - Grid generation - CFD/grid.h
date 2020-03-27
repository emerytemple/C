#ifndef GRID_H
#define GRID_H

#include "geometry.h"

struct Point **transfinite_interpolation(struct Point *top, struct Point *bottom, struct Point *left, struct Point *right, int nx, int ny);
void elliptic(struct Point **grid, int nx, int ny);
void gnuplot_write(char *filename, struct Point **grid, int nx, int ny);

// TODO(emery): add TTM, Thomas Middlecoff, and grape control functions
// also add hyperbolic grid generation and stretching functions

#endif
