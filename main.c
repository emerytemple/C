#include "geometry.h"
#include "grid.h"
#include "cfd.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sor();

int main(int argc, char **argv)
{
/*
	double ni, nj;

	struct Point tl, tr, bl, bm, br;
	struct Line top, left, right, botleft, botright;
	struct Point *top2, *left2, *right2, *botleft2, *botright2, *bot2;
	struct Point **grid;

	// create cad geometry
	set_point(&tl, 0.0, 40.0);
	set_point(&tr, 65.0, 40.0);
	set_point(&bl, 0.0, 0.0);
	set_point(&bm, 10.0, 0.0);
	set_point(&br, 65.0, -9.697984);

	set_line(&top, &tl, &tr);
	set_line(&left, &bl, &tl);
	set_line(&right, &br, &tr);
	set_line(&botleft, &bl, &bm);
	set_line(&botright, &bm, &br);

	ni = 100;
	nj = 41;

	top2 = new_linear_interpolation(&top, ni);
	left2 = new_linear_interpolation(&left, nj);
	right2 = new_linear_interpolation(&right, nj);
	botleft2 = new_linear_interpolation(&botleft, 20);
	botright2 = new_linear_interpolation(&botright, 81);

	bot2 = combine_interpolation(botleft2, botright2, 20, 81);

	// todo: free linear interpolations

	// generate grid
	grid = transfinite_interpolation(top2, bot2, left2, right2, ni, nj);
	gnuplot_write("mygrid.dat", grid, ni, nj);
*/
	sor();

//	cfd(grid, ni, nj);

	return 0;
}

void sor()
{
	int i, j, k;
	int n, N;

	double w, tol;
	double a[4][4], b[4];
	double x[10][4];

	double sum1, sum2;
	double norm;

	n = 4;
	tol = 0.00001;
	N = 11;
	w = 1.0;

	a[0][0] = 10.0;
	a[0][1] = -1.0;
	a[0][2] = 2.0;
	a[0][3] = 0.0;
	a[1][0] = -1.0;
	a[1][1] = 11.0;
	a[1][2] = -1.0;
	a[1][3] = 3.0;
	a[2][0] = 2.0;
	a[2][1] = -1.0;
	a[2][2] = 10.0;
	a[2][3] = -1.0;
	a[3][0] = 0.0;
	a[3][1] = 3.0;
	a[3][2] = -1.0;
	a[3][3] = 8.0;

	x[0][0] = 0.0;
	x[0][1] = 0.0;
	x[0][2] = 0.0;
	x[0][3] = 0.0;

	b[0] = 6.0;
	b[1] = 25.0;
	b[2] = -11.0;
	b[3] = 15.0;

	printf("%d  %f  %f  %f  %f\n", 0, x[0][0], x[0][1], x[0][2], x[0][3]);

	for(k = 1; k < N; ++k) {
		for(i = 0; i < n; ++i) {
			sum1 = 0.0;
			for(j = 0; j < i; ++j) {
				sum1 += a[i][j]*x[k][j];
			}

			sum2 = 0.0;
			for(j = i+1; j < n; ++j) {
				sum2 += a[i][j]*x[k-1][j];
			}

			x[k][i] = ((1.0-w)*x[k-1][i]) + ((w*(-sum1-sum2+b[i])) / a[i][i]);
		}

		norm = 0.0;
		for(i = 0; i < n; ++i) {
			norm += pow(x[k][i]-x[k-1][i], 2);
		}
		norm = sqrt(norm);

		printf("%d,  %f  %f  %f  %f,  %f\n", k, x[k][0], x[k][1], x[k][2], x[k][3], norm);

		if(norm < tol) {
			break;
		}
	}
}
























