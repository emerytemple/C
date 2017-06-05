#include "grid.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Point **transfinite_interpolation(struct Point *top, struct Point *bottom, struct Point *left, struct Point *right, int nx, int ny)
{
	int i, j, k;
	double u, v;
	double F1, F2;

	struct Point **grid;

	grid = malloc(nx*sizeof(struct Point*));
	for(i = 0; i < nx; ++i)
		grid[i] = malloc(ny*sizeof(struct Point));

	for(i = 0; i < nx; ++i) {
		u =  (double)i/(double)(nx-1);
		for(j = 0; j < ny; ++j) {
			v = (double)j/(double)(ny-1);
			for(k = 0; k < 2; ++k) {
				F1 = ((1.0-v)*bottom[i].p[k]) + (v*top[i].p[k]) + ((1.0-u)*left[j].p[k]) + (u*right[j].p[k]);
				F2 = ((1.0-u)*(1.0-v)*bottom[0].p[k]) + (u*v*top[nx-1].p[k]) + (u*(1.0-v)*bottom[nx-1].p[k]) + ((1.0-u)*v*top[0].p[k]);

				grid[i][j].p[k] = F1 - F2;
			}
		}
	}

	return grid;
}

void elliptic(struct Point **grid, int nx, int ny)
// NOTE(emery): winslow smoothing, control functions are zero
{
	int i, j, k, it;
	int max_iter;
	int im1, ip1, jm1, jp1;
	double xs, ys, xn, yn;
	double xss, yss, xsn, ysn, xnn, ynn;
	double a, b, g;
	double p1, p2, p3, bot;
	double tol, rms;
	double xnew, ynew, xtol, ytol;

	tol = 0.001;
	max_iter = 100;

	for(it = 0; it < max_iter; ++it) {
		for(i = 1; i < nx-1; ++i) {
			im1 = i-1;
			ip1 = i+1;
			for(j = 1; j < ny-1; ++j) {
				jm1 = j-1;
				jp1 = j+1;

				// central differences
				xs = 0.5*(grid[ip1][j].x - grid[im1][j].x);
				ys = 0.5*(grid[ip1][j].y - grid[im1][j].y);

				xn = 0.5*(grid[i][jp1].x - grid[i][jm1].x);
				yn = 0.5*(grid[i][jp1].y - grid[i][jm1].y);

				xss = grid[ip1][j].x - (2.0*grid[i][j].x) + grid[im1][j].x;
				yss = grid[ip1][j].y - (2.0*grid[i][j].y) + grid[im1][j].y;

				xsn = 0.25*(grid[ip1][jp1].x - grid[ip1][jm1].x - grid[im1][jp1].x + grid[im1][jm1].x);
				ysn = 0.25*(grid[ip1][jp1].y - grid[ip1][jm1].y - grid[im1][jp1].y + grid[im1][jm1].y);

				xnn = grid[i][jp1].x - (2.0*grid[i][j].x) + grid[i][jm1].x;
				ynn = grid[i][jp1].y - (2.0*grid[i][j].y) + grid[i][jm1].y;

				a = (xn*xn) + (yn*yn);
				b = (xs*xn) + (ys*yn);
				g = (xs*xs) + (ys*ys);

				// (a*rss) - (2.0*b*rsn) + (g*rnn) = 0
				// solve for rij
				p1 = a*(grid[ip1][j].x + grid[im1][j].x);
				p2 = 2.0*b*xsn;
				p3 = g*(grid[i][jp1].x + grid[i][jm1].x);
				bot = 2.0*(a + g);
				xnew = (p1 - p2 + p3) / bot;

				p1 = a*(grid[ip1][j].y + grid[im1][j].y);
				p2 = 2.0*b*xsn;
				p3 = g*(grid[i][jp1].y + grid[i][jm1].y);
				bot = 2.0*(a + g);
				ynew = (p1 - p2 + p3) / bot;

				xtol = fabs(grid[i][j].x - xnew);
				ytol = fabs(grid[i][j].y - ynew);

				rms = rms + (xtol*xtol) + (ytol*ytol);

				grid[i][j].x = xnew;
				grid[i][j].y = ynew;
			}
		}
		rms = sqrt(rms/nx/ny/2.0);
		if(rms < tol)
			break;
	}
}

void gnuplot_write(char *filename, struct Point **grid, int nx, int ny)
{
	// gnuplot
	// plot 'filename' using 1:2
	int i, j;

	FILE *fp;

	fp = fopen(filename,"w");

	for(i = 0; i < nx; ++i) {
		for(j = 0; j < ny; ++j) {
			fprintf(fp, "%f	%f\n", grid[i][j].x, grid[i][j].y);
		}
	}

	fclose(fp);
}
