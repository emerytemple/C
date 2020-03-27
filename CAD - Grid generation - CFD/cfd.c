#include "cfd.h"

#include <math.h>
#include <stdio.h>

void cfd(struct Point **grid, int ni, int nj)
{
	int i, j, k;

	double Mi, pi, ri, ti, ai, ui, vi, ei;
	double Rgas, gamma;
	double gm1, tmg;
	double u[100][41][10][4], f[100][41][10][4], g[100][41][10][4];
	double u1[100][41][10][4], f1[100][41][10][4], g1[100][41][10][4];
	double ubar[100][41][10][4];
	int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2;
	double xs[100][41], xn[100][41], ys[100][41], yn[100][41];
	double J[100][41];
	double dt; // timestep
	double uavg;


	// constants
	Rgas = 287.0; // J/kg.K
	gamma = 1.4;
	gm1 = gamma-1.0;
	tmg = 3.0-gamma;

	// initial conditions
	Mi = 2.0;
	pi = 101000.0; // Pa
	ti = 286.0; // K

	dt = 0.00001;

	ai = sqrt(gamma*Rgas*ti);
	ui = Mi*ai;
	vi = 0.0;

	ri = pi/(Rgas*ti);
	ei = (Rgas*ti)/(gamma-1.0);

	// grid metrics
	for(i = 0; i < ni; ++i) {
		im2 = i-2;
		im1 = i-1;
		ip1 = i+1;
		ip2 = i+2;

		for(j = 0; j < nj; ++j) {
			jm2 = j-2;
			jm1 = j-1;
			jp1 = j+1;
			jp2 = j+2;

			if(0 == i) {
				xs[i][j] = (-(3.0*grid[i][j].x)+(4.0*grid[ip1][j].x)-grid[ip2][j].x)/2.0;
				ys[i][j] = (-(3.0*grid[i][j].y)+(4.0*grid[ip1][j].y)-grid[ip2][j].y)/2.0;
			} else if(ni-1 == i) {
				xs[i][j] = ((3.0*grid[i][j].x)-(4.0*grid[im1][j].x)+grid[im2][j].x)/2.0;
				ys[i][j] = ((3.0*grid[i][j].y)-(4.0*grid[im1][j].y)+grid[im2][j].y)/2.0;
			else if(0 == j) {
				xn[i][j] = (-(3.0*grid[i][j].x)+(4.0*grid[i][jp1].x)-grid[i][jp2].x)/2.0;
				yn[i][j] = (-(3.0*grid[i][j].y)+(4.0*grid[i][jp1].y)-grid[i][jp2].y)/2.0;
			} else if(nj-1 == j) {
				xn[i][j] = ((3.0*grid[i][j].x)-(4.0*grid[i][jm1].x)+grid[i][jm2].x)/2.0;
				yn[i][j] = ((3.0*grid[i][j].y)-(4.0*grid[i][jm1].y)+grid[i][jm2].y)/2.0;
			} else {
				xs[i][j] = (grid[ip1][j].x - grid[im1][j].x)/2.0;
				ys[i][j] = (grid[ip1][j].y - grid[im1][j].y)/2.0;

				xn[i][j] = (grid[i][jp1].x - grid[i][jm1].x)/2.0;
				yn[i][j] = (grid[i][jp1].y - grid[i][jm1].y)/2.0;
			}

			J[i][j] = (xs[i][j]*yn[i][j])-(ys[i][j]*xn[i][j]);
		}
	}

	// initialize flowfield
	for(i = 0; i < ni; ++i) {
		for(j = 0; j < nj; ++j) {
			u[i][j][0][0] = ri;
			u[i][j][0][1] = ri*ui;
			u[i][j][0][2] = ri*vi;
			u[i][j][0][3] = ri*(ei+0.5*((ui*ui)+(vi*vi)));

			f[i][j][0][0] = u[i][j][0][1];
			f[i][j][0][1] = (gm1*u[i][j][0][3])+(tmg*(u[i][j][0][1]*u[i][j][0][1])/(2.0*u[i][j][0][0]))-(gm1*((u[i][j][0][2]*u[i][j][0][2])/(2.0*u[i][j][0][0])));
			f[i][j][0][2] = (u[i][j][0][1]*u[i][j][0][2])/u[i][j][0][0];
			f[i][j][0][3] = (gamma*((u[i][j][0][1]*u[i][j][0][3])/u[i][j][0][0]))-(gm1*((u[i][j][0][1]*((u[i][j][0][1]*u[i][j][0][1])+(u[i][j][0][2]*u[i][j][0][2])))/(2.0*u[i][j][0][0]*u[i][j][0][0])));

			g[i][j][0][0] = u[i][j][0][2];
			g[i][j][0][1] = (u[i][j][0][1]*u[i][j][0][2])/u[i][j][0][0];
			g[i][j][0][2] = (gm1*u[i][j][0][3])+(tmg*(u[i][j][0][2]*u[i][j][0][2])/(2.0*u[i][j][0][0]))-(gm1*((u[i][j][0][1]*u[i][j][0][1])/(2.0*u[i][j][0][0])));
			g[i][j][0][3] = (gamma*((u[i][j][0][1]*u[i][j][0][3])/u[i][j][0][0]))-(gm1*((u[i][j][0][1]*((u[i][j][0][1]*u[i][j][0][1])+(u[i][j][0][2]*u[i][j][0][2])))/(2.0*u[i][j][0][0]*u[i][j][0][0])));

			for(k = 0; k < 4; ++k) {
				u1[i][j][0][k] = J[i][j]*u[i][j][0][k];
				f1[i][j][0][k] = (f[i][j][0][k]*yn[i][j])-(g[i][j][0][k]*xn[i][j]);
				g1[i][j][0][k] = (g[i][j][0][k]*xs[i][j])-(f[i][j][0][k]*ys[i][j]);
			}
		}
	}

	for(n = 0; n > 10; ++n) {
		for(i = 0; i < ni; ++i) {
			for(j = 0; j < nj; ++j) {
				for(k = 0; k < 4; ++k) {
					ubar[i][j][n+1][k] = u[i][j][n][k] - (dt*(f1[i+1][j][n][k]-f1[i][j][n][k])) - (dt*(g1[i][j+1][n][k]-g1[i][j][n][k]));
					uavg = 0.5*(ubar[i][j][n+1][k]+u[i][j][n][k]);
					u[i][j][n+1][k] = uavg - (dt*(f1[i][j][n][k]-f1[i-1][j][n][k])) - (dt*(g1[i][j][n][k]-g1[i][j-1][n][k]));
				}
			}
		}
	}
}


















