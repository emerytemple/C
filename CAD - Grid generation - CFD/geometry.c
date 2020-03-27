#include "geometry.h"
#include "mymath.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void set_point(struct Point *p, double x, double y)
{
	p->x = x;
	p->y = y;
}

void print_point(struct Point *p)
{
	printf("P = (%f, %f)\n", p->x, p->y);
}

void set_line(struct Line *line, struct Point *p1, struct Point *p2)
{
	line->p1 = p1;
	line->p2 = p2;
}

void print_line(struct Line *line)
{
	printf("start = (%f, %f),  ", line->p1->x, line->p1->y);
	printf("end = (%f, %f)\n", line->p2->x, line->p2->y);
}

struct Point *new_linear_interpolation(struct Line *line, int npts)
{
	int i;
	double vx, vy;
	struct Point *interp;

	interp = malloc(npts*sizeof(struct Point));

	vx = line->p2->x - line->p1->x;
	vy = line->p2->y - line->p1->y;

	interp[0].x = line->p1->x;
	interp[0].y = line->p1->y;

	for(i = 1; i < npts-1; ++i) {
		interp[i].x = line->p1->x + (i*vx/(npts-1.0));
		interp[i].y = line->p1->y + (i*vy/(npts-1.0));
	}

	interp[npts-1].x = line->p2->x;
	interp[npts-1].y = line->p2->y;
/*
	for(i = 0; i < npts; ++i) {
		print_point(&interp[i]);
	}
	puts("\n");
*/
	return interp; // does this work ??? ?!?!?
}

struct Point *combine_interpolation(struct Point *interp1, struct Point *interp2, int npts1, int npts2)
{
	int i, total;
	struct Point *retval;

	total = npts1 + npts2 - 1;
//	printf("total = %d\n", total);

	retval = malloc(total*sizeof(struct Point));

	for(i = 0; i < total; ++i) {
		if(i < npts1) {
			retval[i] = interp1[i];
//			printf("%d\n", i);
		} else {
			retval[i] = interp2[i-npts1+1];
//			printf("%d\n", i);
		}
//		print_point(&retval[i]);
	}

	return retval;
}

void print_interpolation(struct Point *interp, int npts)
{
	int i;

	for(i = 0; i < npts; ++i) {
		print_point(&interp[i]);
	}
}

void delete_interpolation(struct Point *interp)
{
	free(interp);
}

void linear_interpolation(struct Line line, int npts)
{
	int i;
	double curr;
	double len;

	line.interp = malloc(npts*sizeof(struct Point));

	for(i = 0; i < npts; ++i) {
		curr = (double)i/(double)(npts-1.0);
		len = length(line.use[0], line.use[1]);

		printf("curr = %f\n", curr);


		if(line.use[0].x < line.use[1].x) {
			line.interp[i].x = line.use[0].x + (curr*len);
		} else {
			line.interp[i].x = line.use[1].x + (curr*len);
		}

		if(line.use[0].y < line.use[1].y) {
			line.interp[i].y = line.use[0].y + (curr*len);
		} else {
			line.interp[i].y = line.use[1].y + (curr*len);
		}
	}
}

#if 0
void print_ctrl_points(struct Point **P, int len)
{
	int i;

	for(i = 0; i < len; ++i)
		printf("P[%i] = (%f, %f)\n", i, P[i]->x, P[i]->y);
	printf("\n");
}

void print_nurbs_curve(struct Point **P, int degree, int num_pt)
{
	int i;

	struct Point **C;

	C = malloc(num_pt*sizeof(struct Point*));

	for(i = 0; i < num_pt; ++i)
		C[i] = new_point(0.0, 0.0);

	for(i = 0; i < num_pt; ++i) {
		decasteljau(C[i], P, degree, i/(num_pt-1.0));
		print_point(C[i]);
	}

	for(i = 0; i < num_pt; ++i)
		free(C[i]);
	free(C);
}

void find_knot(char *choice, double *u, struct Point **P, int num_pt)
{
	int i;
	double el[num_pt];
	double sum;

	if(0 == strcmp(choice, "uniform")) {

		for(i = 0; i < num_pt; ++i)
			u[i] = (double)i/(num_pt-1.0);

	} else if(0 == strcmp(choice, "chordlength")) {

		sum = 0.0;
		for(i = 0; i < num_pt-1; ++i) {
			el[i] = length(P[i+1], P[i]);
			sum += el[i];
		}

		u[0] = 0.0;
		for(i = 1; i < num_pt; ++i)
			u[i] = u[i-1] + (el[i-1]/sum);

	} else if(0 == strcmp(choice, "centripetal")) {

		sum = 0.0;
		for(i = 0; i < num_pt-1; ++i) {
			el[i] = sqrt(length(P[i+1], P[i]));
			sum += el[i];
		}
		u[0] = 0.0;
		for(i = 1; i < num_pt; ++i)
			u[i] = u[i-1] + (el[i-1]/sum);
	}
}

void bezier(struct Point *C, struct Point **P, int n, double u)
{
	int i;
	double B[n];
	double temp1, temp2;

	for(i = 0; i < n; ++i)
		B[i] = choose(n-1.0, i)*pow(u, i)*pow(1.0-u, n-i-1.0);

	temp1 = 0.0;
	temp2 = 0.0;
	for(i = 0; i < n; ++i) {
		temp1 += B[i]*P[i]->x;
		temp2 += B[i]*P[i]->y;
	}

	C->x = temp1;
	C->y = temp2;
}

void decasteljau(struct Point *C, struct Point **P, int n, double u)
{
	int i, j;
	double Bx[n][n], By[n][n];

	for(i = 0; i < n; ++i) {
		for(j = 0; j < n; ++j) {
			Bx[i][j] = 0.0;
			By[i][j] = 0.0;
		}
	}

	for(i = 0; i < n; ++i) {
		Bx[0][i] = P[i]->x;
		By[0][i] = P[i]->y;
	}

	for(i = 1; i < n; ++i) {
		for(j = 0; j < n-i; ++j) {
			Bx[i][j] = ((1-u)*Bx[i-1][j]) + (u*Bx[i-1][j+1]);
			By[i][j] = ((1-u)*By[i-1][j]) + (u*By[i-1][j+1]);
		}
	}

	C->x = Bx[n-1][0];
	C->y = By[n-1][0];
}

void horner(struct Point *C, struct Point **P, int n, double u)
{
	int i;
	double Ctempx[n], Ctempy[n];
	double s;

	s = 1.0 - u;

	Ctempx[0] = P[0]->x;
	Ctempy[0] = P[0]->y;

	for(i = 1; i < n; ++i) {
		Ctempx[i] = (Ctempx[i-1]*s) + (choose(n-1, i)*P[i]->x*pow(u, i));
		Ctempy[i] = (Ctempy[i-1]*s) + (choose(n-1, i)*P[i]->y*pow(u, i));
	}

	C->x = Ctempx[n-1];
	C->y = Ctempy[n-1];
}

struct Point **degree_raise(struct Point **P, int *num_pt)
{
	int i;
	double a;

	// realloc can suck it
	struct Point **Pnew = malloc((*num_pt+1)*sizeof(struct Point*));
	for(i = 0; i < *num_pt+1; ++i)
		Pnew[i] = new_point(0.0, 0.0);

	Pnew[0]->x = P[0]->x;
	Pnew[0]->y = P[0]->y;

	for(i = 1; i < *num_pt; ++i) {
		a = (double)i / (double)*num_pt;

		Pnew[i]->x = (a*P[i-1]->x) + ((1.0-a)*P[i]->x);
		Pnew[i]->y = (a*P[i-1]->y) + ((1.0-a)*P[i]->y);
	}

	Pnew[*num_pt]->x = P[*num_pt-1]->x;
	Pnew[*num_pt]->y = P[*num_pt-1]->y;

	(*num_pt)++;

	for(i = 0; i < *num_pt-1; ++i)
		free(P[i]);
	free(P);

	return Pnew;
}

struct Point **degree_lower(struct Point **P, int *num_pt)
{
	int i;
//	double x1, x2, y1, y2;

	struct Point **Pnew = malloc((*num_pt+1)*sizeof(struct Point*));
	for(i = 0; i < *num_pt+1; ++i)
		Pnew[i] = new_point(0.0, 0.0);

	Pnew[0]->x = P[0]->x;
	Pnew[0]->y = P[0]->y;

	Pnew[*num_pt-2]->x = P[*num_pt-1]->x;
	Pnew[*num_pt-2]->y = P[*num_pt-1]->y;

	for(i = 1; i < ((*num_pt)-1.0)/2.0; ++i) {
		Pnew[i]->x = (((*num_pt-1)*P[i]->x) - ((i-1)*Pnew[i-1]->x)) / (*num_pt-i-1);
		Pnew[i]->y = (((*num_pt-1)*P[i]->y) - ((i-1)*Pnew[i-1]->y)) / (*num_pt-i-1);
	}
	for(i = (*num_pt)-2; i > (((*num_pt)-1.0)/2.0) - 1.0; --i) {
		Pnew[i-1]->x = (((*num_pt-1)*P[i]->x) - ((*num_pt-i-1)*Pnew[i]->x)) / i;
		Pnew[i-1]->y = (((*num_pt-1)*P[i]->y) - ((*num_pt-i-1)*Pnew[i]->y)) / i;
	}

/*
	if(0 == (*num_pt) % 2) {
		i = *num_pt/2;
		printf("iii = %i\n", i);
		x1 = (((*num_pt-1)*P[i]->x) - ((i-1)*Pnew[i-1]->x)) / (*num_pt-i-1);
		x2 = (((*num_pt-1)*P[i+1]->x) - ((*num_pt-i-1)*Pnew[i+1]->x)) / i;
		Pnew[i]->x = (x1+x2)/2.0;

		y1 = (((*num_pt-1)*P[i]->y) - ((i-1)*Pnew[i-1]->y)) / (*num_pt-i-1);
		y2 = (((*num_pt-1)*P[i+1]->y) - ((*num_pt-i-1)*Pnew[i+1]->y)) / i;
		Pnew[i]->y = (x1+x2)/2.0;
		printf("newx = %f, newy = %f\n", Pnew[i]->x, Pnew[i]->y);
	}
*/

	(*num_pt)--;

	for(i = 0; i < *num_pt-1; ++i)
		free(P[i]);
	free(P);

	return Pnew;
}

void rotate(double theta, struct Point **P, int n)
// angle is in degrees
// + is ccw, - is cw
{
	int i;
	double C, S;
	double Px, Py;
	double tmp;

	theta *= M_PI/180.0;

	C = cos(theta);
	S = sin(theta);

	Px = P[0]->x;
	Py = P[0]->y;

	for(i = 1; i < n; ++i) {
		tmp = P[i]->x;
		P[i]->x = tmp*C - P[i]->y*S + Px - Px*C + Py*S;
		P[i]->y = tmp*S + P[i]->y*C + Py - Px*S - Py*C;
	}
}

void scale(struct Point *o, double sx, double sy, struct Point **P, int n)
{
	int i;

	for(i = 0; i < n; ++i) {
		P[i]->x = (sx*P[i]->x) + ((1.0-sx)*o->x);
		P[i]->y = (sy*P[i]->y) + ((1.0-sy)*o->y);
	}
}

void translate(double x, double y, struct Point **P, int n)
{
	int i;

	for(i = 0; i < n; ++i) {
		P[i]->x += x;
		P[i]->y += y;
	}
}

struct Point **cubic_curve(struct Point *P1, struct Point *P2, double theta1, double theta2, double scale1, double scale2)
// angles are in degrees
{
	int i;
	int num_pt = 4;
	double len;

	struct Point **Pnew = malloc(num_pt*sizeof(struct Point*));
	for(i = 0; i < num_pt; ++i)
		Pnew[i] = new_point(0.0, 0.0);

	theta1 *= M_PI/180.0;
	theta2 *= M_PI/180.0;

	Pnew[0]->x = P1->x;
	Pnew[0]->y = P1->y;

	len = sqrt(pow(P2->x - P1->x, 2.0)+pow(P2->y - P1->y, 2.0));

	Pnew[1]->x = P1->x + (scale1*(len/3.0)*cos(theta1));
	Pnew[1]->y = P1->y + (scale1*(len/3.0)*sin(theta2));

	Pnew[2]->x = P2->x - (scale2*(len/3.0)*cos(theta1));
	Pnew[2]->y = P2->y - (scale2*(len/3.0)*sin(theta2));

	Pnew[num_pt-1]->x = P2->x;
	Pnew[num_pt-1]->y = P2->y;

	return Pnew;
}

void c1quadratic_interpolation(struct Point **Pnew, int newlen, struct Point **P, int n)
{
	int i;
	int k = 3; // because quadratic

	double p1, p2;
	double u[n+k], del[n+k-1];

	newlen = n+1;

	struct Point **PL = malloc((n+1)*sizeof(struct Point *));
	for(i = 0; i < n+1; ++i)
		PL[i] = new_point(0.0, 0.0);

	struct Point **PR = malloc((n+1)*sizeof(struct Point *));
	for(i = 0; i < n+1; ++i)
		PR[i] = new_point(0.0, 0.0);

	// uniform
	for(i = k-1; i < n+1; ++i)
		u[i] = (i+1.0-k)/(n+1.0-k);

	for(i = n+1; i < n+k; ++i)
		u[i] = u[i-1] + u[n] - u[n-1];

	for(i = k-2; i >= 0; --i)
		u[i] = u[i+1] - (u[k] - u[k-1]);

	for(i = 0; i < n+k-1; ++i)
		del[i] = u[i+1]-u[i];

	PL[0]->x = P[0]->x - ((1.0/3.0)*(P[1]->x - P[0]->x));
	PL[0]->y = P[0]->y - ((1.0/3.0)*(P[1]->y - P[0]->y));

	for(i = 1; i < n+1; ++i) {
		p1 = (del[i] + del[i+1])/del[i];
		p2 = del[i+1]/del[i];

		PL[i]->x = (p1*P[i-1]->x) - (p2*PL[i-1]->x);
		PL[i]->y = (p1*P[i-1]->y) - (p2*PL[i-1]->y);
	}

	PR[n]->x = P[n-1]->x + ((1.0/3.0)*(P[n-1]->x - P[n-2]->x));
	PR[n]->y = P[n-1]->y + ((1.0/3.0)*(P[n-1]->y - P[n-2]->y));

	for(i = n-1; i >= 0; --i) {
		p1 = (del[i+1] + del[i+2])/del[i+2];
		p2 = del[i+1]/del[i+2];

		PR[i]->x = (p1*P[i]->x) - (p2*PR[i+1]->x);
		PR[i]->y = (p1*P[i]->y) - (p2*PR[i+1]->y);
	}

	for(i = 0; i < n+1; ++i) {
		Pnew[i]->x = (PL[i]->x + PR[i]->x)/2.0;
		Pnew[i]->y = (PL[i]->y + PR[i]->y)/2.0;
	}

	for(i = 0; i < n+1; ++i)
		free(PL[i]);
	free(P);

	for(i = 0; i < n+1; ++i)
		free(PR[i]);
	free(PR);

	// return Pnew;
}

#endif
















