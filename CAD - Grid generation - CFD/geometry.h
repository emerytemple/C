#ifndef GEOMETRY_H
#define GEOMETRY_H

struct Point {
	union {
		double p[2];
		struct {
			double x;
			double y;
		};
	};
};

struct Line {
	union {
		struct Point *l[2];
		struct {
			struct Point *p1;
			struct Point *p2;
		};
	};
};

void set_point(struct Point *p, double x, double y);
void print_point(struct Point *p);

void set_line(struct Line *line, struct Point *p1, struct Point *p2);
void print_line(struct Line *line);

struct Point *new_linear_interpolation(struct Line *line, int npts);
struct Point *combine_interpolation(struct Point *interp1, struct Point *interp2, int npts1, int npts2);
void print_interpolation(struct Point *interp, int npts);
void delete_interpolation(struct Point *interp);

/*
struct NURBS {
	struct Point *cp;
	struct Point uses[2];

	// int degree; ??? do i need
};

void print_ctrl_pts(struct Point **P, int len);
void print_nurbs_curve(struct Point **P, int degree, int num_pt);

void find_knot(char *choice, double *u, struct Point **P, int num_pt);  // do I need?

void bezier(struct Point *C, struct Point **P, int n, double u);
void decasteljau(struct Point *C, struct Point **P, int n, double u);
void horner(struct Point *C, struct Point **P, int n, double u);

struct Point **degree_raise(struct Point **P, int *num_pt);
struct Point **degree_lower(struct Point **P, int *num_pt);

void rotate(double angle, struct Point **P, int n);
void scale(struct Point *o, double sx, double sy, struct Point **P, int n);
void translate(double x, double y, struct Point **P, int n);

struct Point **cubic_curve(struct Point *P1, struct Point *P2, double theta1, double theta2, double scale1, double scale2);

void c1quadratic_interpolation(struct Point **Pnew, int newlen, struct Point **P, int n);
*/
#endif
