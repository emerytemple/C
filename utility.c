#include "utility.h"

#include <math.h>

void swap(double* a, double* b)
{
	// printf("a = %f, b = %f\n", *a, *b);
	double temp = *b;
	*b = *a;
	*a = temp;
	// printf("a = %f, b = %f\n", *a, *b);
}

double brentq(double (*fp)(double), double a, double b)
{
	int mflag, bflag, i;
	double fa, fb, fc, fs, c, d, delta, s;
	double P, Q, R, S, T;
	double low;

	delta = 0.000001;

	fa = (*fp)(a);
	fb = (*fp)(b);

	if(fa*fb >= 0.0) {
		return 0; // exit because root is not bracketed
	}

	if(fabs(fa) < fabs(fb)) {
		// puts("swap");
		swap(&a,&b);
	}

	c = a;
	mflag = 1;

	for(i = 1; (i < 100) && fabs(b-a) > delta; ++i) {
		fa = (*fp)(a);
		fb = (*fp)(b);
		fc = (*fp)(c);

		/*
		printf("i = %d\n", i);
		if(i != 1) {
			printf("	a%d = %f, b%d = %f, b%d = %f, b%d = %f\n", i-1, a, i-1, b, i-2, c, i-3, d);
			printf("	fa = %f, fb = %f, fc = %f\n", fa, fb, fc);
		}
		else {
			printf("	a%d = %f, b%d = %f, b%d = %f\n", i-1, a, i-1, b, i-2, c);
			printf("	fa = %f, fb = %f, fc = %f\n", fa, fb, fc);
		}
		*/

		R = fb/fc;
		S = fb/fa;
		T = fa/fc;
		if((fa != fc) && (fb != fc)) {
			// puts("	inv quad");
			P = S*((T*(R-T)*(c-b))-((1.0-R)*(b-a)));
			Q = (T-1.0)*(R-1.0)*(S-1.0);
			s = b + (P/Q); // inverse quadratic interpolation
			// printf("	s = %f\n", s);
		}
		else {
			// puts("	secant");
			s = b - (S*((b-a)/(S-1.0))); // secant method
			// printf("	s = %f\n", s);
		}

		/*
		printf("	mflag = %d\n", mflag);
		printf("	1	%f <= %f <= %f\n", ((3.0*a)+b)/4.0, s, b);
		if(mflag == 1) {
			printf("	2	mflag = 1, %f < %f\n", fabs(s-b), (fabs(b-c)/2.0));
			printf("	4	mflag = 1, %f < %f\n", fabs(delta), fabs(b-c));
		}
		else {
			printf("	3	mflag = 0, %f < %f\n", fabs(s-b), (fabs(c-d)/2.0));
			printf("	5	mflag = 0, %f < %f\n", fabs(delta), fabs(c-d));
		}
		*/

		low = ((3.0*a)+b)/4.0;
		bflag = 0;
		if((low <= s) && (s <= b)) { // is s in interval
			if(mflag == 1) { // is mflag set
				if(fabs(s-b) >= fabs(b-c)/2.0) {
					// puts("	condition 2");
					bflag = 1;
				}
				if(fabs(b-c) < fabs(delta)) {
					// puts("	condition 4");
					bflag = 1;
				}
			}
			else {
				if(fabs(s-b) >= fabs(c-d)/2.0) {
					// puts("	condition 3");
					bflag = 1;
				}
				if(fabs(c-d) < fabs(delta)) {
					// puts("	condition 5");
					bflag = 1;
				}
			}
		}
		else {
			// puts("	condition 1");
			bflag = 1;
		}

		if(bflag == 1) {
			// puts("	bisection");
			s = (a+b)/2.0; // bisection method
			// printf("	s = %f\n", s);
			mflag = 1;
		}
		else {
			// puts("	no conditions met");
			mflag = 0;
		}

		d = c;
		c = b;
		// printf("	a = %f, b = %f\n	fa = %f, fs = %f\n", a, s, fff(a), fff(s));
		fa = (*fp)(a);
		fs = (*fp)(s);
		if(fa*fs < 0.0) {
			b = s;
			// puts("	b = s");
		}
		else {
			a = s;
			// puts("	a = s");
		}

		// printf("	d = %f, fs = %f\n", d, fff(s));
		// printf("	fa = %f, fb = %f\n", fabs(fff(a)), fabs(fff(b)));
		fa = (*fp)(a);
		fb = (*fp)(b);
		if(fabs(fa) < fabs(fb)) {
			// puts("swap");
			swap(&a,&b);
		}
	}

	return s;
}

