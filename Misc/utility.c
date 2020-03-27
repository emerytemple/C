#include "utility.h"

#include <math.h>

void swap(double* a, double* b)
{
	double temp = *b;
	*b = *a;
	*a = temp;
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
		swap(&a,&b);
	}

	c = a;
	mflag = 1;

	for(i = 1; (i < 100) && fabs(b-a) > delta; ++i) {
		fa = (*fp)(a);
		fb = (*fp)(b);
		fc = (*fp)(c);

		R = fb/fc;
		S = fb/fa;
		T = fa/fc;
		if((fa != fc) && (fb != fc)) {
			P = S*((T*(R-T)*(c-b))-((1.0-R)*(b-a)));
			Q = (T-1.0)*(R-1.0)*(S-1.0);
			s = b + (P/Q); // inverse quadratic interpolation
		}
		else {
			s = b - (S*((b-a)/(S-1.0))); // secant method
		}

		low = ((3.0*a)+b)/4.0;
		bflag = 0;
		if((low <= s) && (s <= b)) { // is s in interval
			if(mflag == 1) {
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
			mflag = 1;
		}
		else {
			// puts("	no conditions met");
			mflag = 0;
		}

		d = c;
		c = b;
		fa = (*fp)(a);
		fs = (*fp)(s);
		if(fa*fs < 0.0) {
			b = s;
		} else {
			a = s;
		}

		fa = (*fp)(a);
		fb = (*fp)(b);
		if(fabs(fa) < fabs(fb)) {
			swap(&a,&b);
		}
	}

	return s;
}

