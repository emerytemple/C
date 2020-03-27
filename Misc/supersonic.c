#include "supersonic.h"
#include "utility.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

double Supersonic(char *string, ...)
{
	int i;
	int count;
	va_list list;

	char str;
	double val;
	double retval;

	double ToT;
	double Mach, gamma, gp1, gm1, PR;

	va_start(list, string);

	count = 2;
	for(i = 0; i < count; ++i) {
		str = va_arg(list, int); // char passed in is promoted to int the cast back to char
		val = va_arg(list, double);
		printf("%c = %f\n", str, val);
		switch(str) {
			case 'M':
				Mach = val;
				break;
			case 'g':
				gamma = val;
				gp1 = gamma+1.0;
				gm1 = gamma-1.0;
				break;
			case 'r':
				PR = val;
				break;
			default:
				break;
		}
	}
	printf("M = %f, gamma = %f\n", Mach, gamma);

	va_end(list);

	// 1-D
	if(strcmp(string, "P*/Po") == 0) {
		retval = pow(2.0/(gamma+1.0), gamma/(gamma-1.0));
	}
	else if(strcmp(string, "r*/ro") == 0) {
		retval = pow(2.0/(gamma+1.0), 1.0/(gamma-1.0));
	}
	else if(strcmp(string, "T*/To") == 0) {
		retval = 2.0/(gamma+1.0);
	}

	else if(strcmp(string, "Po/P") == 0) {
		ToT = 1.0+(((gamma-1.0)/2.0)*pow(Mach, 2.0));
		retval = pow(ToT, gamma/(gamma-1.0));
	}
	else if(strcmp(string, "ro/r") == 0) {
		ToT = 1.0+(((gamma-1.0)/2.0)*pow(Mach, 2.0));
		retval = pow(ToT, 1.0/(gamma-1.0));
	}
	else if(strcmp(string, "To/T") == 0) {
		retval = 1.0+(((gamma-1.0)/2.0)*pow(Mach, 2.0));
	}

	else if(strcmp(string, "Mpr") == 0) {
		retval = sqrt((2.0/(gamma-1.0))*(pow(PR,(gamma-1.0)/gamma)-1.0));
	}

	else {
	}

	return retval;

}

#if 0
// static variables with file scope (global variables)
double gxp, gg, gpr, gth, gb, gnu;

// helper functions

double ss_hMpr(double x)
{
	return gpr - ss_PoP(x,gg);
}

double ss_hMsub(double x)
{
	return gxp - ss_AAs(x,gg);
}

double ss_hMob(double x)
{
	return (tan(gth) - (2.0*(1.0/tan(gb))*((pow((x*sin(gb)), 2.0)-1.0)/((pow(x, 2.0)*(gg+cos(2.0*gb)))+2.0))));
}

double ss_hMnu(double x)
{
	double gp1gm1, gm1gp1;

	gp1gm1 = (gg+1.0)/(gg-1.0);
	gm1gp1 = (gg-1.0)/(gg+1.0);
	return ((sqrt(gp1gm1)*atan(sqrt(gm1gp1*(pow(x, 2.0)-1.0))))-atan(sqrt(pow(x, 2.0)-1.0))-gnu);
}

// 1-D
double ss_TsTo(double gamma)
{
	return 2.0/(gamma+1.0);
}

double ss_PsPo(double gamma)
{
	return pow(ss_TsTo(gamma), gamma/(gamma-1.0));
}

double ss_RsRo(double gamma)
{
	return pow(ss_TsTo(gamma), 1.0/(gamma-1.0));
}

double ss_ToT(double Mach, double gamma)
{
	return 1.0+(((gamma-1.0)/2.0)*pow(Mach, 2.0));
}

double ss_PoP(double Mach, double gamma)
{
	return pow(ss_ToT(Mach,gamma), gamma/(gamma-1.0));
}

double ss_RoR(double Mach, double gamma)
{
	return pow(ss_ToT(Mach,gamma), 1.0/(gamma-1.0));
}

double ss_Mpr(double PR, double gamma)
{
	gpr = PR;
	gg = gamma;

	return brentq(ss_hMpr, 0.0, 50.0);
}

// 1-D w/ heat addition
double ss_PPsh(double Mach, double gamma)
{
	return (1.0+gamma)/(1.0+(gamma*pow(Mach, 2.0)));
}

double ss_TTsh(double Mach, double gamma)
{
	return pow(Mach, 2.0)*pow(ss_PPsh(Mach, gamma), 2.0);
}

double ss_RRsh(double Mach, double gamma)
{
	return (1.0/pow(Mach, 2.0))*(1.0/ss_PPsh(Mach, gamma));
}

double ss_PoPosh(double Mach, double gamma)
{
	double base = (2.0+((gamma-1.0)*pow(Mach, 2.0)))/(gamma+1.0);

	return ss_PPsh(Mach, gamma)*pow(base, gamma/(gamma-1.0));
}

double ss_ToTosh(double Mach, double gamma)
{
	double top, bot;

	top = (gamma+1.0)*pow(Mach, 2.0);
	bot = pow(1.0+(gamma*pow(Mach, 2.0)), 2.0);

	return (top/bot)*(2.0+((gamma-1.0)*pow(Mach, 2.0)));
}

// 1-D w/ friction
double ss_TTsf(double Mach, double gamma)
{
	return (gamma+1.0)/(2.0+((gamma-1.0)*pow(Mach, 2.0)));
}

double ss_PPsf(double Mach, double gamma)
{
	return (1.0/Mach)*sqrt(ss_TTsf(Mach, gamma));
}

double ss_RRsf(double Mach, double gamma)
{
	return (1.0/Mach)*sqrt(1.0/ss_TTsf(Mach, gamma));
}

double ss_PoPosf(double Mach, double gamma)
{
	double exp = (gamma+1.0)/(2.0*(gamma-1.0));

	return (1.0/Mach)*pow(1.0/ss_TTsf(Mach, gamma), exp);
}

double ss_flD(double Mach, double gamma)
{
	double first, second, inside, top, bot;
	first = 6.5;
	top = 1.0-pow(Mach, 2);
	bot = gamma*pow(Mach, 2.0);
	second = (gamma+1.0)/(2.0*gamma);
	inside = ss_TTsf(Mach, gamma)*pow(Mach, 2.0);

	return (top/bot)+(second*log(inside));
}

// quasi 1-D
double ss_AAs(double Mach, double gamma)
{
	double exp, inside;
	exp = (gamma+1.0)/(gamma-1.0);
	inside = 2.0/(gamma+1.0);
	return sqrt(pow(inside*ss_ToT(Mach, gamma), exp)/pow(Mach, 2.0));
}

double ss_Msub(double xp, double gamma)
{
	gg = gamma;
	gxp = xp;

	return brentq(ss_hMsub, 0.0, 1.0);
}

double ss_Msup(double xp, double gamma)
{
	gg = gamma;
	gxp = xp;

	return brentq(ss_hMsub, 1.0, 50.0);
}

// normal shock
double ss_M2(double Mach, double gamma)
{
	return sqrt(ss_ToT(Mach,gamma)/((gamma*pow(Mach, 2.0))-((gamma-1.0)/2.0)));
}

double ss_P2P1(double Mach, double gamma)
{
	return 1.0+(((2.0*gamma)/(gamma+1.0))*(pow(Mach, 2.0)-1.0));
}

double ss_R2R1(double Mach, double gamma)
{
	return ((gamma+1.0)*pow(Mach, 2.0))/(2.0+((gamma-1.0)*pow(Mach, 2.0)));
}

double ss_T2T1(double Mach, double gamma)
{
	return ss_P2P1(Mach, gamma)/ss_R2R1(Mach, gamma);
}

double ss_Po2Po1(double Mach, double gamma, double R)
{
	double cp, ds;
	cp = (gamma*R)/(gamma-1.0);
	ds = (cp*log(ss_T2T1(Mach, gamma)))-(R*log(ss_P2P1(Mach, gamma)));
	return exp(-ds/R);
}

double ss_Po2P1(double Mach, double gamma, double R)
{
	return ss_Po2Po1(Mach, gamma, R)*pow(ss_ToT(Mach, gamma), gamma/(gamma-1.0));
}

// oblique shock
double ss_theta(double Mach, double gamma, double beta)
{
	double tcb, top, bot;

	beta *= (M_PI/180.0);

	tcb = 2.0*(1.0/tan(beta));
	top = pow((Mach*sin(beta)), 2.0)-1.0;
	bot = (pow(Mach, 2.0)*(gamma+cos(2.0*beta)))+2.0;

	return atan(tcb*(top/bot))*(180.0/M_PI);
}

double ss_betas(double Mach, double gamma, double theta)
{
	double m2, m3, gm2, gp2, gp4, t2, lam, xi, delta, cx, part, cp, tb, gm1, gp1;

	theta *= (M_PI/180.0);

	gm1 = gamma-1.0;
	gp1 = gamma+1.0;
	m2 = pow((pow(Mach, 2.0)-1.0), 2.0);
	m3 = pow((pow(Mach, 2.0)-1.0), 3.0);
	gm2 = (gm1/2.0)*pow(Mach, 2.0);
	gp2 = (gp1/2.0)*pow(Mach, 2.0);
	gp4 = (gp1/4.0)*pow(Mach, 4.0);
	t2 = pow(tan(theta), 2.0);
	lam = sqrt(m2-(3.0*(1.0+gm2)*(1.0+gp2)*t2));
	xi = (m3-(9.0*(1.0+gm2)*(1.0+gm2+gp4)*t2))/pow(lam, 3.0);
	delta = 0.0;
	cx = acos(xi);
	part = ((4.0*M_PI*delta)+cx)/3.0;
	cp = cos(part);
	tb = (pow(Mach, 2.0)-1.0+(2.0*lam*cp))/(3.0*(1.0+gm2)*tan(theta));

	return atan(tb)*(180.0/M_PI);
}

double ss_betaw(double Mach, double gamma, double theta)
{
	double m2, m3, gm2, gp2, gp4, t2, lam, xi, delta, cx, part, cp, tb, gm1, gp1;

	theta *= (M_PI/180.0);

	gm1 = gamma-1.0;
	gp1 = gamma+1.0;
	m2 = pow((pow(Mach, 2.0)-1.0), 2.0);
	m3 = pow((pow(Mach, 2.0)-1.0), 3.0);
	gm2 = (gm1/2.0)*pow(Mach, 2.0);
	gp2 = (gp1/2.0)*pow(Mach, 2.0);
	gp4 = (gp1/4.0)*pow(Mach, 4.0);
	t2 = pow(tan(theta), 2.0);
	lam = sqrt(m2-(3.0*(1.0+gm2)*(1.0+gp2)*t2));
	xi = (m3-(9.0*(1.0+gm2)*(1.0+gm2+gp4)*t2))/pow(lam, 3.0);
	delta = 1.0;
	cx = acos(xi);
	part = ((4.0*M_PI*delta)+cx)/3.0;
	cp = cos(part);
	tb = (pow(Mach, 2.0)-1.0+(2.0*lam*cp))/(3.0*(1.0+gm2)*tan(theta));

	return atan(tb)*(180.0/M_PI);
}

double ss_Mob(double theta, double beta, double gamma)
{
	gth = theta*(M_PI/180.0);
	gb = beta*(M_PI/180.0);
	gg = gamma;

	return brentq(ss_hMob, 1.0, 50.0);;
}

// prandtl-meyer expansion waves
double ss_mu(double Mach)
{
	return asin(1.0/Mach)*(180.0/M_PI);
}

double ss_nu(double Mach, double gamma)
{
	double m2, gm1, gp1, nu;

	gp1 = gamma+1.0;
	gm1 = gamma-1.0;

	m2 = pow(Mach, 2.0)-1.0;
	nu = (sqrt(gp1/gm1)*atan(sqrt((gm1/gp1)*m2)))-atan(sqrt(m2));

	return nu*(180.0/M_PI);
}

double ss_Mnu(double nu, double gamma)
{
	gg = gamma;
	gnu = nu;

	return brentq(ss_hMnu, 1.0, 50.0);
}
#endif
