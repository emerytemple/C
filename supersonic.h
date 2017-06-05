#ifndef __SUPERSONIC_H__
#define __SUPERSONIC_H__

double Supersonic(char *string, ...);

#if 0
// 1-D
double ss_TsTo(double gamma); // T*/To
double ss_PsPo(double gamma); // P*/Po
double ss_RsRo(double gamma); // R*/Ro

double ss_ToT(double Mach, double gamma); // To/T
double ss_PoP(double Mach, double gamma); // Po/P
double ss_RoR(double Mach, double gamma); // Ro/R

double ss_Mpr(double PR, double gamma); // find Mach given Po/P

// 1-D w/ heat addition
double ss_PPsh(double Mach, double gamma); // P/P*
double ss_TTsh(double Mach, double gamma); // T/T*
double ss_RRsh(double Mach, double gamma); // R/R*

double ss_PoPosh(double Mach, double gamma); // Po/Po*
double ss_ToTosh(double Mach, double gamma); // To/To*

// 1-D w/ friction
double ss_TTsf(double Mach, double gamma); // T/T*
double ss_PPsf(double Mach, double gamma); // P/P*
double ss_RRsf(double Mach, double gamma); // R/R*

double ss_PoPosf(double Mach, double gamma); // Po/Po*

double ss_flD(double Mach, double gamma); // (4FL*)/D

// quasi 1-D
double ss_AAs(double Mach, double gamma); // A/A*
double ss_Msub(double xp, double gamma); // find M < 1 given A/A*
double ss_Msup(double xp, double gamma); // find M > 1 given A/A*

// normal shock (To1 = To2)
double ss_M2(double Mach, double gamma); // M2
double ss_P2P1(double Mach, double gamma); // P2/P1
double ss_R2R1(double Mach, double gamma); // R2/R1
double ss_T2T1(double Mach, double gamma); // T2/T1
double ss_Po2Po1(double Mach, double gamma, double R); // Po2/Po1
double ss_Po2P1(double Mach, double gamma, double R); // Po2/P1

// oblique shock
double ss_theta(double Mach, double gamma, double beta);
double ss_betas(double Mach, double gamma, double theta);
double ss_betaw(double Mach, double gamma, double theta);

double ss_Mob(double theta, double beta, double gamma);

// prandtl-meyer expansion waves
double ss_mu(double Mach);
double ss_nu(double Mach, double gamma);

double ss_Mnu(double nu, double gamma);
#endif
#endif
