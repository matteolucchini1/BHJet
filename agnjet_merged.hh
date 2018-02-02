// =====================================================================================
// 
//       Filename:  agnjet.hh
// 
//    Description:  header of agnjet.cc
// 
//        Version:  1.0
//        Created:  09/24/2012 17:57:51
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Samia Drappeau (sd), drappeau.samia@gmail.com
//        Company:  
// 
// =====================================================================================
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <stdarg.h>

#define NEBIN	10000
#define kpc	(1e3*GSL_CONST_CGSM_PARSEC)
#define pi	M_PI
#define cee	GSL_CONST_CGSM_SPEED_OF_LIGHT
#define herg	GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hkev	(GSL_CONST_CGSM_PLANCKS_CONSTANT_H*6.2415e8)
#define pmgm	GSL_CONST_CGSM_MASS_PROTON
#define emgm	GSL_CONST_CGSM_MASS_ELECTRON
#define kboltz	GSL_CONST_CGSM_BOLTZMANN // Boltzmann constant in erg/K
#define kboltz_kev 8.617e-8 // Boltzmann constant in keV/K
#define msun	GSL_CONST_CGSM_SOLAR_MASS
#define mjy	1.e-26
#define charg	4.8e-10
#define re0	2.81794e-13
#define sigtom	GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define gconst	GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT
#define sbconst GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT
#define aconst	7.56e-15
#define emerg	(GSL_CONST_CGSM_MASS_ELECTRON*GSL_CONST_CGSM_SPEED_OF_LIGHT*GSL_CONST_CGSM_SPEED_OF_LIGHT)

#define precision 30
#define SIZEMAX 10000

#define NR_END 		1
#define FREE_ARG 	char*
#define e_kin_fact 	9.1094*1e-28*pow(cee,2)*6.24*1e8

/**
 * Structures used in the code
 *
 */
struct synemis_params{
        double freq;
        double bfield;
        gsl_spline *spline_syn;
        gsl_interp_accel *acc_syn;
        gsl_spline *spline_eldis1;
        gsl_interp_accel *acc_eldis1;
};

struct absfnc_params{
        double freq;
        double bfield;
        gsl_spline *spline_syn;
        gsl_interp_accel *acc_syn;
        gsl_spline *spline_derivs;
        gsl_interp_accel *acc_derivs;
};

struct comint_params{
        double eph;
        double ephmin;
        double ephmax;
        gsl_spline *spline_eldis1;
        gsl_interp_accel *acc_eldis1;
        gsl_spline *spline_com1;
        gsl_interp_accel *acc_com1;
};

struct comfnc_params{
        double game;
        double e1;
        gsl_spline *spline_com1;
        gsl_interp_accel *acc_com1;
};

struct bbdiskfnc_params{
        double rin;
        double rout;
        double z;
        double hbb;
        double gamv;
        double tin;
};

struct bbfnc3_params{
        double gamv;
        double tin;
        double nu;
};

struct bbfnc_params{
        double gamv;
        double z;
        double hbb;
        double tin;
        double rin;
        double nu;
};

struct bbearth_params{
        double tin;
        double rin;
        double frq;
        double inclin;
        double dist;
};


void xrbjet(double* ear, int ne, double *params, double *phot_spect, double *photeng);
void xrbinterp(double *ear, double *energ, double *phot, double *photar, int ne, int newne);
void equipartition(double velsw, double mxsw, int nenpsw, double eta, double equip, double pspec, double nprot0, double emin, double emax, double endnsmj, double &cnorm, double &ntot0, double &b_en);
void jetpars(double velsw, double zeta, double z0, double z, double zsh, double r_g, double b0, double n0, double g0, double r0, double h0, gsl_spline *spline, gsl_interp_accel *acc, double &gb, double &dgb, double &b, double &n, double &g, double &r);
void bljetpars(double mxsw, double velsw, double z, double zacc, double r_g, double eta, double n0, double g0, double r0, double h0, double endnsmj, double pspec, double cnorm, double emin, double emax, double ebreak, gsl_spline *spline, gsl_interp_accel *acc, double &gb, double &b, double &n, double &g, double &r, double temp);
void b_prof(double mxsw, double gfin, double eta, double n, double endsmj, double pspec, double cnorm, double emin, double emax, double ebreak, double &field, double g, double temp);
double bbfnc(double thet, void *p);
double bbfnc3(double thet, void *p);
double bbdisk(double lnu, double gamv, double tin, double rin, double rout, double z, double hbb);
double bbdiskfnc(double lnu, void *p);
void bbintegrals(double lnub, double lnut, double gamv, double tin, double rin, double rout, double z, double hbb, double &bbint);
double synemis(double e1e, void *p);
double absfnc(double e2, void *p);
void synintegrals(double freq, double bfield, double elenmn, double elenmx, gsl_spline *spline_syn, gsl_interp_accel *acc_syn, gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_derivs, gsl_interp_accel *acc_derivs, double &esum, double &asum);
void synint(double freq, double r, double angle, double dopfac, double h, double bfield, double dist, double &esum, double &asum, double &fluxa, double &fluxa2);
double comfnc(double ein, void *p);
double comint(double gam, void *p);
double comintegral(double blim, double ulim, double eph, double ephmin, double ephmax, gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_com1, gsl_interp_accel *acc_com1);
double bbearth(double lr, void *p);
void bbearthint(double blim, double ulim, double frq, double tin, double rin, double dist, double inclin, double & bbflx);

double k0_fnc(double x);
double k1_fnc(double x);
double k2_fnc(double x);
double k3_fnc(double x);

extern double enden (double T, double ul);
