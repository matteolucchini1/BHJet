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
//	 Rewritten by:  Chiara Ceccobello
//  Maintained by:  Matteo Lucchini, m.lucchini@uva.nl
//        Company:  API - UvA 
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
#define kpc		(1e3*GSL_CONST_CGSM_PARSEC)
#define pi		M_PI
#define cee		GSL_CONST_CGSM_SPEED_OF_LIGHT
#define herg	GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hkev	(GSL_CONST_CGSM_PLANCKS_CONSTANT_H*6.2415e8)
#define pmgm	GSL_CONST_CGSM_MASS_PROTON
#define emgm	GSL_CONST_CGSM_MASS_ELECTRON
#define kboltz	GSL_CONST_CGSM_BOLTZMANN
#define kboltz_kev 8.617e-8 // Boltzmann constant in keV/K
#define msun	GSL_CONST_CGSM_SOLAR_MASS
#define mjy		1.e-26
#define charg	4.8e-10
#define re0		2.81794e-13
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
 * Structures used by the GSL libraries
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

/**
 * Functions used in the code
 */
 
void bhjet(double* ear,int ne,double *params,double *photeng,double *phot_spect);

//bb.cc
void bbdisk_comp(int disksw,int bbsw,double z,double reff2,double hbb,double tin,double rin,double rout,
double gamv,double tbbeff,double &ucom,double &uphdil,double &uphdil2);

double bbfnc(double thet,void *p);

double bbfnc3(double thet,void *p);

double bbdisk(double lnu,double gamv,double tin,double rin,double rout,double z,double hbb);

double bbdiskfnc(double lnu,void *p);

void bbintegrals(double lnub,double lnut,double gamv,double tin,double rin,double rout,double z,double hbb,
double &bbint);

double bbearth(double lr,void *p);

void bbearthint(double blim,double ulim,double frq,double tin,double rin, double dist,double inclin,
double &bbflx);

//compton.cc
void MICompton(bool isSSC,int Niter,int ncom,int njet,int nz,int nzdum,int k,double r,double ntot,double vol,
double dist,double dopfac[],gsl_spline *spline_eldis1,gsl_interp_accel *acc_eldis1,gsl_spline *spline_com2,
gsl_interp_accel *acc_com2,double eltemp,double elenmn,double elenmx,double &ephmin,double &ephmax,
double ephxr[],double comspc[],double nucom[],double totcom[]);

double comfnc(double ein,void *p);

double comint(double gam,void *p);

double comintegral(double blim,double ulim,double eph,double ephmin,double ephmax,gsl_spline *spline_eldis1,
gsl_interp_accel *acc_eldis1,gsl_spline *spline_com1,gsl_interp_accel *acc_com1);

void ucompton(int nz,int njet,int nsyn,int k,double dist,double r,double nusyn[],double dopfac[],
double synabs[],double &ucom);

void seed_sync_and_disk_phtns(bool isVerbose,int disksw,int bbsw,double tin,double rin,double rout,
double gamv,int nsyn,double snumax, double z,double r,double reff2,double hbb,double tbbeff,
double nphot[],double nurad[],double nubb[],double energ[], double phodis[],double ephot[],
double &ephmax,double &ephmin);

//dynamics.cc
void jetgrid(int k,int infosw,int nz,double r_g,double velsw,double h0,double r,double zmin,double zcut,
double zmax,int &zone_zcut,double &delz,double &z,double zed[]);

void jetpars(double velsw,double zeta,double z0,double z,double zsh,double r_g,double b0,double n0,double g0,
double r0,double h0,gsl_spline *spline,gsl_interp_accel *acc,double &gb,double &b,double &n,double &g,
double &r);

void vel_prof(double z,double mbh,double zmin,double visco,double r0,double h0,gsl_spline *spline,
gsl_interp_accel *acc,double &vel,double &dvel);

void bljetpars(double mxsw,double velsw,double z,double zacc,double r_g,double eta,double n0,double g0,
double r0,double h0,double endnsmj,double pspec,double cnorm,double emin,double emax,double ebreak, 
gsl_spline *spline,gsl_interp_accel *acc,double &gb,double &b,double &n,double &g,double &r,double sigsh);

void b_prof(double mxsw,double gfin,double eta,double n,double endsmj,double pspec,double cnorm,double emin,
double emax,double ebreak,double &field, double g,double sigsh);

void zonepars(int njet,int k,int nz,double z,double r,double delz,double rvel,double reff,double hbb,
double tbb2,double inclin,double gamax0,double gamax,int &nw,double &area,double &vol,double &gamv2,
double &gamv,double &beta,double &theff,double &tbbeff,double &gshift,double dopfac[]);

//electrons.cc
void electrons(int mxsw,int nelec,double pspec,double gamfac,double endnsmj,double ntot0,double cnorm,
double eltemp,double &gamax0,double &endens,double rdlgen[],double rdedn[],double thmbase[],double plbase[],
double &bete);

void pl_and_th_comp(bool isShock,double thmfrac,int &nw,int nelec,double mxsw,double rdlgen[],
double thmbase[],double ebreak,double emin,double emax,double gshift,double gshock,double mjteff,
double ratio_ne,double ntot,double cnorm,double enorm,double pspec,double &pltrm,double etemp[],
double dtemp[],double plcomp[],double thmcomp[],double &bete);

void update_ebreak(int nz,int k,double z,double r,double h0,double brk,double bfield,double beta,
double &ebreak,double &gebreak);

void read_and_shift_eldists(int nelec,int &nw,double ntot,double ntot0,double gshift,double rdlgen[],
double rdedn[],double elen[],double ledens[],double etemp[],double dtemp[]);

void lost_ele(int &nw,int nelec,double etemp[],double dtemp[],double ntot,double einc,double elen[],
double eled[],double lelec[],double ledens[]);

//equipartition.cc
void equipartition(double velsw,double mxsw,int nenpsw,double eta,double equip,double pspec,double nprot0,
double emin,double emax,double endnsmj,double &cnorm,double &ntot0,double &b_en,double &b0);

//max_ene.cc
void maxene(double z,double r,double ucom,double fsc,double bfield,double beta,double delz,double brk,
double &emax,double &gemax,double &ebreak,double &gebreak,double &bemax);

//mj_aven.cc
double K2(double x);

double mj(double g,double T);

double mjed(double g,double T);

void integrate(double ll,double ul,double T,double *res,double *err);

void integrate2(double ll,double ul,double T,double *res,double *err);

double enden(double T,double ul);

void aven(double T,double ul,double res);

double f_int(double g,void *p);

double f_int2(double g,void *p);

double enden (double T,double ul);

//norm_ele.cc
void elnorm(int flagNorm,double mxsw,double ntot,double pspec,double plfrac,double eltemp,double gshift,
double emax,double ebreak,double &mjteff,double &emin,double &cnseed_sync_and_disk_phtnsorm,double &enorm,
double &betat);

void norm_at_shock(double thmfrac,double heat,double pspec,double cnorm,double enorm,double betat,
double ebreak,double emax,double emin,double einc,int nelec,double lelec[],double ledens[],double etemp[],
double dtemp[],double thmshock[],double plshock[],double elen[],double ntot,double eled[],double shden[],
double shelen[],double &bete);

//pair.cc
void pproduction(int ncom,double r,double ephxr[],double totcom[],double &n_pp,double &rate_pp);

void annihilation(double z,double zcut,double r,double ntot,double eltemp,double &rate_pa,double &n_pa);

void pa_rate (double gamma_plus,double gamma_minus,double &rate);

void pa_rate_thermal (double numden,double T_e,double &rate);

//plots.cc
void create_plotFiles(int plotsw,int infosw);

void write_plotFiles(int plotsw,int bbsw,int disksw,double tin,double rin,double rout,double dist,
double inclin,double bbf1,double tbb2,double nutot,double &fplot,double complot,double presyn,
double postsyn,double bbplot);

//spectral_compnts.cc
void spcomponents(int infosw,int plotsw,int ne,int njet,int nz,int nsyn,int ncom,double zed[],double zcut,
double zsh,double ephxr[],double nusyn[],double synabs[],double nucom[],double comspc[],double cflx_array[],
double nutot[],double complot[],double presyn[],double postsyn[],double fplot[]);

//sycnhrotron.cc
void synchrotron(bool &isBreaknjetnsyn,int nsyn,int njet,int nz,int k,double bfield,double delz,double inclin,
double dist,double r,double elenmn,double elenmx,gsl_spline *spline_syn,gsl_interp_accel *acc_syn,
gsl_spline *spline_eldis1,gsl_interp_accel *acc_eldis1,gsl_spline *spline_derivs,gsl_interp_accel *acc_derivs,
double nurad[],double dopfac[],double nphot[],double nusyn[],double synabs[]);

double synemis(double e1e,void *p);

double absfnc(double e2,void *p);

void synintegrals(double freq,double bfield,double elenmn,double elenmx,gsl_spline *spline_syn,
gsl_interp_accel *acc_syn,gsl_spline *spline_eldis1,gsl_interp_accel *acc_eldis1,gsl_spline *spline_derivs,
gsl_interp_accel *acc_derivs,double &esum,double &asum);

void synint(double freq,double r,double angle,double dopfac,double h,double bfield,double dist,double &esum,
double &asum,double &fluxa,double &fluxa2);

//thermal.cc
double bb(double nu,double T);

double star(double nu,double r,double T,double d,double cos_theta,double gamma,int sw);

double mcd_obs(double nu,double Rin,double Tin,double d,double incl_deg);

double mcdobs1 (double nu, double Rin,double Tin,double Rout,double Tout,double d,double incl_deg);

double fjet1(double logr,void *p);

double mcd_jet(double nu,double Rin,double Tin,double z,double gamma);

double fjet11(double logr, void *p);

double mcdjet1(double nu,double Rin,double Tin,double Rout,double Tout,double z,double gamma);

//utilities.cc
void disk_init(int infosw,double jetrat,double rin,double rout,double eddlum,double bbf2,double tbb2,
double &tin,int &disksw,double &thrlum,double &eddrat,double &bbf1,double &reff,double &reff2,double &hbb);

void jet_init(int zfrac,int sizegb,double mxsw,double velsw,double jetrat,double r_g,double r0,double hratio,
double zacc,double zmax,double beta_pl,double &equip,double &h0,double &zsh,double &zcut,double &gad4_3,
double &betas0,double &gam0,double &rvel0,double &vel0,double &zmin,double gbx[],double gby[],
double gbx_vel2[],double gby_vel2[]);

void particles_init(double mxsw,double velsw,double jetrat,double r0,double eltemp,double pspec,double gamfac,
double sigsh,double plfrac,double gad4_3,double vel0,double gam0,double &thmfrac,double &emin,double &gamin,
double &ebreak,double &emax,double &uplim,double &endncom,double &endnsmj,double &betat,double &nprot0,
double &eta,double &equip,double &sig0,double &beta_pl);

double k0_fnc(double x);

double k1_fnc(double x);

double k2_fnc(double x);

double k3_fnc(double x);

void ene_arrays(bool &isSSC,int ne,int nsyn,int ncom,int &zfrac,double mbh,double inclin,int &njet,
double &snumin,double &snumax,double &snuinc,double &cnumin,double &cnumax,double &cnuinc,double ear[],
double nutot[],double nubb[],double nurad[],double energ[],double ephot[],double ephxr[]);

double average_ene(int flaglog,int N,double ndensity[],double energy[],double eldens[]);

void bhinterp(double *ear, double *energ, double *phot, double *photar, int ne, int newne);


/*velprof_func.cc, not used but included if one wants to play with the Euler equation

void rhs0(const double y,double &dydt,const double t);

void rhs1(const double y,double &dydt,const double t);

void rhs2(const double y,double &dydt,const double t);

void rhs3(const double y,double &dydt,const double t);

void write_cout(const double &y,double t);

void VelProf(int velsw,double Gam,double gbx[],double gby[],int &sizegb);
*/
