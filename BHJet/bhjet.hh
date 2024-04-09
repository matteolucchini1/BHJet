#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <fenv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

#include "../Kariba/Particles.hpp"
#include "../Kariba/Thermal.hpp"
#include "../Kariba/Powerlaw.hpp"
#include "../Kariba/Mixed.hpp"
#include "../Kariba/Kappa.hpp"
#include "../Kariba/Bknpower.hpp"

#include "../Kariba/Radiation.hpp"
#include "../Kariba/ShSDisk.hpp"
#include "../Kariba/Cyclosyn.hpp"
#include "../Kariba/Compton.hpp"
#include "../Kariba/BBody.hpp"


#define kpc             (1e3*GSL_CONST_CGSM_PARSEC)
#define cee             GSL_CONST_CGSM_SPEED_OF_LIGHT
#define emgm            GSL_CONST_CGSM_MASS_ELECTRON
#define pmgm            GSL_CONST_CGSM_MASS_PROTON
#define kboltz          GSL_CONST_CGSM_BOLTZMANN
#define kboltz_kev2erg  1.6022e-9 		//Boltzman constant in keV/erg 
#define emerg           (GSL_CONST_CGSM_MASS_ELECTRON*pow(GSL_CONST_CGSM_SPEED_OF_LIGHT,2.))
#define pi              M_PI
#define charg           4.8e-10
#define sigtom          GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define herg            GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hkev            (GSL_CONST_CGSM_PLANCKS_CONSTANT_H*6.2415e8)
#define mjy             1.e-26
#define re0             2.81794e-13
#define gconst          GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT
#define sbconst         GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT
#define aconst          7.56e-15
#define msun            GSL_CONST_CGSM_SOLAR_MASS

//Most functions in the code use input parameters arranged in a structure rather than passed as a long list of
//multiple int/double variables. The reason for this is imply to make the code easier to read and understand.
//Performance is not impacted.

//Structure including dynamical jet parameters 
typedef struct jet_dynpars{
    double min;				//jet launching point
    double max;				//max distance for jet calculations
    double h0;				//jet nozzle/corona height
    double r0;				//jet initial radius
    double acc;				//jet magnetic acceleration end location
    double beta0;			//jet initial speed in units of c
    double gam0;			//jet initial Lorentz factor
    double gamf;			//jet final Lorentz factor (only used with magnetic acceleration)		
    double Rg;				//gravitational radius
} jet_dynpars;

//Structure including parameters of jet energetics
typedef struct jet_enpars{
    double av_gamma;		//average Lorentz factor of electrons 
    double pbeta;			//plasma beta (Ue/Ub)
    double Nj;				//injected jet power
    double bfield;			//magnetic field strength 
    double lepdens;			//lepton number density
    double protdens;		//proton number density	
    double eta;				//pair content of the jet, ne/np
    double sig0;			//initial magnetization (Ub+Pb)/Up
    double sig_acc;			//final magnetization; values of sigma only used for magnetic acceleration
} jet_enpars;

//Structure including distance grid calculations
typedef struct grid_pars{
    int nz;					//total number of zones
    int cut;				//zone counter after which grid switched to logarithmic spacing
    double zcut;			//distance at which grid switched to log spacing
} grid_pars;

//Structure with parameters of each zone, used to check whether to calculate IC
typedef struct zone_pars{
    double gamma;			//lorentz factor of zone	
    double beta;			//speed of zone in units of c
    double delta;			//doppler factor of zone
    double r;				//radius of zone
    double delz;			//height of zone
    double bfield;			//magnetic field in zone
    double lepdens;			//number density of zone
    double avgammasq;		//avg lorentz factor squared in zone
    double eltemp;			//particle temperature in zone
    double nth_frac;        //fraction of non-thermal particles in the zone
} zone_pars;

//Structure with parameters needed for external inverse Compton photon fields
typedef struct com_pars{
    double lblr;			//luminosity of BLR
    double ublr;			//energy density of BLR
    double tblr;			//temperature of BLR
    double rblr;			//radius of BLR
    double ldt;				//luminosity of torus
    double udt;				//energy density of torus
    double tdt;				//temperature of torus
    double rdt;				//radius of torus
    double urad_total;		//total energy density
} com_pars;


void jetmain(double *ear,int ne,double *param,double *photeng,double *photspec);

void param_write(const double *par,std::string path);
void plot_write(int size,double *en,double *lum,std::string path,double dist,double redshift);					
void plot_write(int size,const double *en,const double *lum,std::string path,double dist,double redshift);	
void plot_write(int size,const double *p,const double *g,const double *pdens,const double *gdens,
                std::string path);
									
bool Compton_check(bool IsShock,int i,double Mbh,double Nj,double Ucom,double velsw,zone_pars &zone);

void sum_counterjet(int size,const double* input_en,const double* input_lum,double* en,double* lum);
void output_spectrum(int size,double* en,double* lum,double* spec,double redsh,double dist);
void sum_zones(int size_in,int size_out,double* input_en,double* input_lum,double* en,double* lum);
void sum_ext(int size_in,int size_out,const double* input_en,const double* input_lum,double* en,double* lum);
double integrate_lum(int size,double numin,double numax,const double* input_en,const double* input_lum);
double photon_index(int size,double numin,double numax,const double* input_en,const double* input_lum);

void velprof_ad(gsl_spline *spline);
void velprof_iso(gsl_spline *spline);
void velprof_mag(jet_dynpars &dyn,gsl_spline *spline);

void equipartition(int npsw,jet_dynpars &dyn,jet_enpars &en);
void equipartition(double Nj,jet_dynpars &dyn,jet_enpars &en);

void jetgrid(int i,grid_pars &grid,jet_dynpars &dyn,double r,double &delz,double &z);
void isojetpars(double z,jet_dynpars &dyn,jet_enpars &en,double &t,zone_pars &zone,gsl_spline *spline,
                gsl_interp_accel *acc);
void adjetpars(double z,jet_dynpars &dyn,jet_enpars &en,double &t,zone_pars &zone,gsl_spline *spline,
               gsl_interp_accel *acc);
void bljetpars(double z,double brk,jet_dynpars &dyn,jet_enpars &en,double &t,zone_pars &zone,gsl_spline 
               *spline,gsl_interp_accel *acc);
void b_profile(double g,double n,jet_dynpars &dyn,jet_enpars &en,double &field);

void agn_photons_init(double lum,double f1,double f2,com_pars &agn_com);
void zone_agn_phfields(double z,zone_pars &zone,double &ublr_zone,double &udt_zone,com_pars &agn_com);

void clean_file(std::string path,int check);
void jetinterp(double *ear,double *energ,double *phot,double *photar,int ne,int newne);

void ebl_atten_gil(int size,double* en,double* lum,double redsh);
