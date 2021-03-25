#include <stdarg.h>
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
#include <omp.h>

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
#include "../Kariba/Bknpower.hpp"
#include "../Kariba/Mixed.hpp"
#include "../Kariba/Kappa.hpp"

#include "../Kariba/Radiation.hpp"
#include "../Kariba/ShSDisk.hpp"
#include "../Kariba/Cyclosyn.hpp"
#include "../Kariba/Compton.hpp"
#include "../Kariba/BBody.hpp"

#define kpc	(1e3*GSL_CONST_CGSM_PARSEC)
#define cee				GSL_CONST_CGSM_SPEED_OF_LIGHT
#define emgm			GSL_CONST_CGSM_MASS_ELECTRON
#define pmgm			GSL_CONST_CGSM_MASS_PROTON
#define kboltz	GSL_CONST_CGSM_BOLTZMANN
#define kboltz_kev2erg	1.6022e-9 		//Boltzman constant in keV/erg 
#define me_kev			511 			//electron mass in kev
#define emerg			(GSL_CONST_CGSM_MASS_ELECTRON*pow(GSL_CONST_CGSM_SPEED_OF_LIGHT,2.))
#define pi		M_PI
#define charg	4.8e-10
#define sigtom	GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define herg	GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hkev	(GSL_CONST_CGSM_PLANCKS_CONSTANT_H*6.2415e8)
#define mjy		1.e-26
#define re0		2.81794e-13
#define gconst	GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT
#define sbconst GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT
#define aconst	7.56e-15
#define msun	GSL_CONST_CGSM_SOLAR_MASS

using namespace std;

void plot_write(int size,double *en,double *lum,char path[],double dist,double redshift);					
void plot_write(int size,const double *en,const double *lum,char path[],double dist,double redshift);	
void plot_write(int size,const double *p,const double *g,const double *pdens,const double *gdens,
				char path[]);

void sum_zones(int size_in,int size_out,double* input_en,double* input_lum,double* en,double* lum);
void sum_ext(int size_in,int size_out,const double* input_en,const double* input_lum,double* en,double* lum);
double integrate_lum(int size,double numin,double numax,const double* input_en,const double* input_lum);
double photon_index(int size,double numin,double numax,const double* input_en,const double* input_lum);

void clean_file(char path[],bool check);
void read_params(string file,double *pars);





