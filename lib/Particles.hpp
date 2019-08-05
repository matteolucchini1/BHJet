#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

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

//Template class for particle distributions
//This class contains members and methods that are used for thermal, non-thermal and mixed distributions
//note: ndens is number density per unit momentum

//Structures used for GSL integration
typedef struct pl_params { 
	double s; 
	double n;
} pl_params;

typedef struct bkn_params { 
	double s1;
	double s2;
	double brk;
	double max; 
	double m;
} bkn_params;

typedef struct th_params { 
	double t; 
	double n;
	double m; 
} th_params;

typedef struct k_params { 
	double t;
	double k; 
} k_params;

typedef struct injection_mixed_params {
	double s; 
	double t; 
	double nth;
	double npl;
	double m; 
	double min;
	double max;
	double cutoff;
} injection_mixed_params;

typedef struct injection_kappa_params {
	double t; 
	double k;
	double n; 
	double m;
} injection_kappa_params;

typedef struct injection_pl_params {
	double s;
	double n; 
	double m; 
	double max;
} injection_pl_params;

typedef struct injection_bkn_params {
	double s1;
	double s2;
	double brk;
	double max; 
	double m;
	double n;
} injection_bkn_params;

class Particles {
	protected:
		double mass;
		int size, type;
		double *p;				//array of particle momenta
		double *ndens;			//array of number density per unit volume, per unit momentum
		double *gamma;			//array of particle kinetic energies for each momentum
		double *gdens;		//array of number density per unit volume, per unit gamma
		double *gdens_diff;//array with differential of number density for radiation calculation
								//to treat non-relativistic particles the cooling has to be calculated in
								//momentum space, but relativistic radiation codes require gamma/energy space
								//In this class both can be accessed after calculating the appropriate cooling
		bool FP_flag;			//This flag regulates whether the particle distribution is going to include a
								//FP treatement for or not and sets up arrays needed for radiation accordingly
		gsl_integration_workspace *w1;
	public:
		~Particles();
		
		void initialize_gdens();	
		void initialize_pdens();
		void gdens_differentiate();
		
		const double *get_p()const				{ return p; }
		const double *get_pdens()const			{ return ndens; }
		const double *get_gamma()const			{ return gamma; }
		const double *get_gdens()const			{ return gdens; }
		const double *get_gdens_diff()const		{ return gdens_diff; }
					
		double beta(int i);							
		
		double count_particles();
		double count_particles_energy();
		double av_p();
		double av_gamma();
		double av_psq();
		double av_gammasq();
		
		void array_test();				
};
#endif
