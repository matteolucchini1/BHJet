/** Blackbody function: for a blackbody at temperature T (K),
  returns the blackbody flux at frequency nu (Hz).
  The spectrum peaks at h*nu_max=2.82*kB*T or nu_max=5.88e10*T. 
  At frequencies sufficiently below the peak (Rayleigh-Jeans part), 
  bb(nu,T)=(2*kB/c^2)*T*nu**2.

  Also we assume that there is very little emission at frequencies 
  above 50*nu_max and below 1e-3*nu_max.

  Unit: bb(nu,T) returns ergs/s/cm**2/Hz/str
*/

#include <cmath>

#define h       6.626068e-27      ///< Planck's const       (erg*s)
#define kB      1.380650e-16      ///< Boltzmann constant   (erg/K)
#define c1      1.474499e-47      ///< 2*h/c/c in CGS
#define c2      3.072361e-37      ///< 2*kB/c/c in CGS

struct fobs_params {
	double nu;
};

struct fobs1_params{
	double nu; 
	double t1; 
	double t2; 
};

struct fjet_params{
	double nu;
	double t1; 
	double z; 
	double beta; 
	double gamma;
};

struct fjet1_params{
	double nu; 
	double t1; 
	double t2; 
	double beta; 
	double z; 
	double gamma;
};

double bb(double nu,double T){
	double nu_max, x=0;
	nu_max=5.88e10*T;
	if (nu < 1e-3*nu_max || nu > 5e1*nu_max) x=1e-99;
	else
		if	(nu < nu_max/1e2) x=c2*nu*nu*T;
	else {
		x=c1*nu*nu*nu/(exp(x)-1.);
	}

	return x;
}
#undef h
#undef kB
#undef c1
#undef c2

/** 
 *Computes the observed flux from a spherical blackbody emitter, like the companion star,  which has a radius
 *'r' (cm) and temperature of 'T' (K).  The observer is located at a distance of 'd' (cm) away from the center
 *of the star. It can return either the flux observed by the observer (cgs units) or the photon number density 
 *dN/dE/dV at observer (#/cm^3/erg).  
 *
 *Description: The flux at the position of the observer is calculated using the result of eq. (12) of Chia 
 *             1976, Astrophysics and Space Science, 46, 239. Also see Rybicki and Lightman eq. (1.13).      
 *
 *Requirements: function bb() to compute blackbody fluxes.
 *
 *Input: nu (Hz);
 *       r = radius of star (cm);
 *	 	 T = BB temp on stellar surface (K);
 *       d = distance from center of star to the observer (cm);
 *	 	 cos_theta = cos(theta) such that pi-theta is the angle between the vectors  V and JS where V is the 
 *	   				velocity vector of the observer and JS is a vector from the observer to the BB emitting         
 *	        		region
 *		 gamma = Lorentz factor of the observer.
 *	 	 switch = 1 or 2 (See output below)
 *	
 *Output: If switch = 1 : Flux observed by the observer (mJy).
 *                  = 2 : Photon density dN/dE/dV at observer (#/cm^3/erg).	 .
*/

#include <math.h>
#define pi M_PI
#define hhc 1.31623211e-42                     ///< h*h*c in CGS units

double star(double nu,double r,double T,double d,double cos_theta,double gamma,int sw){
	double bb (double nu, double T);
	double beta = sqrt(gamma*gamma-1.)/gamma;
	double result,Teff;
	Teff=T/gamma/(1.+beta*cos_theta);
	result=bb(nu,Teff)*pi*r*r/d/d;
	if (sw == 1) return 1e26*result;    // Return mJy
	if (sw == 2) return result/nu/hhc;  // Return #/cm^3/erg
	else return 0.0;
}
#undef pi
#undef hhc
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/** Spectra of a Multi-color disk blackbody at nu (Hz) with Tin at Rin,
  at a distance 'd', and inclination 'incl_deg'. 
  
  Requirements: 
    (1) GNU Scientific Library (http://www.gnu.org/software/gsl/) 
    must be installed.
    (2) Routine bb() to calculate the blackbody fluxes.

  Compiling: 
   To create the object file mcd_obs.o for this routine:
     gcc -c mcd_obs.c
   For the final compilation step, link properly with GSL routines, e.g.
     -lm -lgsl -lgslcblas mcd_obs.o

   If GSL is not installed in the default location, then in the final step
     replace "-lm -lgsl -lgslcblas" by
     "-I/path/to/gsl/include -L/path/to/gsl/lib -lm -lgsl -lgslcblas".

  Description:
  Following SS73, T~R^(-3/4).
  We assume Tout=1000 K. 
  Normalization, eq (4) of Mitsuda et al. (1984) PASC, 36, 741 is
	norm=8.*pi*Rin*Rin*cos(incl)/3./d/d;
  However we don't want Tin inside the integral, so we have
	norm*=pow(Tin, 8./3.);
  Then the integral becomes
	\int_{Tout}^}{Tin} dT T^(-11/3)*bb(nu,T) 

  Units: 
    Input:  Rin (cm), Tin (K), d (cm), incl_deg (deg).
    Output: Observed flux (mJy).
*/

#include <math.h>
#include <gsl/gsl_integration.h>

#define pi      M_PI
/** GSL QAG Integration parameters */
#define WORKSZ 200
#define EPSABS 0.0
#define EPSREL 1e-2
#define KEY    1

double fobs1 (double T, void * p){
	struct fobs_params * params = (struct fobs_params *)p;
	double nu = (params->nu);
	return pow (T,-11./3.) * bb(nu,T);
}

double mcd_obs(double nu,double Rin,double Tin,double d,double incl_deg){
	double incl=incl_deg*pi/180.;  // incl in rad
	double Tout=1000.;
	double norm, result, error;
	double bb (double nu, double T);

	struct fobs_params params;
	
	params.nu = nu;
	gsl_integration_workspace * w = 
	gsl_integration_workspace_alloc (WORKSZ);
	gsl_function F1;
	F1.function = &fobs1;
	F1.params   = &params;

	// Do the integration from Tout to Tin
	gsl_integration_qag(&F1,Tout,Tin,EPSABS,EPSREL,WORKSZ,KEY,w,&result,&error);
	gsl_integration_workspace_free(w);
	// Return observed flux in mJy
	norm=1e26*pow(Tin, 8./3.)*8.*pi*Rin*Rin*cos(incl)/3./d/d;
	result*=norm;
	return result;
}

// Includes irradiation + viscous disk
double fobs11(double logr,void *p){
	struct fobs1_params * params = (struct fobs1_params *)p;
	double nu    = (params->nu);
	double t1	= (params->t1);
	double t2	= (params->t2);
	double r     = exp(logr);
	double Teff;
	Teff = pow (t1*pow(r,-0.75), 4.) + 
	pow (t2*pow(r,-3./7.), 4.);
	Teff = pow(Teff,0.25);
	return bb(nu,Teff)*r*r;
}

double mcdobs1(double nu,double Rin,double Tin,double Rout,double Tout,double d,double incl_deg){
	double incl=incl_deg*pi/180.;  // incl in rad
	double logrin=log(Rin);
	double logrout=log(Rout);
	double t1 = Tin/pow(Rin,-0.75);      ///< so that Tvis(r)=t1*r^(-3/4)
	double t2 = Tout/pow(Rout,-3./7.);   ///< so that Tirr(r)=t2*r^(-3/7)
	double norm, result, error;
	double bb (double nu, double T);


	struct fobs1_params params;

	params.nu = nu;
	params.t1 = t1;
	params.t2 = t2;
	gsl_integration_workspace * w = 
		gsl_integration_workspace_alloc (WORKSZ);
	gsl_function F1;
	F1.function = &fobs11;
	F1.params   = &params;

        // Do the integration from Rin to Rout
	gsl_integration_qag  (&F1, logrin, logrout, EPSABS, EPSREL, WORKSZ, 
			KEY, w, &result, &error);
	gsl_integration_workspace_free(w);
	// Return observed flux in mJy
	norm=2.*pi*cos(incl)/d/d;
	result*=1e26*norm;
	return result;
}
#undef pi
#undef WORKSZ
#undef EPSABS
#undef EPSREL
#undef KEY
/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
/** Spectra of a multi-color disk blackbody as seen by a jet component.

  Requirements: 
    (1) GNU Scientific Library (http://www.gnu.org/software/gsl/) 
    must be installed.
    (2) Routine bb() to calculate the blackbody fluxes.

  Compiling: 
   To create the object file mcd_jet.o for this routine:
     gcc -c mcd_jet.c
   For the final compilation step, link properly with GSL routines, e.g.
     -lm -lgsl -lgslcblas mcd_jet.o

   If GSL is not installed in the default location, then in the final step
     replace "-lm -lgsl -lgslcblas" by
     "-I/path/to/gsl/include -L/path/to/gsl/lib -lm -lgsl -lgslcblas".

  Description:
  The jet is assumed to be along z axis, and is perpendicular to the 
  disk (on x-y plane). The Lorentz factor of the jet at 'z' is gamma.
  Rout of the disk is assumed to such that Tout(Rout)=1000 K.

  We do this integral over 'r', i.e. the total flux at z, at frequency nu is:

  f(z;nu)= int_{Rin}^{Rout} bb(nu,Tobs(r,z))*cos(theta)^2*(2*pi*r*dr)/d^2
         = int_{Rin}^{Rout} f(r) dr  (say)

  where:
    f(r) = bb(nu,Tobs(r,z))*cos(theta)*(2*pi*r)/d^2 ,
  and
    d^2=z^2 + r^2 ,
    theta=arctan(r/z) .
    Tobs(r,z)=(Tin/Rin^(-3/4)) * r^(-3/4) / (gamma*(1+beta*cos(theta)))

    The first couple of terms in the definition of Tobs(r,z) ensure a
    simple incarnation of the Shakura-Sunyaev (1973) disk, and the last 
    relativistic factor accounts for the Doppler shift in the diskBB
    radiation as observed by the moving jet.

  A simple integration scheme fails because typically the integration limits
  are 3 or so orders of magnitude different and it is important to integrate
  over the entire range with equal importance everywhere. Simple schemes fail
  to divide the intervals logarithmically. So we instead integrate over 
  log(r). Let z=log(r) => dz = dr/r or dr = r*dz = exp(z)dz. Then our 
  integration becomes:

  f(z;nu) = int_{logRin}^{logRout} f(exp(z)) * r dz

  and this function when supplied to regular integration routines, gives
  equal importance over the entire range of integration.

  The way the integral in the routine below is set up,
  2*pi*(result of integration) is the flux observed at 'z'.
  So the number density of seed disk photons at 'z'
  is dn(z)/d(nu) = f(nu;z)/c/h/nu
  However the Comptonization routine asks for dn/dE whereas this
  number is dn/d(nu). So we need to divide by another factor of h
  dn(z)/dE = f(nu;z)/c/h/nu/h to get the number density to supply to
  Comptonization routine.

  Units: 
    Input: Frequency (Hz), Rin (cm), Tin (K), z of jet (cm), 
           Lorentz factor of jet.
    Output: Observed photon density (#/cm^3/erg).
*/

#include <math.h>
#include <gsl/gsl_integration.h>

#define pi M_PI
#define hhc 1.31623211e-42                     ///< h*h*c in CGS units
/** GSL QAG Integration parameters */
#define WORKSZ 200
#define EPSABS 0.0
#define EPSREL 1e-2
#define KEY    1

double fjet1 (double logr, void * p){
	struct fjet_params * params = (struct fjet_params *)p;
	double nu    = (params->nu);
	double t1	= (params->t1);
	double z	= (params->z);
	double beta	= (params->beta);
	double gamma	= (params->gamma);
	double r     = exp(logr);
	double dd,cos_theta,cos2_theta,Teff,sin2_theta;
	dd=r*r+z*z;
	cos2_theta = z*z/dd;
	sin2_theta = r*r/dd;
	cos_theta = sqrt(cos2_theta);
	// Temp of BB observed by jet at z
	Teff=t1*pow(r,-0.75)/gamma/(1.+beta*cos_theta);
	return bb(nu,Teff)*cos2_theta*sin2_theta;
}

double mcd_jet(double nu,double Rin,double Tin,double z,double gamma){
	double Tout = 1e3;
	double Rout = Rin*pow(Tout/Tin, -4./3.);
	double logrin=log(Rin);
	double logrout=log(Rout);
	double beta = sqrt(gamma*gamma-1.)/gamma;
	double t1   = Tin/pow(Rin,-0.75);      ///< so that T(r)=t1*r^(-3/4)
	double result, error;
	double bb (double nu, double T);


	struct fjet_params params;
	params.nu=nu;
	params.t1=t1;
	params.z= z;
	params.beta=beta;
	params.gamma=gamma;

	gsl_integration_workspace * w = 
	gsl_integration_workspace_alloc (WORKSZ);
    gsl_function F1;
    F1.function = &fjet1;
    F1.params = &params;
    // Do the integration from Rin to Rout
    gsl_integration_qag(&F1,logrin,logrout,EPSABS,EPSREL,WORKSZ,KEY,w,&result,&error);
    gsl_integration_workspace_free(w);
    result*=2.*pi/nu/hhc;   // Photon density
    
    return result;
}

// Includes irradiation + viscous disk
double fjet11 (double logr, void * p){
	struct fjet1_params * params = (struct fjet1_params *)p;
	double nu	= (params->nu);
	double t1	= (params->t1);
	double t2	= (params->t2);
	double beta	= (params->beta);
	double z	= (params->z);
	double gamma	= (params->gamma);
	double r     = exp(logr);
	double dd,cos_theta,cos2_theta,Teff,sin2_theta;
	dd=r*r+z*z;
	cos2_theta = z*z/dd;
	sin2_theta = r*r/dd;
	cos_theta = sqrt(cos2_theta);

	// Effective temperature due to irr+visc
	Teff = pow (t1*pow(r,-0.75), 4.) + 
	pow (t2*pow(r,-3./7.), 4.);
	Teff = pow(Teff,0.25);

	// Temp of BB observed by jet at z
	Teff=Teff/gamma/(1.+beta*cos_theta);
	return bb(nu,Teff)*cos2_theta*sin2_theta;
}

double mcdjet1(double nu,double Rin,double Tin,double Rout,double Tout,double z,double gamma){
	double logrin=log(Rin);
	double logrout=log(Rout);
	double beta = sqrt(gamma*gamma-1.)/gamma;
	double t1 = Tin/pow(Rin,-0.75);      ///< so that Tvis(r)=t1*r^(-3/4)
	double t2 = Tout/pow(Rout,-3./7.);   ///< so that Tirr(r)=t2*r^(-3/7)
	double result, error;
	double bb (double nu, double T);

	struct fjet1_params params;
	params.nu=nu;
	params.t1=t1;
	params.t2=t2;
	params.beta=beta;
	params.z=z;
	params.gamma=gamma;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (WORKSZ);
    gsl_function F1;
    F1.function = &fjet11;
    F1.params = &params;
    // Do the integration from Rin to Rout
    gsl_integration_qag(&F1,logrin,logrout,EPSABS,EPSREL,WORKSZ,KEY,w,&result,&error);
    gsl_integration_workspace_free(w);
    result*=2.*pi/nu/hhc;   // Photon density
    //result*=2.*pi*1e26;   // Flux in mJy
    return result;
}

#undef pi
#undef hhc
#undef WORKSZ
#undef EPSABS
#undef EPSREL
#undef KEY
/////////////////////////////////////////////////////////////////////////

