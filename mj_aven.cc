/* Purpose: Find average energy for Maxwell Juttner distribution

	Version August 5, 2010 (SM)

 Compile: 
   gcc -m32 -lm -lgsl -lgslcblas mj_aven.c -o mj_aven

*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

/* Relativistic Maxwellian or Maxwell-Juttner distribution */ /*{{{*/
/* for electrons/positrons.
 *   mj(g,T) = x*g*g*b*exp(-g*x)/K2(x) , where
 *     g  = gamma, the Lorentz factor
 *     b  = beta, g=1/sqrt(1-b*b)
 *     x  = m*c*c/k/T = 5.92988981e9/T
 *     K2 = Modified  Bessel function of the second kind
 *     
 * E.g. Svensson (1982), ApJ, 258, 321, eq.64.
 * */
 
struct f_params {
	double T;
}; 
 
double K2(double x){

	double res;
	if (x < 0.1){
		res = 2./x/x;
	}
	else {
		res = gsl_sf_bessel_Kn(2,x);
	}
	return res;
}
 
double mj(double g,double T){

	double b,x,res;

	if (g!=1.){
		b = sqrt(1.-1./g/g);
		x = 5.92988981e9/T;
		res=x*g*g*b*exp(-g*x)/K2(x);
	}
	else res=0.;

	return res;
}

double mjed(double g,double T){

	double b,x,res;
	
	if (g!=1.){
		b = sqrt(1.-1./g/g);
		x = 5.92988981e9/T;
		res=x*g*g*g*b*exp(-g*x)/K2(x);
	}
	else res=0.;
	
	return res;
}

double f_int(double g,void * p) {

	struct f_params * params = (struct f_params *)p;
	double T=(params->T);
	return mj (g, T) ;
}


double f_int2(double g,void * p) {

	struct f_params * params = (struct f_params *)p;
	double T=(params->T);
	return mjed (g, T) ;
}


void integrate(double ll, double ul,double T,double *res,double *err){

	double result,error;
	size_t neval;
	struct f_params params;
	gsl_function F;

	params.T=T;
	F.function = &f_int;
	F.params = &params;

	gsl_integration_qng (&F,ll,ul,1e-6,1e-3,&result,&error,&neval);
	*res = result;
	*err = error;

	return;
}

void integrate2(double ll,double ul,double T,double *res,double *err){

	double result,error;
	size_t neval;
	struct f_params params;
	gsl_function F;
	
	params.T=T;
	F.function = &f_int2;
	F.params = &params;
	
	gsl_integration_qng (&F,ll,ul,1e-6,1e-3,&result,&error,&neval);
	*res = result;
	*err = error;
	
	return;
}

void aven(double T,double ul,double res){ 

  	double ll, area, error, area2, error2;	

	ll=1. + 1.e-6;
	integrate(ll,ul,T,&area,&error);
	integrate2(ll,ul,T,&area2,&error2);

	res = area2*5.11e2/(area*8.617e-8*T);
} 

double enden(double T, double ul){

  	double ll, area2, error2;
  	double res;

	ll=1. + 1.e-6;
	integrate2(ll,ul,T,&area2,&error2);	
	res = area2;
	
	return res;
} 
