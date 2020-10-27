#include "Radiation.hpp"
#include "Cyclosyn.hpp"

#include <iostream>

//Synchrotron tables for F(nu/nuc) for calculation of single particle spectrum
static double arg[47]	= {0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.03, 0.05, 0.07, 0.1,
0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5,
8., 8.5, 9., 9.5, 10., 12., 14., 16., 18., 20., 25., 30., 40., 50.};
    
static double var[47]	= {-1.002e+0, -9.031e-1, -7.696e-1, -6.716e-1, -5.686e-1, -4.461e-1, -3.516e-1,
-2.125e-1, -1.537e-1, -1.192e-1, -8.725e-2, -4.383e-2, -3.716e-2, -4.528e-2, -5.948e-2, -7.988e-2,
-1.035e-1, -1.296e-1, -1.586e-1, -1.838e-1, -3.507e-1, -5.214e-1, -6.990e-1, -8.861e-1, -1.073e+0,
-1.267e+0, -1.470e+0, -1.670e+0, -1.870e+0, -2.073e+0, -2.279e+0, -2.483e+0, -2.686e+0, -2.893e+0,
-3.097e+0, -3.303e+0, -3.510e+0, -3.717e+0, -4.550e+0, -5.388e+0, -6.230e+0, -7.075e+0, -7.921e+0,
-1.005e+1, -1.218e+1, -1.646e+1, -2.076e+1};

Cyclosyn::~Cyclosyn(){
	delete[] en_phot;
	delete[] num_phot;
	delete[] en_phot_obs;
	delete[] num_phot_obs;
		
	gsl_integration_workspace_free (w1);
	
	gsl_spline_free(syn_f), gsl_interp_accel_free(syn_acc);
}

//This constructor initializes the arrays and interpolations. In this case, calculations are done in frequency
//space, not in photon energies.
Cyclosyn::Cyclosyn(int s,double numin,double numax,double m){
	size = s;

	en_phot = new double[size];
	num_phot = new double[size];
	en_phot_obs = new double[2*size];
	num_phot_obs = new double[2*size];
	
	counterjet = false;
	
	w1 = gsl_integration_workspace_alloc(100);
	
	syn_acc =  gsl_interp_accel_alloc();
	syn_f = gsl_spline_alloc(gsl_interp_cspline,47);
    
    gsl_spline_init(syn_f,arg,var,47); 
	
	double nuinc = (log10(numax)-log10(numin))/(size-1);
	
	for(int i=0;i<size;i++){
		en_phot[i] = pow(10.,log10(numin)+i*nuinc)*herg;
		num_phot[i] = 0;
	}	
	
	mass = m;
}

//Single particle emissivity/absorption coefficient calculations
double cyclosyn_emis(double gamma,void *p){
	struct cyclosyn_emis_params *params = (struct cyclosyn_emis_params *)p;
    double nu = (params -> nu);
    double b = (params -> b);
    double mass = (params -> mass);
    gsl_spline *syn = (params -> syn);
    gsl_interp_accel *acc_syn = (params -> acc_syn);
    gsl_spline *eldis = (params -> eldis);
    gsl_interp_accel *acc_eldis = (params -> acc_eldis);
    
  	double nu_c, x, emisfunc, nu_larmor, psquared, ngamma;

	gamma = exp(gamma); 
    
    //this is in the synchrotron regime
 	if (gamma > 2.) {
		nu_c = (3.*charg*b*pow(gamma,2.))/(4.*pi*mass*cee);
		x = nu/nu_c;
	    if(x <= 1.e-4){
        	emisfunc = 4.*pi*pow(x/2.,(1./3.))/(sqrt(3.)*2.68);
    	}
    	else if(x > 50.){
        	emisfunc = sqrt(pi*x/2.)*exp(-x);
    	}
    	else{
        	emisfunc = pow(10.,gsl_spline_eval(syn, x, acc_syn));
    	}
    } else { //cyclotron regime
 		nu_larmor = (charg*b)/(2.*pi*mass*cee);
 		x = nu/nu_larmor;
 		psquared = pow(gamma,2.)-1.;
 		emisfunc = (2.*psquared)/(1.+3.*psquared)*exp((2.*(1.-x))/(1.+3.*psquared));
 	}
    ngamma = gsl_spline_eval(eldis,gamma,acc_eldis);

 	return ngamma*gamma*emisfunc;
}

double cyclosyn_abs(double gamma,void *p){
	struct cyclosyn_abs_params *params = (struct cyclosyn_abs_params *)p;
	double nu = (params -> nu);
    double b = (params -> b);
    double mass = (params -> mass);
    gsl_spline *syn = (params -> syn);
    gsl_interp_accel *acc_syn = (params -> acc_syn);
    gsl_spline *derivs = (params -> derivs);
    gsl_interp_accel *acc_derivs = (params -> acc_derivs);
    
    double nu_c, x, emisfunc, nu_larmor, psquared, ngamma_diff;

	gamma = exp(gamma);
	
	//this is in the synchrotron regime
 	if (gamma > 2.) {
		nu_c = (3.*charg*b*pow(gamma,2.))/(4.*pi*mass*cee);
		x = nu/nu_c;
	    if(x <= 1.e-4){
        	emisfunc = 4.*pi*pow(x/2.,(1./3.))/(sqrt(3.)*2.68);
    	}
    	else if(x > 50.){
        	emisfunc = sqrt(pi*x/2.)*exp(-x);
    	}
    	else{
        	emisfunc = pow(10.,gsl_spline_eval(syn, x, acc_syn));
    	}
    } else { //cyclotron regime
 		nu_larmor = (charg*b)/(2.*pi*mass*cee);
 		x = nu/nu_larmor;
 		psquared = pow(gamma,2.)-1.;
 		emisfunc = (2.*psquared)/(1.+3.*psquared)*exp((2.*(1.-x))/(1.+3.*psquared));
 	}
 	
 	ngamma_diff = gsl_spline_eval(derivs,gamma,acc_derivs);

	return ngamma_diff*pow(gamma,2.)*emisfunc;
}

//Integrals of single particle emissivity/absorption coefficient over particle distribution 
double Cyclosyn::emis_integral(double nu,double gmin,double gmax,gsl_spline *eldis,gsl_interp_accel
	 						   *acc_eldis){					 	    
    double result1, error1;
    gsl_function F1;
    struct cyclosyn_emis_params F1params = {nu,bfield,mass,syn_f,syn_acc,eldis,acc_eldis};
    F1.function     = &cyclosyn_emis;
    F1.params       = &F1params;
    gsl_integration_qag(&F1,log(gmin),log(gmax),1e1,1e1,100,2,w1,&result1,&error1);

    return result1;
}

double Cyclosyn::abs_integral(double nu,double gmin,double gmax,gsl_spline *derivs,gsl_interp_accel 
							  *acc_derivs){							   
    double result1, error1;
    gsl_function F1;
    struct cyclosyn_abs_params F1params = {nu,bfield,mass,syn_f,syn_acc,derivs,acc_derivs};
    F1.function     = &cyclosyn_abs;
    F1.params       = &F1params;
    gsl_integration_qag(&F1,log(gmin),log(gmax),1e1,1e1,100,2,w1,&result1,&error1);

    return result1;						 
}

//Comoving and observed specific luminosity for frequency bin k, from integrated coefficients
void Cyclosyn::specific_luminosity(int k,double gmin,double gmax, gsl_spline *eldis,gsl_interp_accel 
								   *acc_eldis,gsl_spline *eldis_diff,gsl_interp_accel *acc_eldis_diff){
	double emis, abs;
	double pitch = 0.73;
    double acons, elcons, asyn, epsasyn;
    double absfac, tsyn, tsyn_obs, absfac_obs;
    double dopfac_cj;
    
    dopfac_cj = dopfac*(1.-beta*cos(angle))/(1.+beta*cos(angle));
    en_phot_obs[k] = en_phot[k]*dopfac;
    if(counterjet == true){ 
    	en_phot_obs[k+size] = en_phot[k]*dopfac_cj;
    }
    
    emis = emis_integral(en_phot[k]/herg,gmin,gmax,eldis,acc_eldis);
    abs = abs_integral(en_phot[k]/herg,gmin,gmax,eldis_diff,acc_eldis_diff);
    
   	if (log10(emis) < -50. || log10(abs) < -50.){
    	num_phot_obs[k] = 0;
    	if(counterjet == true){ 
    		num_phot_obs[k+size] = 0;
    	}
    	return;
    }
    
    elcons = sqrt(3.)*(charg*charg*charg)*bfield*sin(pitch)/emerg;
    acons = -cee*cee/(8.*pi*pow(en_phot[k]/herg,2.));
    asyn = acons*elcons*abs;
    epsasyn	= emis/(acons*abs);
    
    tsyn	= pi/2. * asyn*r;
    if(tsyn >= 1.){
        absfac	= (1.-exp(-tsyn));
    }
    else{
        absfac	= tsyn-pow(tsyn,2.)/2.+pow(tsyn,3.)/6.;
    }    
    
    //This includes skin depth/viewing angle effects; the calculation differs from above by ~a few at most
    tsyn_obs = pi/2.*asyn*r/(dopfac*sin(angle));    
    if(tsyn_obs >= 1.){
        absfac_obs	= (1.-exp(-tsyn_obs));
    }
    else{
        absfac_obs	= tsyn_obs-pow(tsyn_obs,2.)/2.+pow(tsyn_obs,3.)/6.;
    }    
    
    num_phot[k]	= pi*r*r*absfac*epsasyn;    
    num_phot_obs[k]	= 2.*r*z*absfac_obs*epsasyn*pow(dopfac,dopnum);
    
    if(counterjet == true){    	
   		tsyn_obs = pi/2.*asyn*r/(dopfac_cj*sin(angle));    
   		if(tsyn_obs >= 1.){
        	absfac_obs	= (1.-exp(-tsyn_obs));
    	}
    	else{
        	absfac_obs	= tsyn_obs-pow(tsyn_obs,2.)/2.+pow(tsyn_obs,3.)/6.;
    	}
    	num_phot_obs[k+size] = 2.*r*z*absfac_obs*epsasyn*pow(dopfac_cj,dopnum);
    } else {
    	num_phot_obs[k+size] = 0;
    }
}

//Cyclosyonchrotron spectrum calculated looping over the entire frequency array
void Cyclosyn::cycsyn_spectrum(double gmin,double gmax,gsl_spline *eldis,gsl_interp_accel *acc_eldis,
							   gsl_spline *eldis_diff,gsl_interp_accel *acc_eldis_diff){			
    for(int k=0;k<size;k++){
    	specific_luminosity(k,gmin,gmax,eldis,acc_eldis,eldis_diff,acc_eldis_diff);
    }
}

//Methods to return the cyclosyn scale frequencies for a given Lorentz factor. Returns synchrotron or larmor
//frequencies depending on whether input gamma is greater or smaller than 2
//The second method does the same thing, but it checks the computed arrays and looking for the maximum value
//of L_nu. Unlike the simple way, this accounts for the fact that the peak may be caused by synchrotron self
//absorption rather than coinciding with the scale frequency
double Cyclosyn::nu_syn(double gamma){
		return (3.*charg*bfield*pow(gamma,2.))/(4.*pi*mass*cee);
}

double Cyclosyn::nu_syn(){
	double temp_lum = 0.;
	int temp = 0;
	for (int i=0;i<size;i++){
		if (num_phot[i] > temp_lum){
			temp_lum = num_phot[i];
			temp = i;
		}
	}
	return en_phot[temp]/herg;
}


//Method to set magnetic field
void Cyclosyn::set_bfield(double b){
	bfield = b;
}

void Cyclosyn::test(){
	std::cout<< "Bfield: " << bfield << " r: " << r << " z: " << z << " v.angle: " << angle << " speed: " 
	<< beta	<< " delta: " << dopfac << std::endl;
}

