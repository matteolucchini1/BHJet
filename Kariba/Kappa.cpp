#include "Particles.hpp"
#include "Kappa.hpp"

#include <iostream>

//Class constructor to initialize object
Kappa::Kappa(int s){
    size = s;

    p = new double[size];
    ndens = new double[size];	
    gamma = new double[size];
    gdens = new double[size];
    gdens_diff = new double[size];

    knorm = 1.;	

    mass_gr = emgm;
    mass_kev = emgm*gr_to_kev;

    for (int i=0;i<size;i++){
        p[i] = 0;
        ndens[i] = 0;
        gamma[i] = 0;
        gdens[i] = 0;
        gdens_diff[i] = 0;
    }
}

//Method to set the temperature, using ergs as input
void Kappa::set_temp_kev(double T){
    theta = T*kboltz_kev2erg/(mass_gr*cee*cee);

    double emin = (1./100.)*T;	                //minimum energy in kev, 1/100 lower than peak
    double emax = 20.*T; 		//maximum energy in kev, 20 higher than peak
    double gmin, gmax;

    gmin = emin/mass_kev+1.;
    gmax = emax/mass_kev+1.;
    pmin = pow(pow(gmin,2.)-1.,1./2.)*mass_gr*cee;
    pmax = pow(pow(gmax,2.)-1.,1./2.)*mass_gr*cee;	
}

void Kappa::set_kappa(double k){
    kappa = k;
}

//Methods to set momentum/energy arrays and number density arrays
void Kappa::set_p(double ucom,double bfield,double betaeff,double r,double fsc){	
    pmax = std::max(max_p(ucom,bfield,betaeff,r,fsc),pmax);	

    double pinc = (log10(pmax)-log10(pmin))/(size-1);

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
    }
}

//Same as above, but assuming a fixed maximum Lorentz factor
void Kappa::set_p(double gmax){
    pmax = pow(pow(gmax,2.)-1.,1./2.)*mass_gr*cee;

    double pinc = (log10(pmax)-log10(pmin))/(size-1);

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
    }	
}

void Kappa::set_ndens(){
    for (int i=0;i<size;i++){
        gdens[i] = knorm*gamma[i]*pow(pow(gamma[i],2.)-1.,1./2.)*pow(1.+(gamma[i]-1.)/(kappa*theta),-kappa-1.);
				
    }
    initialize_pdens();
    gdens_differentiate();	
}

//Methods to calculate the normalization of the function
double norm_kappa_int(double x,void *p){
    struct k_params * params = (struct k_params*)p; 	 

    double t = (params->t); 
    double k = (params->k);  	

    return  x*pow(pow(x,2.)-1.,1./2.)*pow(1.+(x-1.)/(k*t),-k-1.);
}

void Kappa::set_norm(double n){
    double norm_integral, error, min, max;

    min = pow(pow(pmin/(mass_gr*cee),2.)+1.,1./2.);
    max = pow(pow(pmax/(mass_gr*cee),2.)+1.,1./2.);

    gsl_function F1;	
    struct k_params params = {theta,kappa};
    gsl_integration_workspace *w1;    
    w1 = gsl_integration_workspace_alloc (100);
    F1.function = &norm_kappa_int;
    F1.params   = &params;
    gsl_integration_qag(&F1,min,max,0,1e-7,100,1,w1,&norm_integral,&error);
    gsl_integration_workspace_free (w1);

    knorm = n/norm_integral;
}

//Method to solve steady state continuity equation. NOTE: KN cross section not included in IC cooling
double injection_kappa_int(double x,void *p){
    struct injection_kappa_params * params = (struct injection_kappa_params*)p; 	 

    double t = (params->t); 
    double k = (params->k); 
    double n = (params->n);
    double m = (params->m);	

    double mom = pow(pow(x,2.)-1.,1./2.)*m*cee;
    double diff = mom/(pow(m*cee,2.)*pow(pow(mom/(m*cee),2.) +1.,1./2.));

    return  diff*n*x*pow(pow(x,2.)-1.,1./2.)*pow(1.+(x-1.)/(k*t),-k-1.);
}

void Kappa::cooling_steadystate(double ucom, double n0,double bfield,double r,double betaeff){
    double Urad = pow(bfield,2.)/(8.*pi)+ucom;
    double pdot_ad = betaeff*cee/r;
    double pdot_rad = (4.*sigtom*cee*Urad)/(3.*mass_gr*pow(cee,2.));
    double tinj = r/(cee);

    double integral, error;
    gsl_function F1;	
    struct injection_kappa_params params = {theta,kappa,knorm,mass_gr};
    F1.function = &injection_kappa_int;
    F1.params   = &params;

    for (int i=0;i<size;i++){
        if (i < size-1) {
            gsl_integration_workspace *w1;
            w1 = gsl_integration_workspace_alloc (100);
    	    gsl_integration_qag(&F1,gamma[i],gamma[i+1],1e1,1e1,100,1,w1,&integral,&error);
            gsl_integration_workspace_free (w1);

    	    ndens[i] = (integral/tinj)/(pdot_ad*p[i]/(mass_gr*cee)+pdot_rad*(gamma[i]*p[i]/(mass_gr*cee)));
        }
        else {
            ndens[size-1] = ndens[size-2]*pow(p[size-1]/p[size-2],-kappa);
        }
    }		
    // the last bin is set by arbitrarily assuming cooled distribution; this is necessary because the integral 
    //above is undefined for the last bin

    //The last step requires a renormalization. The reason is that the result of gsl_integration_qag strongly
    //depends on the value of "size". Without doing anything fancy, this can be fixed simply by ensuring that
    //the total integrated number of density equals n0 (which we know), and rescaling the array ndens[i] by
    //the appropriate constant.
    double renorm = count_particles()/n0;
	
    for (int i=0;i<size;i++){
        ndens[i] = ndens[i]/renorm;
    }

    initialize_gdens();
    gdens_differentiate();	
}

//Method to calculate maximum momentum of non thermal particles based on acceleration and cooling timescales
double Kappa::max_p(double ucom,double bfield,double betaeff,double r,double fsc){
    double Urad, escom, accon, syncon, b, c, gmax;
    Urad = pow(bfield,2.)/(8.*pi)+ucom;
    escom = betaeff*cee/r;
    syncon = (4.*sigtom*Urad)/(3.*mass_gr*cee);
    accon = (3.*fsc*charg*bfield)/(4.*mass_gr*cee);

    b = escom/syncon;
    c = accon/syncon;

    gmax = (-b+pow(pow(b,2.)+4.*c,1./2.))/2.;

    return pow(pow(gmax,2.)-1.,1./2.)*mass_gr*cee;
}

void Kappa::test(){
    std::cout << "Kappa distribution;" << std::endl;
    std::cout << "kappa index: " << kappa << std::endl;
    std::cout << "dimensionless temperature: " << theta << std::endl;
    std::cout << "Array size: " << size << std::endl;
    std::cout << "Default normalization: " << knorm << std::endl;
    std::cout << "Particle mass in grams: " << mass_gr << std::endl;
}
