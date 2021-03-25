#include "Particles.hpp"
#include "Powerlaw.hpp"

#include <iostream>

//Class constructor to initialize object
Powerlaw::Powerlaw(int s){
    size = s;

    p = new double[size];
    ndens = new double[size];	
    gamma = new double[size];
    gdens = new double[size];
    gdens_diff = new double[size];

    plnorm = 1.;
    
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

//Methods to set momentum/energy arrays
void Powerlaw::set_p(double min,double ucom,double bfield,double betaeff,double r,double fsc){	
    pmin = min;
    pmax = max_p(ucom,bfield,betaeff,r,fsc);	

    double pinc = (log10(pmax)-log10(pmin))/(size-1);

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
    }
}

void Powerlaw::set_p(double min,double gmax){
    pmin = min;
    pmax = pow(pow(gmax,2.)-1.,1./2.)*mass_gr*cee;
    
    double pinc = (log10(pmax)-log10(pmin))/(size-1);

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
    }	
}

//Method to set differential electron number density from known pspec, normalization, and momentum array
void Powerlaw::set_ndens(){
    for (int i=0;i<size;i++){
        ndens[i] = plnorm*pow(p[i],-pspec)*exp(-p[i]/pmax);
    }
    initialize_gdens();
    gdens_differentiate();		
}	

//methods to set the slope and normalization
void Powerlaw::set_pspec(double s1){
    pspec = s1;
}

void Powerlaw::set_norm(double n){
    plnorm = n*(1.-pspec)/(pow(p[size-1],(1.-pspec))-pow(p[0],(1.-pspec)));
}

//Injection function to be integrated in cooling
double injection_pl_int(double x,void *p){	
    struct injection_pl_params * params = (struct injection_pl_params*)p; 	
    double s = (params->s);
    double n = (params->n);
    double m = (params->m);
    double max = (params->max);

    double mom_int = pow(pow(x,2.)-1.,1./2.)*m*cee;	

    return n*pow(mom_int,-s)*exp(-mom_int/max);
}

//Method to solve steady state continuity equation. NOTE: KN cross section not included in IC cooling
void Powerlaw::cooling_steadystate(double ucom, double n0,double bfield,double r,double betaeff){
    double Urad = pow(bfield,2.)/(8.*pi)+ucom;
    double pdot_ad = betaeff*cee/r;
    double pdot_rad = (4.*sigtom*cee*Urad)/(3.*mass_gr*pow(cee,2.));
    double tinj = r/(cee);

    double integral, error;
    gsl_function F1;	
    struct injection_pl_params params = {pspec,plnorm,mass_gr,pmax};
    F1.function = &injection_pl_int;
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
            ndens[size-1] = ndens[size-2]*pow(p[size-1]/p[size-2],-pspec-1)*exp(-1.);
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
//The estimate is identical to the old agnjet but in momentum space; see Lucchini et al. 2019 for the math of
//the old version
double Powerlaw::max_p(double ucom,double bfield,double betaeff,double r,double fsc){
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

//simple method to check quantities.
void Powerlaw::test(){
    std::cout << "Power-law distribution;" << std::endl;
    std::cout << "pspec: " << pspec << std::endl;
    std::cout << "Array size: " << size << std::endl;
    std::cout << "Default normalization: " << plnorm << std::endl;
    std::cout << "Particle mass in grams: " << mass_gr << std::endl;
}

