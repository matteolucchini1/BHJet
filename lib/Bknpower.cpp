#include "Particles.hpp"
#include "Bknpower.hpp"

#include <iostream>

Bknpower::Bknpower(int s){
    size = s;

    p = new double[size];
    ndens = new double[size];	
    gamma = new double[size];
    gdens = new double[size];
    gdens_diff = new double[size];

    norm = 1.;

    for (int i=0;i<size;i++){
        p[i] = 0;
        ndens[i] = 0;
        gamma[i] = 0;
        gdens[i] = 0;
        gdens_diff[i] = 0;
    }
}

//Methods to set momentum/energy arrays
void Bknpower::set_p(double min,double brk,double ucom,double bfield,double betaeff,double r,double fsc){	
    pmin = min;
    pbrk = brk;
    pmax = max_p(ucom,bfield,betaeff,r,fsc);	

    double pinc = (log10(pmax)-log10(pmin))/size;

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
    }
}

void Bknpower::set_p(double min,double brk,double gmax){
    pmin = min;
    pbrk = brk;
    pmax = pow(pow(gmax,2.)-1.,1./2.)*mass_gr*cee;

    double pinc = (log10(pmax)-log10(pmin))/size;

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
    }	
}

//Method to set differential electron number density from known pspec, normalization, and momentum array
void Bknpower::set_ndens(){
    for (int i=0;i<size;i++){
        ndens[i] = norm*pow(p[i]/pbrk,-pspec1)/(1.+pow(p[i]/pbrk,-pspec1+pspec2))*exp(-p[i]/pmax);
    }
    initialize_gdens();
    gdens_differentiate();		
}

//methods to set the slopes, break and normalization
void Bknpower::set_pspec1(double s1){
    pspec1 = s1;
}

void Bknpower::set_pspec2(double s2){
    pspec2 = s2;
}

void Bknpower::set_brk(double brk){
    pbrk = brk;
}

//Methods to calculate the normalization of the function
double norm_bkn_int(double x,void *p){	
    struct bkn_params * params = (struct bkn_params*)p; 	 

    double s1 = (params->s1);
    double s2 = (params->s2);
    double brk = (params->brk);
    double max = (params->max);
    double m = (params->m);	

    double mom_int = pow(pow(x,2.)-1.,1./2.)*m*cee;	

    return pow(mom_int/brk,-s1)/(1.+pow(mom_int/brk,-s1+s2))*exp(-mom_int/max);
}

void Bknpower::set_norm(double n){
    double norm_integral, error, min, max;

    min = pow(pow(pmin/(mass_gr*cee),2.)+1.,1./2.);
    max = pow(pow(pmax/(mass_gr*cee),2.)+1.,1./2.);

    gsl_integration_workspace *w1;
    w1 = gsl_integration_workspace_alloc (100);
    gsl_function F1;	
    struct bkn_params params = {pspec1,pspec2,pbrk,pmax,mass_gr};
    F1.function = &norm_bkn_int;
    F1.params   = &params;
    gsl_integration_qag(&F1,min,max,0,1e-7,100,1,w1,&norm_integral,&error);
    gsl_integration_workspace_free (w1);

    norm = n/(norm_integral*mass_gr*cee);
}

//Injection function to be integrated in cooling
double injection_bkn_int(double x,void *p){	
    struct injection_bkn_params * params = (struct injection_bkn_params*)p; 	 

    double s1 = (params->s1);
    double s2 = (params->s2);
    double brk = (params->brk);
    double max = (params->max);
    double m = (params->m);	
    double n = (params->n);

    double mom_int = pow(pow(x,2.)-1.,1./2.)*m*cee;	

    return n*pow(mom_int/brk,-s1)/(1.+pow(mom_int/brk,-s1+s2))*exp(-mom_int/max);
}


//Method to solve steady state continuity equation. NOTE: KN cross section not included in IC cooling
void Bknpower::cooling_steadystate(double ucom, double n0,double bfield,double r,double betaeff){
    double Urad = pow(bfield,2.)/(8.*pi)+ucom;
    double pdot_ad = betaeff*cee/r;
    double pdot_rad = (4.*sigtom*cee*Urad)/(3.*mass_gr*pow(cee,2.));
    double tinj = r/(cee);

    double integral, error;
    gsl_function F1;	
    struct injection_bkn_params params = {pspec1,pspec2,pbrk,pmax,mass_gr,n0};
    F1.function = &injection_bkn_int;
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
            ndens[size-1] = ndens[size-2]*pow(p[size-1]/p[size-2],-pspec2)*exp(-1.);
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
double Bknpower::max_p(double ucom,double bfield,double betaeff,double r,double fsc){
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
void Bknpower::test(){
    std::cout << "Broken power-law distribution;" << std::endl;
    std::cout << "pspec1: " << pspec1 << std::endl;
    std::cout << "pspec2: " << pspec2 << std::endl;
    std::cout << "pbreak: " << pbrk << std::endl;
    std::cout << "Array size: " << size << std::endl;
    std::cout << "Default normalization: " << norm << std::endl;
    std::cout << "Particle mass: " << mass_gr << std::endl;
}
