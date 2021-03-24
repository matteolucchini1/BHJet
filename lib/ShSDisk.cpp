#include "ShSDisk.hpp"

ShSDisk::~ShSDisk(){
    delete[] en_phot;
    delete[] num_phot;
    delete[] en_phot_obs;	
    delete[] num_phot_obs;
}

//Constructor for the disk. No input parameters because the sized of the arrays is set automatically, and 
//every other property needs to be handled by the setters below
ShSDisk::ShSDisk(){
    size = 50;

    en_phot = new double[size];
    num_phot = new double [size];
    en_phot_obs = new double [size];
    num_phot_obs = new double [size];
    
    for (int i=0;i<size;i++) {
        en_phot[i] = 0.;
        en_phot_obs[i] = 0.;
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    } 
}	

//return SD spectrum over a given radius, frequency to be integrated over radius
double disk_int(double lr,void *p){
    struct disk_obs_params *params = (struct disk_obs_params *)p;
    double tin      = (params -> tin);
    double rin      = (params -> rin);
    double nu       = (params -> nu);

    double r,temp,fac,bb;
    r = exp(lr);
    temp = tin*pow(rin/r,0.75);
    fac	= herg*nu/(kboltz*temp);

    if(fac<1.e-3){
        bb	= 2.*herg*pow(nu,3.)/(pow(cee,2.)*fac);
    }
    else{
        bb	= 2.*herg*pow(nu,3.)/(pow(cee,2.)*(exp(fac)-1.));
    }

    return 2.*pi*pow(r,2.)*bb;
}

void ShSDisk::disk_spectrum(){
    double result, error;

    for(int k=0;k<size;k++){
        gsl_integration_workspace *w1;        
        w1 = gsl_integration_workspace_alloc (100);
        gsl_function F1;
        struct disk_obs_params F1params = {Tin,r,en_phot_obs[k]/herg};
        F1.function = &disk_int;
        F1.params = &F1params;
        gsl_integration_qag(&F1,log(r),log(z),0,1e-2,100,2,w1,&result,&error);
        gsl_integration_workspace_free (w1);
        
        num_phot[k] = result;
        num_phot_obs[k] = cos(angle)*result;        
    }
}

//this method removes a fraction f from the observed disk luminosity, assuming that it was absorbed and 
//reprocessed into some other unspecified radiative mechanism
void ShSDisk::cover_disk(double f){
    for (int k=0;k<size;k++) {
        num_phot_obs[k] = num_phot_obs[k]*(1.-f);
    }
}

//Simple integration method to integrate num_phot_obs and get the total luminosity of the disk
double ShSDisk::total_luminosity(){
    double temp = 0.;
    for (int i=0;i<size-1;i++){
        temp = temp+(1./2.)*(en_phot_obs[i+1]/herg-en_phot_obs[i]/herg)*
               (num_phot_obs[i+1]/cos(angle)+num_phot_obs[i]/cos(angle));
    }
    return temp;
}

void ShSDisk::set_mbh(double M){
    Mbh = M;
    Rg = gconst*Mbh*msun/(cee*cee);
}

//note: for disk r is inner radius, z is the outer radius, both are input in cm
void ShSDisk::set_rin(double R){
    r = R;					
}

void ShSDisk::set_rout(double R){
    z = R;					
}

//NOTE: the constant factor to go between inner temperature Tin and disk lunminosity Ldisk is
//2 because the model uses the zero torque inner boundary condition, Kubota et al. 1998, hence the factor 2
//rather than 4pi when converting between luminosity and temperature
void ShSDisk::set_luminosity(double L){
    double emin,emax,einc;

    Ldisk = L;
    Tin =  pow(Ldisk*1.25e38*Mbh/(2.*sbconst*pow(r,2.)),0.25);
    Hratio = std::max(0.1,Ldisk);    
    emin = 0.0001*kboltz*Tin;
    emax = 30.*kboltz*Tin;
    einc = (log10(emax)-log10(emin))/(size-1);

    for(int i=0;i<size;i++){
        en_phot_obs[i] = pow(10.,log10(emin)+i*einc);
        en_phot[i] = en_phot_obs[i];
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    }
}

void ShSDisk::set_tin_kev(double T){
    double emin,emax,einc;

    //note: 1 keV = kboltz_kev2erg/kboltz keV
    Tin = T*kboltz_kev2erg/kboltz;
    Ldisk = 2.*sbconst*pow(Tin,4.)*pow(r,2.)/(1.25e38*Mbh);    
    Hratio = std::max(0.1,Ldisk);
    emin = 0.0001*kboltz*Tin;
    emax = 30.*kboltz*Tin;
    einc = (log10(emax)-log10(emin))/(size-1);

    for(int i=0;i<size;i++){
        en_phot_obs[i] = pow(10.,log10(emin)+i*einc);
        en_phot[i] = en_phot_obs[i];
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    }
}

void ShSDisk::set_tin_k(double T){
    double emin,emax,einc;    
    
    Tin = T;
    Ldisk = 2.*sbconst*pow(Tin,4.)*pow(r,2.)/(1.25e38*Mbh);    
    Hratio = std::max(0.1,Ldisk);
    emin = 0.0001*kboltz*Tin;
    emax = 30.*kboltz*Tin;
    einc = (log10(emax)-log10(emin))/(size-1);

    for(int i=0;i<size;i++){
        en_phot_obs[i] = pow(10.,log10(emin)+i*einc);
        en_phot[i] = en_phot_obs[i];
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    }
}


void ShSDisk::test(){
    std::cout << "Inner disk radius: " << r << " cm, " << r/Rg << " rg; outer disk radius: " << z <<
                 " cm, " << z/Rg << " rg; disk scale height: " << Hratio << std::endl;	
    std::cout << "Inner disk temperature: " << Tin*kboltz/kboltz_kev2erg << 
                 " kev; emitted  disk luminosity in Eddington units: " << Ldisk << " and erg s^-1: " << 
                 total_luminosity() << std::endl;
}
