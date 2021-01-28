#include "Thermal.hpp"

//Class constructor to initialize object
Thermal::Thermal(int s,int type, double T) {
    size = s;

    p = new double[size];
    ndens = new double[size];	
    gamma = new double[size];
    gdens = new double[size];
    gdens_diff = new double[size];

    w1 = gsl_integration_workspace_alloc (100);
	
    if (type==1) {mass = emgm;}
    else  {mass = pmgm;}
	
    Temp = T;
    theta = Temp/(mass*cee*cee);
    thnorm = 1.;

    for (int i=0;i<size;i++){
        p[i] = 0;
        ndens[i] = 0;
    }
}

//Method to initialize momentum array a with default interval
void Thermal::set_p(){						//TODO: check that the proton case is handled correctly
    double emin = (1./100.)*(Temp/kboltz_kev2erg);	//minimum energy in kev, 1/50 lower than peak
    double emax = 20.*(Temp/kboltz_kev2erg); 		//maximum energy in kev, 20 higher than peak
    double gmin, gmax, pmin, pmax, pinc;

    if (mass == emgm){
        gmin = emin/me_kev+1.;
        gmax = emax/me_kev+1.;
    } else if (mass == pmgm) {
        std::cout << "Error! Protons currently unsupported!" <<std:: endl;
    }

    pmin = pow(pow(gmin,2.)-1.,1./2.)*mass*cee;
    pmax = pow(pow(gmax,2.)-1.,1./2.)*mass*cee;	
    pinc = (log10(pmax)-log10(pmin))/size;

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass*cee),2.)+1.,1./2.);
    }	
}

//Method to set differential electron number density from known temperature, normalization, and momentum array
void Thermal::set_ndens(){
    for (int i=0;i<size;i++){
        ndens[i] = thnorm*pow(p[i],2.)*exp(-gamma[i]/theta);
    }
    initialize_gdens();
    gdens_differentiate();
}	

//methods to set the temperature and normalization. NOTE: temperature must be in ergs, no factor kb
void Thermal::set_temp(double T){
    Temp = T;
    theta = Temp/(mass*cee*cee);
}

void Thermal::set_norm(double n){
    thnorm = n/(pow(mass*cee,3.)*theta*K2(1./theta));
}

//Evaluate Bessel function as in old agnjet
double Thermal::K2(double x){
    double res;

    if (x < 0.1){
        res = 2./x/x;
    }
    else {
        res = gsl_sf_bessel_Kn(2,x);
    }

    return res;
}

//simple method to check quantities.
void Thermal::test(){
    std::cout << "Thermal distribution;" << std::endl;
    std::cout << "Temperature: " << Temp << " erg, " << Temp/kboltz_kev2erg << " kev" << std::endl;
    std::cout << "Array size: " << size <<std:: endl;
    std::cout << "Normalization: " << thnorm << std::endl;
    std::cout << "Particle mass: " << mass <<std:: endl;
    std::cout << "kT/mc^2: " << theta << std::endl;	
}
