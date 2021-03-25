#include "Thermal.hpp"

//Class constructor to initialize object
Thermal::Thermal(int s) {
    size = s;

    p = new double[size];
    ndens = new double[size];	
    gamma = new double[size];
    gdens = new double[size];
    gdens_diff = new double[size];
	
    thnorm = 1.;
    
    mass_gr = emgm;
    mass_kev = emgm*gr_to_kev;

    for (int i=0;i<size;i++){
        p[i] = 0;
        ndens[i] = 0;
    }
}

//Method to initialize momentum array a with default interval
void Thermal::set_p(){						        //
    double emin = (1./100.)*Temp;	    //minimum energy in kev, 1/100 lower than peak
    double emax = 20.*Temp; 		    //maximum energy in kev, 20 higher than peak
    double gmin, gmax, pmin, pmax, pinc;

    gmin = emin/mass_kev+1.;
    gmax = emax/mass_kev+1.;
   
    pmin = pow(pow(gmin,2.)-1.,1./2.)*mass_gr*cee;
    pmax = pow(pow(gmax,2.)-1.,1./2.)*mass_gr*cee;	
    pinc = (log10(pmax)-log10(pmin))/(size-1);

    for (int i=0;i<size;i++){
        p[i] = pow(10.,log10(pmin)+i*pinc);
        gamma[i] = pow(pow(p[i]/(mass_gr*cee),2.)+1.,1./2.);
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
void Thermal::set_temp_kev(double T){
    Temp = T;
    theta = (Temp*kboltz_kev2erg)/(mass_gr*cee*cee);
}

void Thermal::set_norm(double n){
    thnorm = n/(pow(mass_gr*cee,3.)*theta*K2(1./theta));
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
    std::cout << "Particle mass in grams: " << mass_gr <<std:: endl;
    std::cout << "Particle mass in keV: " << mass_kev <<std:: endl;
    std::cout << "kT/mc^2: " << theta << std::endl;	
}
