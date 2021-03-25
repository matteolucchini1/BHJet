#include "BBody.hpp"

BBody::~BBody(){
    delete[] en_phot_obs;
    delete[] num_phot_obs;
}

BBody::BBody(){
    size = 40;

    en_phot = new double[size];
    num_phot = new double[size];
    en_phot_obs = new double[size];
    num_phot_obs = new double[size];

    for(int i=0;i<size;i++){
        en_phot[i] = 0.;
        num_phot[i] = 0.;  
        en_phot_obs[i] = 0.;
        num_phot_obs[i] = 0.;       
    }	
}

//Methods to set BB quantities
void BBody::set_temp_kev(double T){
    double emin,emax,einc;

    Tbb = T*kboltz_kev2erg/kboltz;

    emin = 0.02*kboltz*Tbb;
    emax = 30.*kboltz*Tbb;

    einc = (log10(emax)-log10(emin))/(size-1);

    for(int i=0;i<size;i++){
        en_phot[i] = pow(10.,log10(emin)+i*einc);
        en_phot_obs[i] = en_phot[i];
    }	
}

void BBody::set_temp_k(double T){
    double emin,emax,einc;

    Tbb = T;

    emin = 0.02*kboltz*Tbb;
    emax = 30.*kboltz*Tbb;

    einc = (log10(emax)-log10(emin))/(size-1);

    for(int i=0;i<size;i++){
        en_phot[i] = pow(10.,log10(emin)+i*einc);
        en_phot_obs[i] = en_phot[i];
    }	
}

void BBody::set_temp_hz(double nu){
    double emin,emax,einc;

    Tbb = (herg*nu)/(2.82*kboltz);

    emin = 0.02*kboltz*Tbb;
    emax = 30.*kboltz*Tbb;

    einc = (log10(emax)-log10(emin))/(size-1);

    for(int i=0;i<size;i++){
        en_phot[i] = pow(10.,log10(emin)+i*einc);
        en_phot_obs[i] = en_phot[i];
    }	
}

void BBody::set_lum(double L){
    Lbb = L;
    normbb = Lbb/(sbconst*pow(Tbb,4.));
}

//Method to set BB spectrum
void BBody::bb_spectrum(){
    for (int i=0;i<size;i++){
	    num_phot[i] = normbb*2.*herg*pow(en_phot_obs[i]/herg,3.)/
                      (pow(cee,2.)*(exp(en_phot_obs[i]/(Tbb*kboltz))-1.));
        num_phot_obs[i] = num_phot[i];
    }
}


//Methods to return BB temperature, luminosity, energy density at a given distance d (or for a given radius d 
//of the source)
const double BBody::temp_kev(){
    return Tbb*kboltz/kboltz_kev2erg;
}

const double BBody::temp_k(){
    return Tbb;
}

const double BBody::temp_hz(){
    return 2.82*kboltz*Tbb/herg;
}

const double BBody::lum(){
    return Lbb;
}

const double BBody::norm(){
    return normbb;
}

const double BBody::Urad(double d){
    return Lbb/(4.*pi*pow(d,2.)*cee);
}

void BBody::test(){
    std::cout << "Black body temperature in K: " << temp_k() << " in keV: " << temp_kev() << std::endl;
    std::cout << "Black body bolometric luminosity in erg s^-1 " << Lbb << std::endl;
}
