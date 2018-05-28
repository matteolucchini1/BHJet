/* -----------------------------------------
 *
 * From Maitra's fortran version of agnjet
 * Described in Maitra et al. 2009
 * Adapted version by Chiara Ceccobello
 * -----------------------------------------
 */
#include "agnjet.hh"
using namespace std;

void pproduction(int ncom,double r,double ephxr[],double totcom[],double &n_pp,double &rate_pp){
    //Dump local IC-photon-number-density-spectrum at jet base 
    
    for (int ii=0; ii<ncom; ii++){
        // cDM Get n(alpha) at alpha~1 and roughly estimate PP rate
        if (ii>2 &&((ephxr[ii-1] < 1.23559E+20) && (ephxr[ii] > 1.23559E+20))) {
            n_pp=log10(totcom[ii])-6.1+log10(r/cee);
            rate_pp = 2.0*n_pp-15.0;
            
            cout << "log_10 [ PP rate (1/s/cm^3) ] ~ " << rate_pp << endl;
            cout << "log_10 [ PP number density (1/cm^3) ] ~ " << n_pp << endl;
        }
    }
}


void annihilation(double z,double zcut,double r,double ntot,double eltemp,double &rate_pa,double &n_pa){
    //Estimate pair annihilation rate for thermal plasma at base */
    if(z < zcut){
        pa_rate_thermal (ntot, eltemp, rate_pa);
        n_pa = rate_pa-log10(r/cee);
        cout << "log_10 [ ntot (1/cm^3) ]     = " << log10(ntot) << endl;
        cout << "log_10 [ PA rate (1/s/cm^3) ] = " << rate_pa << endl;
        cout << "log_10 [ PA number density (1/s/cm^3) ] = " << n_pa << endl;
    } else {
        n_pa = 1.0;
        rate_pa = 1.0;
    }
}

/* 
 * Pair annihilation rate (eq.3.7 of Coppi & Blandford 1990)
 */
void pa_rate(double gamma_plus,double gamma_minus,double &rate) {
    double constant,x;

    constant=7.48125E-15;
    x=gamma_plus*gamma_minus;
    rate=constant*(log(x)+1./sqrt(x))/x;
}

/*
 *  Pair annihilation rate for one temperature thermal plasma (eq.68 of Svensson 1982)
 */
void pa_rate_thermal(double numden,double T_e,double &rate){
    
    double pi_c_re_re, k_mcc, theta, denom;

    pi_c_re_re=7.48433E-15;
    k_mcc=1.6863E-10;
    theta=k_mcc*T_e;
    denom=1.0 + 2.0*theta*theta/(log(1.12291*theta+1.3));
    rate=log10(numden*numden*pi_c_re_re/denom);
}
