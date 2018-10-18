#include "agnjet.hh"

using namespace std;

void synchrotron(bool &isBreaknjetnsyn,int nsyn,int njet,int nz,int k,double bfield,double delz,double inclin,
double dist,double r,double elenmn,double elenmx,gsl_spline *spline_syn, gsl_interp_accel *acc_syn, 
gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_derivs, 
gsl_interp_accel *acc_derivs,double nurad[],double dopfac[],double nphot[],double nusyn[],double synabs[]){

    int i,m,j,l;
    double esum,asum,phoden,absd;
    
    //4.\times \pi in synabs is for assumed isotropic source. Synint returns units of ergs/cm^2/s/Hz for 
    //absorbed flux; need to convert to Jansky
    for(i=0; i<nsyn; i++){
        for(m=0;m<njet; m++){
            nusyn[(m*nz+k)*nsyn+i]= log10(nurad[i]*dopfac[m*nz+k]);           
        }
    }
    isBreaknjetnsyn = false;

    for(i=0; i<nsyn; i++){
        // First performing integrations for synchrotron emissivity and absorption at current frequency. 
        synintegrals(nurad[i],bfield,elenmn,elenmx,spline_syn,acc_syn,spline_eldis1,acc_eldis1,spline_derivs,
        acc_derivs,esum,asum);
        
        for(m=0; m<njet; m++){
            synint(nurad[i],r,inclin,dopfac[m*nz+k],delz,bfield,dist,esum,asum,absd,phoden);            
            nphot[i]= phoden;        
            if(absd < 1.e-100){
                for(j=i; j<nsyn; j++){ //CHECK boundaries loop !!!
                    for(l=0; l<njet; l++){
                        synabs[(l*nz+k)*nsyn+j]=-100.;
                    }
                    nphot[j]= 0;
                }
                isBreaknjetnsyn = true;
                break; //out of both m(njet) and i(nsyn) loops. Needs isBreak flag.
            }
            synabs[(m*nz+k)*nsyn+i]=log10(absd*dopfac[m*nz+k]*dopfac[m*nz+k]/mjy);
        }//END for-loop njet
        if(isBreaknjetnsyn){
            break;
        }
    }//END for-loop nsyn
}

/**
 * This function is for synchrotron emissivity to be integrated over particle distribution in synint
 */
double synemis(double e1e, void *p){
    struct synemis_params *params = (struct synemis_params *)p;
    double freq     = (params -> freq);
    double bfield   = (params -> bfield);
    gsl_spline *spline_syn = (params -> spline_syn);
    gsl_interp_accel *acc_syn = (params -> acc_syn);
    gsl_spline *spline_eldis1 = (params -> spline_eldis1);
    gsl_interp_accel *acc_eldis1 = (params -> acc_eldis1);
    
    double ele, game, x, emisfunc,  ellog;
    double emisfuncl= 0;
    double eden	= 0;
    
    ele	= exp(e1e);
    game	= ele/emerg;
    x	= freq*4.*pi*emgm*cee/(3.*charg*bfield*game*game);
    
    if(x <= 1.e-4){
        emisfunc	= 4.*pi*pow(x/2.,(1./3.))/(sqrt(3.)*2.68);
    }
    else if(x > 50.){
        emisfunc	= sqrt(pi*x/2.)*exp(-x);
    }
    else{
        emisfuncl       = gsl_spline_eval(spline_syn, x, acc_syn);
        emisfunc	= pow(10,emisfuncl);
    }
    
    ellog	= e1e * 0.434294481903252; //0.434294481903252 = log10(exp(1))
    eden    = gsl_spline_eval(spline_eldis1, ellog, acc_eldis1);
    eden	= pow(10, eden);
    
    //gsl_spline_free(spline_syn), gsl_interp_accel_free(acc_syn);
    //gsl_spline_free(spline_eldis1), gsl_interp_accel_free(acc_eldis1);
    return ele*eden*emisfunc;
}

/**
 * This function is for self-absorption to be intergrated over the derived function of the particle 
 * distribution in synint
 */
double absfnc(double e2, void *p){
    struct absfnc_params *params = (struct absfnc_params *)p;
    double freq     = (params -> freq);
    double bfield   = (params -> bfield);
    gsl_spline *spline_syn = (params -> spline_syn);
    gsl_interp_accel *acc_syn = (params -> acc_syn);
    gsl_spline *spline_derivs = (params -> spline_derivs);
    gsl_interp_accel *acc_derivs = (params -> acc_derivs);
    
    double elec, game, x, absfunc;
    double absfuncl = 0;
    double deden	= 0;
    
    elec	= exp(e2);
    game	= elec/emerg;
    x	= freq*4.*pi*emgm*cee/(3.*charg*bfield*game*game);
    
    if(x <= 1.e-4){
        absfunc	= 4.*pi*pow(x/2.,1./3.)/(sqrt(3.)*2.68);
    }
    else if(x > 50.){
        absfunc	= sqrt(pi*x/2.)*exp(-x);
    }
    else{
        absfuncl= gsl_spline_eval(spline_syn, x, acc_syn);
        absfunc	= pow(10, absfuncl);
    }
    
    deden   = gsl_spline_eval(spline_derivs, elec, acc_derivs);
    
    //gsl_spline_free(spline_syn), gsl_interp_accel_free(acc_syn);
    //gsl_spline_free(spline_derivs), gsl_interp_accel_free(acc_derivs);
    return absfunc*(elec*elec*elec)*deden;
}

/**
 * Synintegrals does frequency and electron spectrum dependent calc, single component. Returns necessary for 
 * emissivity and self-absorption for use in synint, which calculates the actual spectrum.
 * @param freq	Frequency of emission
 *
 * @return esum	Synchrotron emission kernal at frequency 'freq'
 * @return asum Synchrotron absorption kernal at frequency 'freq'
 */
void synintegrals(double freq, double bfield, double elenmn, double elenmx, gsl_spline *spline_syn, 
gsl_interp_accel *acc_syn, gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_derivs, 
gsl_interp_accel *acc_derivs, double &esum, double &asum){ 
    //esum
    gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(1000);
    double result1, error1;
    gsl_function F1;
    struct synemis_params F1params = {freq, bfield, spline_syn, acc_syn, spline_eldis1, acc_eldis1};
    F1.function     = &synemis;
    F1.params       = &F1params;
    gsl_integration_qag(&F1, log(elenmn), log(elenmx), 1e-1, 1e-1, 1000, 2, w1, &result1, &error1);
    esum    = result1;
    gsl_integration_workspace_free(w1);     
    
    //asum
    gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(1000);
    double result2, error2;
    gsl_function F2;
    struct absfnc_params F2params = {freq, bfield, spline_syn, acc_syn, spline_derivs, acc_derivs};
    F2.function     = &absfnc;
    F2.params       = &F2params;
    gsl_integration_qag(&F2, log(elenmn), log(elenmx), 1e-1, 1e-1, 1000, 2, w2, &result2, &error2);
    asum    = result2;
    gsl_integration_workspace_free(w2);
}

/**
 * For selfabsorption, single component returns specific intensity in erg/s/cm^2/st/Hz integrates 
 * distributions over electron spectrum. This takes sin(pitch) = 2/3 equiv to angle averaging over isotropic
 *
 * @param freq		Frequency of emission
 * @param r			Radius of the zone
 * @param angle		Source viewing angle
 * @param dopfac	Doppler factor in the zone
 * @param h			Height in the jet along z-axis
 * @param bfield    Magnetic field magnitude
 * @param dist      Distance to the source
 * @param esum		Synchrotron emission kernel
 * @param asum		Synchrotron absportion kernel
 *
 * @return fluxa	Synchrotron self-absorption emission
 * @return fluxa2	Resulting photon density
 *
 */
void synint(double freq, double r, double angle, double dopfac, double h, double bfield, double dist, 
double &esum, double &asum, double &fluxa, double &fluxa2){
    double pitch = 0.73;
    double area, acons, tsyn;
    double absfac, tsyn2, absfac2;
    double asyn;
    double elcons;
    double epsasyn;    
    
    elcons	= sqrt(3.)*(charg*charg*charg)*bfield*sin(pitch)/emerg;
    acons	= -cee*cee/(8.*pi*freq*freq);
    asyn	= acons*elcons*asum;
    epsasyn	= esum/(acons*asum);
    area	= 4.*pi*dist*dist;
    
    //This first term is for what is seen externally. It includes skin depth and angle effects. See pg 12,
    //ntbk 4 for consideration.
    tsyn	= pi/2. * asyn*r/(dopfac*sin(angle));
    if(tsyn >= 1.){
        absfac	= (1.-exp(-tsyn));
    }
    else{
        absfac	= tsyn*(1.+tsyn*(-0.5+tsyn*(1./6.+tsyn*(-1./24.+tsyn/120))));
    }    
    fluxa	= (2.*r*h*sin(angle)*dopfac*absfac*epsasyn)/area;
    
    //This second term is the same as above, but assuming what is "seen" locally by the particles for compton
    //scattering. Pass 'specific luminosity', erg/s/Hz to main code
    tsyn2	= pi/2.*asyn*r;
    if(tsyn2 >= 1.){
        absfac2	= (1.-exp(-tsyn2));
    }
    else{
        absfac2	= tsyn2*(1.+tsyn2*(-0.5+tsyn2*(1./6.+tsyn2*(-1./24.+tsyn2/120))));
    }    
    fluxa2	= pi*r*r*absfac2*epsasyn;
}
