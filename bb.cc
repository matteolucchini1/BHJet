#include "agnjet.hh"

using namespace std;

void ec_comp(int disksw,int compsw, double z,double tbbeff,double normbb,double rbb,double tin,double rin,
double rout,double hbb,double gamv,double &ucom,double &uphdil,double &uphdil2){
    
    double nubot,nutop,tst,bbint;
    
    //Estimate disk contribution to cooling
    if(disksw==1){
		nubot	= 1.5e-5*kboltz*tin/herg;
        nutop	= 8.*kboltz*tin/herg;
        if(bbdisk(nutop,gamv,tin,rin,rout,z,hbb)*nubot/(bbdisk(nubot,gamv,tin,rin,rout,z,hbb)*nutop) > 1.e30){
			cerr << "nutop not high enough, something is wrong with the disk" << endl;
		}
        tst	= bbdisk(nutop,gamv,tin,rin,rout,z,hbb);
        if(tst <= 1.e-20){
        	while(tst < 1.e-20){
            	nutop	*= 0.75;
                tst	= bbdisk(nutop,gamv,tin,rin,rout,z,hbb);
            }
		}
        bbintegrals(log(nubot),log(nutop),gamv,tin,rin,rout,z,hbb,bbint);
        uphdil	= bbint/cee; //this has units of erg cm^-3, it's the disk integrated energy density
		ucom	= ucom + uphdil;     
    }//END if(disksw==1)

    if (compsw==1){
 		uphdil2 = (pow(gamv,2)*normbb*sbconst*pow(tbbeff/gamv,4))/(4.*pi*pow(rbb,2)*cee);
 		ucom = ucom + uphdil2;
    }        
}

double bbfnc(double thet, void *p){
    struct bbfnc_params *params = (struct bbfnc_params *)p;
    double gamv     = (params -> gamv);
    double z        = (params -> z);
    double hbb      = (params -> hbb);
    double tin      = (params -> tin);
    double rin      = (params -> rin);
    double nu       = (params -> nu);
    
    double beta, r, tineff, fac, temp, bbfunc;
    
    //This is what the jet sees from BB in stationary frame. It's the effective temperature as input.
    beta	= sqrt(gamv*gamv-1.)/gamv;
    r	= (z - hbb/2.)*tan(thet);
    tineff	= tin/(gamv*(1.+beta*cos(thet)));
    temp	= tineff*pow(rin/r, 0.75);
    fac	= herg*nu/(kboltz*temp);
    
    if(fac < 1e-3){
        bbfunc	= sin(thet)*cos(thet)/fac;
    }
    else{
        bbfunc	= sin(thet)*cos(thet)/(exp(fac)-1.);
    }
    
    //Since we are interested only in the photon energy density magnitude and not direction,F dot n is 
    //equivalent to absolute value
    return abs(bbfunc);
}

/**
 * bb_func3 is for inside the disk
 */
double bbfnc3(double thet, void *p){
    struct bbfnc3_params *params = (struct bbfnc3_params *)p;
    double gamv     = (params -> gamv);
    double tin      = (params -> tin);
    double nu       = (params -> nu);
    
    double beta, tineff, fac, bbfunc3;
    
    //This is what the jet sees from BB in stationary frame. see pg 256, iv and R&L pg 153. It's the effective
    //temperature as input.
    beta	= sqrt(gamv*gamv-1.)/gamv;
    tineff	= tin/(gamv*(1.+beta*cos(thet)));
    fac	= herg*nu/(kboltz*tineff);
    
    if(fac < 1.e-3){
        bbfunc3	= sin(thet)*cos(thet)/fac;
    }
    else{
        bbfunc3	= sin(thet)*cos(thet)/(exp(fac)-1.);
    }
    
    //Since we are interested only in the photon energy density magnitude and not direction,F dot n is 
    //equivalent to absolute value
    bbfunc3	= abs(bbfunc3);
    
    return bbfunc3;
}

/**
 * bbdisk gives flux of the top of the disk as a function of the frequency and assumes that hbb is constant, 
 * set by inner radius
 */
double bbdisk(double nu,double gamv,double tin,double rin,double rout,double z,double hbb){
    double blim,ulim,inside,top;
    
    blim	= atan(rin/(z+hbb/2.));
    if(z > hbb/2.){
        ulim	= atan(rin/(z-hbb/2.));
    }
    else{
        ulim	= pi/2. + atan((hbb/2-z)/rin);
    }

    gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(1000);
    double result2, error2;
    gsl_function F2;
    struct bbfnc3_params F2params   = {gamv, tin, nu};
    F2.function     = &bbfnc3;
    F2.params       = &F2params;
    gsl_integration_qag(&F2, blim, ulim, 1e-1, 1e-1, 1000, 1, w2, &result2, &error2);
    inside  = result2;
    gsl_integration_workspace_free(w2);
    
    if(z <= hbb/2.){
        top	= 0;
    }
    else{
        ulim	= atan(rout/(z-hbb/2.));
        blim	= atan(rin/(z-hbb/2.));        
        gsl_integration_workspace *w3   = gsl_integration_workspace_alloc(1000);
        double result3, error3;
        gsl_function F3;
        struct bbfnc_params F3params    = {gamv, z, hbb, tin, rin, nu};
        F3.function     = &bbfnc;
        F3.params       = &F3params;
        gsl_integration_qag(&F3, blim, ulim, 1e-1, 1e-1, 1000, 1, w3, &result3, &error3);
        top     = result3;
        gsl_integration_workspace_free(w3);
    }
    
    return nu*(inside+top)*2.*herg * nu*nu*nu /(cee*cee);
}

/**
 * bbdiskfnc for integration of bbdisk in bbintegral
 */
double bbdiskfnc(double nu, void *p){
    struct bbdiskfnc_params *params = (struct bbdiskfnc_params *)p;
    double rin      = (params -> rin);
    double rout     = (params -> rout);
    double z        = (params -> z);
    double hbb      = (params -> hbb);
    double gamv     = (params -> gamv);
    double tin      = (params -> tin);
    
    return bbdisk(nu, gamv, tin, rin, rout, z, hbb);
}

/**
 * bbintegrals calculates the integration over log(nubot)-log(nutop) of bbdisk
 */
void bbintegrals(double lnub,double lnut,double gamv,double tin,double rin,double rout,double z, double hbb,
double &bbint){
    gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(1000);
    double result1, error1;
    gsl_function F1;
    struct bbdiskfnc_params F1params = {rin, rout, z, hbb, gamv, tin};
    F1.function     = &bbdiskfnc;
    F1.params       = &F1params;
    gsl_integration_qag(&F1, lnub, lnut, 1e-1, 1e-1, 1000, 1, w1, &result1, &error1);
    bbint   = result1;
    gsl_integration_workspace_free(w1);
}


/**
 * bbearth gives the flux at Earth as a function of r(T) and freq at distance dist
 */
double bbearth(double lr, void *p){
    struct bbearth_params *params = (struct bbearth_params *)p;
    double tin      = (params -> tin);
    double rin      = (params -> rin);
    double frq      = (params -> frq);
    double inclin   = (params -> inclin);
    double dist     = (params -> dist);
    
    double r, temp, fac, bb;
    
    r	= exp(lr);
    temp	= tin*pow(rin/r, 0.75);
    fac	= herg*frq/(kboltz*temp);
    
    if(fac<1.e-3){
        bb	= 2.*herg*frq*frq*frq/(cee*cee*fac);
    }
    else{
        bb	= 2.*herg*frq*frq*frq/(cee*cee*(exp(fac)-1.));
    }
    
    return r*cos(inclin)*2.*pi*r*bb/(dist*dist);
}
/**
 * bbearthint is the integration of bbearth function
 */
void bbearthint(double blim,double ulim,double frq,double tin,double rin,double dist,double inclin, 
double & bbflx){
    gsl_integration_workspace *w   = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F1;
    struct bbearth_params F1params  = {tin, rin, frq, inclin, dist};
    F1.function     = &bbearth;
    F1.params       = &F1params;
    gsl_integration_qag(&F1, blim, ulim, 1e-1, 1e-1, 1000, 1, w, &result, &error);
    bbflx   = result;
    gsl_integration_workspace_free(w);
}
