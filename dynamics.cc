#include "agnjet.hh"

using namespace std;

//Defines zonal increments in the code to avoid resolution dependance

void jetgrid(int k,int infosw,int nz,double r_g,double velsw,double h0,double r,double zmin,double zcut,
double zmax,int &zone_zcut,double &delz,double &z,double zed[]){

	double zinc;										//increment in logarithmic grid 

	if (velsw <=1){		
		if(k==0){
			z = pow(10.,zmin)*(1.+pow(10.,-6));
			delz = h0;
			zone_zcut = zone_zcut + 1;
		}
		else if (z + 2.*r < zcut){
			z = z + delz;
			delz = 2.*r;
			zone_zcut = zone_zcut + 1;
		}
		else{
       	 	zinc = (log10(zmax) - log10(zcut))/(nz - zone_zcut);
			z = pow(10.,log10(zcut) + zinc*(k-zone_zcut)) * (1. + pow(10.,-6));
			delz = z - z*pow(10.,-zinc);
			if (infosw == 1) {
				cout << "Grid shifted; distance: " << z/r_g << " r_g; zone number " << k  << endl;
				if (z >= zcut) {
					cout << "Out of the Comptonization region; distance: " << z/r_g << endl;
				}	
			}
		}
	} else {
		if(k==0){
        	z = pow(10.,zmin)*(1.+pow(10.,-6));
            delz = h0;
            zone_zcut = zone_zcut + 1;
		}
        else if (z + 2.*r < zcut/80.){
        	z = z + delz;
			delz = 2.*r;
            zone_zcut = zone_zcut + 1;
            if (infosw ==1 ){
            	cout << "Zone number: " << k << endl;
            	if (z >= zcut) {
					cout << "Out of the Comptonization region; distance: " << z/r_g << endl;
				}
            }            
		}
        else{
        	zinc = (log10(zmax) - log10(zcut/80.))/(nz - zone_zcut);
            z = pow(10.,log10(zcut/80.) + zinc*(k-zone_zcut)) * (1. + pow(10.,-6));
            delz = z - z*pow(10.,-zinc);
			if (infosw == 1) {
				cout << "Grid shifted; zone number " << k  << endl;
				if (z >= zcut) {
					cout << "Out of the Comptonization region; distance: " << z/r_g << endl;
				}	
			}
		}
	}
	zed[k] = z;
}

/**
 * Main function of the jet properties
 * **************************************
 *
 * @param z		distance from origin along jet axis (same units as r0, h0
 * @param b0		magnetic field in nozzle
 * @param n0		density in nozzle
 * @param g0		electron Lorentz factor in nozzle
 * @param r0		width of nozzle
 * @param h0		height of nozzle
 * @param nl		number of lines of the parameters array
 * @param par       parameters array
 * @param flag_sol 	flag for selection the jet profile (0 = old fixed profile, 1 = first rel sol)
 *
 * @return gb	bulk velocity of jet (gamme*beta)
 * @return dgb	analytical derivative of bulk velocity
 * @return b	magnetic field at z
 * @return n	density at z
 * @return g	electron Lorentz factor at z
 * @return r	width at z
 * @return k_equip  equipartition parameter @ z
 *
 * NOTE: mj is the Mach number but relative to soundspeed in nozzle! Local Mach number will be slightly higher
 * due to cooling
 */
void jetpars(double velsw,double zeta,double z0,double z,double zsh,double r_g,double b0,double n0,double g0, 
double r0,double h0,gsl_spline *spline,gsl_interp_accel *acc,double &gb,double &b,double &n,double &g,
double &r){
    
	double y	= 0;
	double mj, gbs0, Gammas, betas0;    
	double x;
        
	Gammas	= 4./3.;
    betas0	= sqrt(Gammas*(Gammas-1.)/(Gammas+1.));
    gbs0	= betas0/sqrt((1.-betas0)*(1.+betas0));
    x	= log10((max(z-h0,0.)+r0)/r0);
        
    if(z < h0){
    	y	= gbs0;
	} else {
    	x = pow(10,x);
        y = gsl_spline_eval(spline, x, acc);
	}
        
	gb	= y;
    mj	= gb/gbs0;        
    r	= r0+max(z-h0, 0.)/mj;        
    n	= zeta*n0*(r0/r)*(r0/r)/mj;

    if(z < h0){
    	b = b0*(r0/r)/pow(mj/zeta,0.5);
        g = g0;
	} else {
		if (velsw== 0) {
    		b = b0*(r0/r)/pow(mj/zeta, 0.5+1./6.);
        	g = g0/pow(mj, 1./3.)/pow(r/r0, 2./3.);
		} else {
    		b = b0*sqrt(zeta)*(r0/r)/pow(mj, 0.5+1./6.);
        	g = g0/pow(mj, 1./3.);
		}
	} 
}

/*
 * VELOCITY PROFILE FUNCTION
 * -------------------------
 * @param z		distance from origin along jet axis (same units as r0, h0
 * @param r0		width of nozzle
 * @param h0		height of nozzle
 * @param spline	first parameter for interpolation
 * @param acc 		second parameter for interpolation (see GSL manual)
 *
 * @return vel		velocity profile
 * @return dvel		derivative of the velocity profile
 *
 * NOTE: this is a slightly different version of the old mini routine that was present
 * in jpars, but expanded to include the derivative. The derivative is needed in the Compton routine
 */

void vel_prof(double z, double mbh, double zmin, double visco, double r0, double h0, gsl_spline *spline, 
gsl_interp_accel *acc, double &vel, double &dvel){
    
    double x, betas0, Gammas, zff, gbin, bff, gbs0;    
    
    Gammas	= 4./3.;
    betas0	= sqrt((Gammas-1.)/(Gammas+1.));
    gbs0	= betas0/sqrt((1.-betas0)*(1.+betas0));    
    x	= log10((max(z-h0,0.)+r0)/r0);    
    zff	= pow(10, zmin);
    
    if(z < h0){
        if(visco <= 0.){
            vel	= gbs0;
            dvel	= 0;
        }
        else{
            bff	= sqrt(2.*6.67e-8*2.e33*mbh/r0)/cee*visco;
            gbin	= bff/sqrt(1.-bff*bff);
            vel	= gbin+(gbs0-gbin)/(h0-zff)*(z-zff);
            dvel	= -gbs0/(h0-zff);
        }
    }
    else{
        vel       = gsl_spline_eval(spline, x, acc);
        dvel	  = gsl_spline_eval_deriv(spline, x, acc);
    }
}

void bljetpars(double mxsw,double velsw,double z,double zacc,double r_g,double eta,double n0,double g0,
double r0,double h0,double endnsmj,double pspec,double cnorm,double emin,double emax,double ebreak, 
gsl_spline *spline,gsl_interp_accel *acc,double &gb,double &b,double &n,double &g,double &r,double sigsh){

	double y	= 0;
	double mj, gam0, gbs0, Gammas, betas0, theta, tshock, nshock, bshock, gam;
             
	Gammas	= 4./3.;
	betas0	= sqrt(Gammas*(Gammas-1.)/(Gammas+1.));
	gam0 = 1./sqrt((1.-betas0)*(1.+betas0));
	gbs0	= betas0/sqrt((1.-betas0)*(1.+betas0));
               
	if (z < h0){
		y = gbs0;
	} else if (z < zacc){
		y = gsl_spline_eval(spline, z/r_g, acc);
	} else {
		y = sqrt(velsw*velsw-1.);
	}
	
	gb = y;
	mj = gb/gbs0;
	gam = sqrt(y*y+1.);
	
	theta = 0.15/gam;
	r = r0+max(z-h0,0.)*tan(theta);
	tshock = 0.15/velsw;
	nshock = eta*n0*pow(h0/zacc,2.)*(gbs0/sqrt(velsw*velsw-1.))*pow(sin(0.15/gam0)/sin(tshock),2.);
	n = eta*n0*pow(r0/r,2.)/mj;
	
	if(z < h0){
		b_prof(mxsw, velsw, eta, n, endnsmj, pspec, cnorm, emin, emax, ebreak, b, gam, sigsh);
		g = g0;
	} else if (z < zacc){
		g = g0;;
		b_prof(mxsw, velsw, eta, n, endnsmj, pspec, cnorm, emin, emax, ebreak, b, gam, sigsh);
	} else {
		g = g0;
		b_prof(mxsw, velsw, eta, nshock, endnsmj, pspec, cnorm, emin, emax, ebreak, bshock, velsw, sigsh);
		b = bshock*(zacc/z);
	}
		
}

/* BL Lac magnetic field function
 * ------------------------------
 * This function computes the magnetic field as it is dissipated to reach high velocities (velsw > 1)
 * @param mxsw 		type of particle distribution
 * @param z		distance along the jet in r_g
 * @param zacc		length in r_g of the jet acceleration zone
 * @param r_g		gravitational radius of the BH
 * @param gfin		maximum Lorentz factor reached
 * @param n 		electron number density along the jet
 * @param endnsmj	electron M-J energy density
 * @param pspec		non thermal power-law index
 * @param cnorm		non thermal power-law normalization
 * @param gmin 		non thermal power-law minimum Lorentz factor
 * @param gmax 		non thermal power-law maximum Lorentz factor
 * @param spline	input Lorentz factor as a function of distance
 * @param acc		second parameter of interp
 *
 * @return field	magnetic field in the acceleration zone  
 */

void b_prof(double mxsw,double gfin,double eta,double n,double endnsmj,double pspec,double cnorm,double emin,
double emax,double ebreak,double &field,double g,double sigsh){

	double betas0, g0, sig0, sig, w;
	double Gammas = 4./3.;
	betas0 = sqrt(Gammas*(Gammas-1.)/(Gammas+1.));
   	g0 = 1./sqrt((1.-betas0)*(1.+betas0));

	//Only accounts for thermal distribution
	w = Gammas*eta*n*endnsmj;
	sig0 = (1.+sigsh)*gfin/g0 - 1.;
	sig = (g0/g)*(1.+sig0)-1.;
	field = sqrt(sig*4.*pi*(n*pmgm*cee*cee+w));
}

//parameters of each jet section
void zonepars(int njet,int k,int nz,double z,double r,double delz,double rvel,double tbb,double inclin,
double gamax0,double gamax,int &nw,double &area,double &vol,double &gamv2,double &gamv,double &beta,
double &tbbeff,double &gshift,double dopfac[]){

	int l;
	
 	//Area and volume of jet's segment
	area = pi*r*r;
	vol = delz*area;	

	//Express jet velocity in the zone
	gamv2 = 1.+rvel*rvel;
	gamv = sqrt(gamv2);
	beta = sqrt(gamv2-1.)/gamv;

	//Setting Doppler factor in jet segment
    for(l=0; l< njet; l++){
		dopfac[l*nz+k]	= 1./(gamv*(1. - beta*cos(inclin)*pow(-1.,l)));
	}
	
	//Calculate relevant effective temperature (see R&L 153 or 156, iv) for secondary photon fields
    tbbeff = tbb*gamv;
	//tbbeff = tbb;
	
    //Express gshift and nw, used in shifting down energy distribution
    gshift  = gamax/gamax0;
    nw      = 0;
}
