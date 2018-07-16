#include "agnjet.hh"

using namespace std;

/* Initialize disk quantities */

void ext_init(int infosw,int compsw,double jetrat,double compar1,double compar2,double compar3,double rin,
double rout,double eddlum,double &tin,int &disksw,double &thrlum,double &eddrat,double &tbb,double &normbb,
double &rbb, double &hbb){

	double outfac;
	outfac = rout/rin;
	
    if(outfac < 1){
    	disksw  = 0;
        if (infosw ==1) {
           	cout << "rin > rout, no disk component!" << endl;
    	}
	}

	//Setting disk temperature/checking if it's super Eddington, normalizations of disk and second blackbody
	if (tin > 0) { 
    	thrlum  = 4*pi*sbconst*pow(tin,4)*pow(rin,2)*(1-1/outfac);
        eddrat = thrlum/eddlum;
        hbb = rin*eddrat;
        if(eddrat > 1){
        	cout << "Disk flux super-Eddington, please readjust rin/tin. eddrat larger than 1: " << eddrat <<
        	endl;                
    	}       
	}
	else {
    	eddrat = - tin * kboltz_kev;
        thrlum = eddrat*eddlum;
		tin = pow(thrlum/(4*pi*sbconst*pow(rin,2)*(1-1/outfac)),0.25);	
		hbb = rin*eddrat;
	}	
	
	//setting external photon fields;
	/*
	 compsw == 1: homogeneus external black body
	 			  compar1 is the bb temperature in the observer frame
	 			  compar2 is the (desired or estimated) bb energy density in the observer frame
	 			  compar3 is the radius of the bb
	 			  the normalization/luminosity are calculated from the desired urad (ie what's needed for IC)
	*/
	if (compsw == 1) {
		tbb = compar1; //bbody temperature (from kev to kelvin?)
		rbb = compar3; 	//bbody radius (U_bb = normbb/4*pi*rbb^2)
		normbb = 4*pow(pi,2.)*pow(rbb,2.)*cee*compar2/(sbconst*pow(tbb,4.)); //bbody normalization    	
    }
};

void jet_init(int zfrac,int sizegb,double mxsw,double velsw,double jetrat,double r_g,double r0,double hratio,
double zacc,double zmax,double beta_pl,double &equip,double &h0,double &zsh,double &zcut,double &gad4_3,
double &betas0,double &gam0,double &rvel0,double &vel0,double &zmin,double gbx[],double gby[],
double gbx_vel2[],double gby_vel2[]){

	double step;
	
	equip = 1./beta_pl;
	h0 = hratio*r0;	
	
	if (mxsw != 1){ //if non thermal particles injected then set large zsh automatically
			zsh = zmax+r0+h0;
	}
	
	//Setting distance in the jet above which comptonisation is neglected for computational speed
	zcut    = zfrac*r0;
	if (velsw > 1 && r0 < 10.*r_g){
		zcut = 5.*zcut;
	}
	
	//Setting initial jet velocity
	gad4_3  = 4./3.;
	betas0  = sqrt(gad4_3*(gad4_3 - 1.)/(gad4_3 + 1.)); 
	gam0    = 1./sqrt(1.-(betas0*betas0));
	rvel0   = gam0*betas0;
	vel0    = cee*rvel0;
	
	//Set up the updated (2017) grid; start with delz = 2r before going to a log grid
	zmin = int(log10(r0))*1.0; //zmin should have no impact on the SED because we neglect gravity
	
	if (velsw == 0) {
		cout << "The model is running with adiabatic agnjet profile, careful!" << endl;
	}
	
	if (velsw > 1){
		step = (log10(zmax/r_g) + 1. -log10(pow(10.,zmin)/r_g))/(sizegb-3.);				
		for (int i = 0; i<sizegb-1; i++){
			gbx[i] = (log10(pow(10.,zmin)/r_g)+i*step);
			gbx[i] = pow(10.,gbx[i]);
		if (gbx[i]*r_g < h0) {
			gby[i] = gam0;
		} else if (gbx[i]*r_g < zacc){
			gby[i] = gam0 + ((velsw-gam0)/((pow(zacc,0.5)-pow(h0,0.5))))*(pow(gbx[i]*r_g,0.5)-pow(h0,0.5));
		} else {
			gby[i] = velsw;				
		}
		gby[i] = sqrt(pow(gby[i],2.)-1.);
		}
	} else {
	    for (int i = 0; i<sizegb-1; i++) {
			gbx[i] = gbx_vel2[i];
		    gby[i] = gby_vel2[i];
		}	
	}  	
};

void particles_init(double mxsw,double velsw,double jetrat,double r0,double eltemp,double pspec,double gamfac,
double sigsh,double plfrac,double gad4_3,double vel0,double gam0,double &thmfrac,double &emin,double &gamin,
double &ebreak,double &emax,double &uplim,double &endncom,double &endnsmj,double &betat,double &nprot0,
double &eta,double &equip,double &sig0,double &beta_pl){

	double numcorr;										//normalization factor in number density calculation
	double avgen;										//average energy of injected particles
										
    thmfrac = 1.-plfrac;     
    emin    = 2.23*kboltz*eltemp;
    gamin   = emin/emerg;
    emax    = gamfac*gamin*emerg;
    ebreak = emax + 1.;
    uplim   = pow(10, -54.32 + 9.416*log10(eltemp) - 0.376*log10(eltemp)*log10(eltemp));
    endncom = enden(eltemp, uplim);
    endnsmj = endncom*emerg;
    betat   = emerg/(kboltz*eltemp);
    
	if (velsw <= 1){
        	nprot0  = (jetrat/4.)/(vel0*pmgm*cee*cee*pi*r0*r0);
	} else {
		if (mxsw == 1){
			avgen = endncom;
		} else if (mxsw == 0){
			if (pspec >= 1.99 && pspec <= 2.01){
				avgen = log(emax/emin)*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,1.-pspec))*emerg);
			} else {
				avgen = (pow(emax,2.-pspec)-pow(emin,2.-pspec))*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,
				1.-pspec)*(2.-pspec))*emerg);
			}
		} else {
			if (pspec >= 1.99 && pspec <= 2.01){
				avgen = mxsw*endncom+(1.-mxsw)*log(emax/emin)*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,
				1.-pspec))*emerg);
			} else {
				avgen = mxsw*endncom+(1.-mxsw)*(pow(emax,2.-pspec)-pow(emin,2.-pspec))*(1.-pspec)/
				((pow(emax,1.-pspec)-pow(emin,1.-pspec)*(2.-pspec))*emerg);
			}
		}
		sig0 = ((1.+sigsh)*velsw)/gam0-1.;
		if (beta_pl == 0){
			eta = 1.;
			equip = (sig0/2.)*(gad4_3+(pmgm)/(avgen*emgm));
		} else {
			eta = (sig0*pmgm)/(avgen*emgm*(2.*equip-sig0*gad4_3));
		}
		numcorr = 1.+eta*(1.+equip)*(avgen*emgm)/(pmgm);
		nprot0  = jetrat/(2.*numcorr*vel0*pi*r0*r0*pmgm*cee*cee); 
	}
	
	if (velsw > 1 && eta < 1) {
		equip = equip/3.;
		beta_pl = 1./equip;
		eta = (sig0*pmgm)/(avgen*emgm*(2.*equip-sig0*gad4_3));
		numcorr = 1.+eta*(1.+equip)*(avgen*emgm)/(pmgm);
		nprot0  = jetrat/(2.*numcorr*vel0*pi*r0*r0*pmgm*cee*cee); 
		cout << "Changed equip to: " << equip << " and plasma beta to: " << beta_pl << endl;
		cout << "Check plasma beta value! Pair content now is: " << eta << endl;
	}
	if (nprot0 < 0){
		cout << "Proton number density: " << nprot0 << endl;
		cout << "Too many pairs, the jet has a negative proton number density!" << endl;		
	}
};



/* Rough estimate of average energy density */

double average_ene(int flaglog,int N,double ndensity[],double energy[],double eldens[]){

    double sum1,sum2,tmp1,tmp2;
    
    double *elexdens	= new double[N]();
    
    sum1	= 0;
    sum2	= 0;
    tmp1    = 0;
    tmp2    = 0;
    for(int i=0; i<N; i++){
        if (flaglog == 1) {
            elexdens[i]	= pow(10,energy[i]+ndensity[i]);
            energy[i]   = pow(10,energy[i]);
            ndensity[i]   = pow(10,ndensity[i]);
        } else {
            elexdens[i]	= eldens[i];
        }
    }
    for(int i=0; i<N-1; i++){
        
        tmp1	= (energy[i+1]-energy[i])*0.5*(elexdens[i]+elexdens[i+1]);
        tmp2	= (energy[i+1]-energy[i])*0.5*(ndensity[i]+ndensity[i+1]);
        sum1	+=tmp1;
        sum2	+=tmp2;
    }
    for(int i=0; i<N; i++){
        if (flaglog == 1) {
            elexdens[i]	= log10(elexdens[i]);
            energy[i]   = log10(energy[i]);
            ndensity[i]   = log10(ndensity[i]);
        } else {
            elexdens[i]	= eldens[i];
        }
    }
//    cout << i << " " << tmp1 << " " << tmp2 << " " << sum1 << " " << sum2 << endl;
    return sum1/sum2;
}

/* Energy arrays */  

void ene_arrays(bool &isSSC,int ne,int nsyn,int ncom,int &zfrac,double mbh,double inclin,int &njet,
double &snumin,double &snumax,double &snuinc,double &cnumin,double &cnumax,double &cnuinc,double ear[],
double nutot[],double nubb[],double nurad[],double energ[],double ephot[],double ephxr[]){
    
    int i;
    
    double *ebin	= new double[NEBIN]();
    
    snumin = log10(0.5*ear[0]/hkev);				//lower boundary Synchrotron grid
    snumax = log10(ear[ne-1]/hkev);					//upper boundary Synchrotron grid
   								
    
    //If aligned agn, extend Compton calculation and grid, disable multiple Compton, disable second jet.
    if (mbh > 1e4 || inclin*(180./pi) <= 20.){
    	cnumin = 14;								//lower boundary Compton grid
    	cnumax = 27;								//upper boundary Compton grid
    	cnuinc = (cnumax-cnumin)/ncom;
		zfrac = 100*zfrac;
		njet = 1;
		isSSC = true;
	} else {		
	    cnumin = 12;								//lower boundary Compton grid
	    cnumax = 23;								//upper boundary Compton grid
	    cnuinc = (cnumax-cnumin)/ncom;		
	}

    snuinc = (snumax-snumin)/nsyn;					//increment of Synchrotron grid
    cnuinc = (cnumax-cnumin)/ncom;					//increment of Compton grid
   
    for(i=0; i<(ne-1); i++){
        ebin[i] = ear[i] + (ear[i+1]-ear[i])/2;
        nutot[i]= log10(ebin[i]/hkev);
    }
    
    //Define radiation grids

    //Synchrotron
    for(i=0; i<nsyn; i++){
        nubb[i] = snumin + (i+0.5)*snuinc;
        nurad[i]= pow(10., nubb[i]);
        energ[i]= herg*nurad[i];
        ephot[i]   = log10(energ[i]);
    }
    //Compton
    for(i=0; i<ncom; i++){
        ephxr[i]= pow(10, cnumin+(i+0.5)*cnuinc);
    }
    delete[] ebin;
}

/**
 * Modified Bessel function
 *
 */
double k0_fnc(double x){
    double k0, y, t, i0, t2, y2, y1;
    
    /**
     * k0 function from A&S good to 2.e-7
     *
     */
    y	= x/2.;
    y2	= y*y;
    y1	= 1./y;
    
    if(y < 1.){
        t	= x/3.75;
        t2	= t*t;
        i0	= 1. + t2*(3.5156229+t2*(3.0899424+t2*(1.2067492*pow(t,6)+t2*(0.2659732*pow(t,8)+t2*(0.0360768+
        t2*0.0045813)))));
        k0	= log(2./x)*i0 - 0.5772156649 + y2*(0.42278420+y2*(0.23069756+y2*(0.03488590+y2*(0.00262698+y2*
        (0.00010750+y2*0.00000740)))));
    }
    else{
        k0	= exp(-x)/sqrt(x) * (1.25331414+y1*(-0.07832358+y1*(0.02189568+y1*(-0.01062446+y1*(0.00587872+y1*
        (-0.00251540+y1*0.00053208))))));
    }
    return k0;
}


double k1_fnc(double x){
    double k1, y, t, i1, t2, y2, y1;
    
    /**
     * k1 function from A&S good to 2.e-7
     *
     */
    y	= x/2;
    y2	= y*y;
    y1	= 1./y;
    
    if(y<1){
        t	= x/3.75;
        t2	= t*t;
        i1=x*(0.5e0+t2*(0.87890594e0 + t2*(0.51498869e0 + t2*(0.15084934e0 + t2*(0.02658733e0 + t2*
        (0.00301532e0 + t2*0.00032411e0))))));
        k1 = 1.e0/x*(x*log(y)*i1+1.e0 + y2*(0.15443144e0-y2*(0.67278579e0 - y2*(0.18156897e0-y2*(0.01919402e0
        - y2*(0.00110404e0-y2*0.00004686e0))))));
    }
    else{
        k1 = exp(-x)/sqrt(x)*(1.25331414e0+y1*(0.23498618e0 - y1*(0.03655620e0+y1*(0.01504268e0 - y1*
        (0.00780353e0+y1*(0.00325614e0 - y1*0.00068245e0))))));
    }
    return k1;
}

double k2_fnc(double x){
    double k2, k0, k1;
    
    /**
     * k2 function from A&S is just k2 = 2.*k1/x+k0
     *
     */
    k0	= k0_fnc(x);
    k1	= k1_fnc(x);
    
    k2	= 2.*k1/x+k0;
    
    return k2;
}


double k3_fnc(double x){
    double k3, k2, k1;
    
    /**
     * k3 function from A&S is just k3 = 4.*k2/x+k1
     *
     */
    k1	= k1_fnc(x);
    k2	= k2_fnc(x);
    
    k3	= 4.*k2/x+k1;
    
    return k3;
}

/**
 * Interpolate from grid sized with ne bins to one of newne bins
 * Adapted from S-lang version, for greater speed (Mike Noble)
 *
 */
void bhinterp(double *ear, double *energ, double *phot, double *photar, int ne, int newne){
    int i, iplus1, j, jplus1;
    double emid, phflux;
    
    j = 0; 		//j = 1 or j = 0??
    
    for(i=0; i<newne; i++){
        
        // Middle of bin
        iplus1	= i+1;
        emid	= (ear[i] + ear[iplus1])/2.;
        
        // Linear search if we don't bracket yet
        if(j == -1){
            j = 1;
        }
        while(j <= ne && energ[j] < emid){
            j++;
        }
        
        jplus1	= j;
        j	= j - 1;
        
        if(j < 1 || j > ne){
            photar[i]	= 0.;
        }
        else{
            // ph/cm^2/s/keV
            phflux	= phot[j] + (phot[jplus1] - phot[j])*(emid-energ[j])/(energ[jplus1]-energ[j]);
            photar[i]= phflux * (ear[iplus1] - ear[i]);
        }
    }
}
