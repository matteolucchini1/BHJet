// =====================================================================================
// 
//       Filename:  agnjet_merged.cc
// 
//    Description:  Main file of the leptonic model of a X-ray binary jet.
//                  Work based on the agnjet.f code from S.B. Markoff.
// 
//        Version:  1.0
//        Created:  09/15/2012 13:36:59
//       Revision:  none
//       Compiler:  
// 
//         Author:  Samia Drappeau (sd), drappeau.samia@gmail.com
//  Maintained by:  Matteo Lucchini, m.lucchini@uva.nl
//        Company:  API - UvA/
//
//         
// 
// =====================================================================================

#include "agnjet_merged.hh"

using namespace std;

/**
 * xrbjet (ear, ne, params, phot_spect)
 * Main function of the code
 * ************************************
 * @param ear           energy array
 * @param ne            size of ear
 * @param params        model's parameters
 *
 * @return phot_spect   photon's spectra
 * @return photend      photon's energy array
 *
 */
void xrbjet(double* ear, int ne, double *param, double *phot_spect, double *photeng){
        /**
         * Declare for-loop, sum and tmp variables
         *
         */
	int i, j, k, l, m;
	double sum1, sum2, tmp1, tmp2;

        /**
         * Flags
         * ****************************************************
         * isShock              Flag for shock
         * isInitialized        Flag for initialization of jet
         * isBreaknjetnsyn      Flag for synchrotron for loop <!-- TODO better name -- >
         * isVerbose            Flag for verbose debug mode
         *
         */
	bool isShock         = false;
	bool isHeated	     = false;
	bool isInitialized   = false;
	bool isBreaknjetnsyn = false;
	bool isVerbose       = false;
    
        //---------------------- PROTOTYPE CHANGE --------------------
        //        if 1, Peter's solution is used
	int flag_sol             = 0;
        //------------------------------------------------------------
    

        /**
         * Define Synchrotron - Compton - Gamma-ray grid [Hz]
         * ************************************************
         * nsyn         number of bins in Synchrotron grid
         * snumin       lower boundary Synchrotron grid
         * snumax       upper boundary Synchrotron grid
         * snuinc       incrementation Synchrotron grid
         * ncom         number of bins in Compton grid
         * cnumin       lower boundary Compton grid
         * cnumax       upper boundary Compton grid
         * cnuinc       incrementation Compton grid
         *
         */
	int nsyn    	= 100;
	int ncom    	= 60;
        double snumin  	= log10(0.5*ear[0]/hkev);
        double snumax  	= log10(ear[ne-1]/hkev);
        double cnumin  	= 16;
        double cnumax  	= 26;

        double snuinc  	= (snumax-snumin)/nsyn;
        double cnuinc  	= (cnumax-cnumin)/ncom;
        int Niter = 500;

        /**
         * Parameters that used to be changeable but are now fixed
         * ********************************************************
         * zfrac        multiplier to set Comptonization region; zcut = zfrac*r0
         * bbsw         switch to get black-body radiation
         * disksw       switch to get disk radiation
         * zinc         size of segment with logarithmic grid
         * nzdum        size of jet's segment array
         * nz           number of zones in jet, is always set by zmax and zinc
         * outfac       factor so that rin !> rout
         * thrlum       thermal luminosity
         * eddrat       Eddington ratio 
         * thmfrac      1-plfrac (see jet's parameter plfrac down)
         * zcut         distance in the jet above which no more irradiation from disk <!-- TODO -->
         * njet         number of jets produced
         * zmin		base of the nozzle
         */
        int nz;
        double outfac, thrlum, eddrat, thmfrac, zcut, zmin;
        int zfrac       = 100;
        int bbsw    	= 1;
        int disksw  	= 1;
        int nzdum   	= 100;   
        int njet    	= 2;

        /**
         * Variables related to the speed of sound in nozzle <!-- TODO -->
         * **************************************************
         * gad4_3       gamma ideal gas <!-- TODO -->
         * betas0       beta speed of sound in nozzle <!-- TODO -->
         *              sqrt((gad4_3 - 1.)/(gad4_3 + 1.))
         * gam0         1./sqrt(1.-(betas0*betas0))
         * rvel0        gam0*betas0
         * vel0         cee*rvel0
         * gamax0       [set to unity]
         * gbs0         <!-- TODO insert definition -->
         *
         */
        double gad4_3, betas0, gam0, rvel0, vel0, gamax0, gbs0;

        /**
         * Variables related to particles and energies
         * ************************************************
         * emin         minimum energy of eletron
         * emax         maximum energy of electron
	 * ebreak       energy of synchrotron cooling break
         * gamin        emin/emerg
         * blp          bottom limit of the mj integration in pc
         * ulp          upper limit of the mj integration in pc
         * endnsmj      Maxwell-J\¨uttner energy density
         * nprot0       initial number of protons <!-- TODO add units & expression -->
         * ntot0        initial number of electrons
         * b_en         magnetic energy in nozzle <!-- TODO check if in nozzle -->
         * cnorm        normalization factor of the power-law distribution <!-- TODO check whether correct or not -->
         * b0           magnetic field magnitude in nozzle
         * endens       electron's energy density
         * betat        emerg/(kboltz*eltemp)
         * enorm        ntot0*betat/(emerg*k2_fnc(betat));
         * elmin        electron minimum energy in thermal distribution
         * elmax        electron maximum energy in thermal distribution
         * einc         electron energy incrementation in thermal distribution
         * nelec        electron array size in one jet's segment
         * game         electron's gamma (rden/emerg)
         * bete         linked to game <!-- TODO better definition -->
         * reled        enorm*game*game*bete*exp(-game*betat)
         * elenmn       <!-- TODO insert definition -->
         * elenmx       <!-- TODO insert definition -->
	 * zacc		distance up until which the jet dissipates magnetic field and accelerates, if velsw > 1
         * avgen	average energy of input particle distribution; necessary for velsw > 1
	 * numcorr 	correction to the initial number of particles if velsw > 1
         * sig0		magnetization at the base of the jet for magnetically accelerated jet
	 * eta		number of pairs per proton when velsw > 1
         */
        double emin, emax, ebreak, gamin, uplim, endncom, endnsmj, nprot0, avgen, numcorr, sig0, eta;
        double ntot0, b_en, b0, cnorm;
        double endens, betat, enorm, elmin, elmax, einc, game, bete, reled;
        int nelec   = 200;
        double elenmn, elenmx;
        int nenpsw =1; // either 1 or 0, 1 if using correct normalisation in param file, 0 to calculate new normalisation given old param file - need this for cooling model
        
        /**
         * Variables related to jet's segment
         * *************************************
         * z            position of the segment
         * delz         height of segment
         * rvel         velocity in segment
         * ntot         electron number density in segment
         * gamax        <!-- TODO -->
         * r            radius of segment
         * area         area of the jet at z
         * vol          volume of segment
         * gamv2        <!-- TODO -->
         * beta         <!-- TODO -->
         * gshift       gamax/gamax0
         * nw           keep track of lost particle
         *              when shifting down energy distribution
         * renorm       renormalize the energy distribibution after shift down
         * totlos       total particle lost when shifting down
         * tmlos        temporary total particle lost
         * oldnum       store previous total particle number
         * gshock       gshift factor at the shock
         * ub           magnetic energy in jet's segment <!-- TODO check -->
         * ucom         rough Compton energy density
         * bbdil        <!-- TODO insert definition -->
         * nubot        <!-- TODO insert definition -->
         * nutop        <!-- TODO insert definition -->
         * tst          <!-- TODO insert definition -->
         * bbint        <!-- TODO insert definition -->
         * uphdil       <!-- TODO insert definition -->
         * uphdil2      <!-- TODO insert definition -->
         * accon        constant describing acceleration rate of shocked electrons
         * syncon       synchrotron cooling rate constant
         * comcon       compton cooling rate constant
         * escom        escape rate constant (travel time of electrons in jet)
         * tdyn         dynamical time of electrons in a jet segment (= delz/(beta*c))
         * qutrmb       escom/(syncon+comcon) - 'b' in quadratic solution to find emax
         * qutrmc       accon/(syncon+comcon) - 'c' in quadratic solution to find emax
         * gemax        gamma of maximum electron energy (gemax = emax/emerg)
         * delebreak    Change in break energy between two zones in the jet
         * gebreak      gamma of synchrotron cooling break energy (gebreak = ebreak/emerg) - ** ADDED JUNE 2015 **
         * bemax        maximum electron velocity of power law distribution (calculated from gemax)
         * mjteff       adiabatically cooled MJ temperature of electron distribution
         * mjpeak       peak of MJ distribution - used to calculate min energy of power law distribution
         * edtrm        electron distribution defined at shock
         * pltrm        power law term at shock, per iteration of energy array
         * shocknd      number density at shock
         * shockr       radius at shock
         * shockgb      rvelocity at shock
         * avelec       <!-- TODO insert definition -->
         * doppler      <!-- TODO insert definition -->
         * direct       <!-- TODO insert definition -->
         * blim         <!-- TODO insert definition -->
         * ulim         <!-- TODO insert definition -->
         * com          <!-- TODO insert definition -->
         * eph          <!-- TODO insert definition -->
         * ephmin       <!-- TODO insert definition -->
         * ephmax       <!-- TODO insert definition -->
         * absd         <!-- TODO insert definition -->
         * phoden       <!-- TODO insert definition -->
         * toten        <!-- TODO insert definition -->
         * bbnumax      <!-- TODO insert definition -->
         * photmx       <!-- TODO insert definition -->
         * peaksw       <!-- TODO insert definition -->
         * bbcon        <!-- TODO insert definition -->
         * hbb          <!-- TODO insert definition -->
         * bbirad       <!-- TODO insert definition -->
         * bbnrm        <!-- TODO insert definition -->
         * bbtrm        <!-- TODO insert definition -->
         * bfield       magnetic field magnitude
         * gamv         <!-- TODO insert definition -->
	 * t		number of zones before logarithmic grid, in which delz = r
         *
         */
        double z, delz, area, vol;
        double rvel     = 0;
        double drvel = 0;
        double gamax    = 0;
        double r        = 0;
        double ntot     = 0;
        double gamv2, beta, gshift, renorm, oldnum, gshock, ub, ucom;
        int nw;
        double totlos   = 0.;
        double tmlos    = 0.;
        double bbdil, nubot, nutop, tst, bbint, uphdil, uphdil2;
        double accon, syncon, comcon, escom, tdyn, qutrmb, qutrmc;
        double gemax, bemax, gebreak, delebreak;
        double mjteff, mjpeak, edtrm, pltrm;
        double shocknd = 0.;
        double shockr  = 0.;
        double shockgb = 0.;
        double avelec;
        double doppler  = 0.;
        double direct   = 0.;
        double blim, ulim, com, eph, ephmin, ephmax;
        double absd, phoden, toten, bbnumax;
        double photmx, peaksw;
        double bbcon    = 0.;
        double hbb;
        double bbirad, bbnrm, bbtrm, theff, tbbeff;
       	double reff, reff2;        
        double bfield, gamv;
	int t = 0;
	/**
         * Variables related to fluxes
         * ****************************
         * sflx         synchrotron flux
         * cflx         compton flux
         * bbflx        blackbody flux
         * annuarea     <!-- TODO insert definition -->
         * frq          frequency (used to be erdpar.frq)
         * inclin       inclination
         * esum         <!-- TODO insert definition -->
         * asum         <!-- TODO insert definition -->
         * avphot       <!-- TODO insert definition -->
         *
         */
        double sflx, cflx, bbflx, annuarea;
        double frq, inclin;
        double esum, asum, avphot;

        /**
         * Define all pointers used in code
         * *********************************
         * ebin         energy bins
         * nutot        frequency bins
         * nubb         frequency bins for
         *              black-body <!-- TODO -->
         * nurad        pow(10., nubb)
         * energ        herg*nurad
         * ephot        log10(energ)
         * ephxr        x-ray photon energy array
         *              (for Compton)
         * rdlgen       log10(electorn energy) <!-- TODO find a better variable name -->
         * rden         electron energy
         * rdend        electron energy density: reled*rden <!-- TODO used nowhere in the code... -->
         * rdedn        electron energy distribution
         * dopfac       Doppler factor in each jet's segment
         * lelec        <!-- TODO insert definition -->
         * elen         electron energy defined for use in z loop
         * ledens       log10(electron number density)
         * etemp        temporary electron energy array <!-- TODO change name -->
         * dtemp        temporary electron energy density array <!-- TODO change name -->
         * eled         <!-- TODO insert definition -->
         * phofrq       frequency of photon <!-- TODO better definition -->
         * nusyn        frequency of the synchrotron photon <!-- TODO better definition -->
         * nucom        <!-- TODO insert definition -->
         * comspc       <!-- TODO insert definition -->
         * phoint       <!-- TODO insert definition -->
         * synabs       <!-- TODO insert definition -->
         * shelen       elen at shock
         * shden        <!-- TODO insert definition -->
         * thmshock     thermal component array at the shock
         * thmcomp      thermal component array in each jet segment
         * plshock      power law component array at the shock
         * plcomp       power law component array in each jet segment
         * tdyntot      Total dynamical time of electron population ( = sum(delz/(beta*c)))
         * tsyntot      Total characteristic synchrotron cooling time, summed up over each segment
         * tstrm        <!-- TODO insert definition -->
         * drtrm        <!-- TODO insert definition -->
	 * phodis       <!-- TODO insert definition -->
         * nphot        <!-- TODO insert definition -->
         * jetu         <!-- TODO insert definition -->
         * disku        <!-- TODO insert definition -->
         * snu          synchrotron frequency array
         * sdump        synchrotron emission array
         * cnu          compton frequency array
         * cdump        compton emission array
         * complot      <!-- TODO insert definition -->
         * presyn       <!-- TODO insert definition -->
         * postsyn      <!-- TODO insert definition -->
         * fplot        <!-- TODO insert definition -->
         * bbplot       <!-- TODO insert definition -->
         * total        <!-- TODO insert definition -->
	 * zed		distance array along the jet to be passed to radiation routines
         *
         */
        double *ebin	= new double[NEBIN]();
        double *nutot	= new double[NEBIN]();
        double *nubb	= new double[nsyn]();
        double *nurad	= new double[nsyn]();
        double *energ	= new double[NEBIN]();
        double *ephot	= new double[nsyn]();
        double *ephxr	= new double[ncom]();
        double *rdlgen	= new double[nelec]();
        double *rden	= new double[nelec]();
        double *rdend	= new double[nelec]();
        double *rdedn	= new double[nelec]();
        double *dopfac	= new double[njet*nzdum]();
        double *lelec	= new double[nelec]();
        double *elen	= new double[nelec]();
        double *ledens	= new double[nelec]();
        double *etemp	= new double[nelec]();
        double *dtemp	= new double[nelec]();
        double *eled	= new double[nelec]();
        double *phofrq	= new double[nsyn]();
        double *nusyn	= new double[njet*nzdum*nsyn]();
        double *nucom	= new double[njet*nzdum*ncom]();
        double *comspc	= new double[njet*nzdum*ncom]();
        double *comspc_src = new double[ncom]();
        double *nucom_src = new double[ncom]();
        double *comi = new double[ncom]();
        double *comspc_it = new double[njet*nzdum*ncom]();
        double *phoint	= new double[nsyn]();
        double *synabs	= new double[njet*nzdum*nsyn]();
        double *shelen	= new double[nelec]();
        double *shden	= new double[nelec]();
        double *thmshock = new double[nelec]();
        double *thmcomp = new double[nelec]();
        double *thmbase = new double[nelec]();
        double *plshock = new double[nelec]();
        double *plcomp = new double[nelec]();
        double *plbase = new double[nelec]();
        double *tstrm	= new double[nelec]();
        double *drtrm	= new double[nelec]();
        double *phodis	= new double[nsyn]();
        double *nphot	= new double[nsyn]();
        double *jetu	= new double[nsyn]();
        double *disku	= new double[nsyn]();
        double *snu	= new double[nsyn]();
        double *sdump	= new double[nsyn]();
        double *cnu	= new double[ncom]();
        double *cdump	= new double[ncom]();
        double *complot	= new double[NEBIN]();
        double *presyn	= new double[NEBIN]();
        double *postsyn	= new double[NEBIN]();
        double *fplot	= new double[NEBIN]();       
        double *bbplot	= new double[NEBIN]();
        double *total	= new double[NEBIN]();

        /**
         * Synchrotron tables for interpolations
         * **************************************
         * arg[47]
         * var[47]
         *
         */
        static double arg[47]	= {0.0001,0.0002, 0.0005, 0.001, 0.002, 0.005,
				   0.01,  0.03,  0.05,   0.07,  0.1,   0.2,
 		 		   0.3,   0.4,   0.5,    0.6,   0.7,   0.8,
				   0.9,   1.,    1.5,    2.,    2.5,   3.,
				   3.5,   4.,    4.5,    5.,    5.5,   6.,
				   6.5,   7.,    7.5,    8.,    8.5,   9.,
				   9.5,   10.,   12.,    14.,   16.,   18.,
				   20.,   25.,   30.,    40.,   50.};
        static double var[47]	= {-1.002e+0, -9.031e-1, -7.696e-1, -6.716e-1,
				   -5.686e-1, -4.461e-1, -3.516e-1, -2.125e-1,
				   -1.537e-1, -1.192e-1, -8.725e-2, -4.383e-2,
				   -3.716e-2, -4.528e-2, -5.948e-2, -7.988e-2,
				   -1.035e-1, -1.296e-1, -1.586e-1, -1.838e-1,
				   -3.507e-1, -5.214e-1, -6.990e-1, -8.861e-1,
				   -1.073e+0, -1.267e+0, -1.470e+0, -1.670e+0,
				   -1.870e+0, -2.073e+0, -2.279e+0, -2.483e+0,
				   -2.686e+0, -2.893e+0, -3.097e+0, -3.303e+0,
				   -3.510e+0, -3.717e+0, -4.550e+0, -5.388e+0,
				   -6.230e+0, -7.075e+0, -7.921e+0, -1.005e+1,
				   -1.218e+1, -1.646e+1, -2.076e+1};

       /**
         * Jet's parameters
         * ********************************************************
         * mbh		mass of the black hole
         * r_g          gravitational radius
         * eddlum       Eddington luminosity of the black hole
         * inclin	viewing angle wrt to the jet axis
         * dist		distance from the source in kpc
         * jetrat       total power injected in jets
         * r0           radius of the nozzle
         * hratio       aspect ratio of nozzle
         * zsh          z at the shock
	 * zacc		distance travelled by the jet before it
	 * 		stops accelerating with magnetically
	 * 		accelerated profile
	 * zmax         end-point of the jets
	 * eltemp       electron's temperature in the nozzle
         * pspec        non-thermal particle index
         * heat		shock heating; at shock eltem = heat*eltemp
         * brk		break energy of non-thermal particles
         * fsc          acceleration timescale
         * gamfac       maximum Lorentz factor if injection is non-thermal; 
         *		gamax = gamfac*eltemp
         * beta_pl      equipartition factor at base Ue/Ub
         * sigsh	final magnetization at the shock, if velsw > 1
         * tin		temperature of inner disk if > 1, accretion rate 
         *		at inner radius if < 1
         * rin		inner accretion disk radius
         * rout		outer accretion disk radius
         * plfrac       fraction of particles being
         *              accelerated into the power-law
         * mxsw         if 1 injected leptons are thermal, if 0 non-thermal, 
         * 		if between 0 and 1 hybrid
         * velsw	Velocity profile used; change between 
	 *              magnetic and pressure accelerated jets
	 * 		velsw = 0: adiabatic agnjet (never used)	
         *		velsw = 1: isothermal agnjet        		
         *              velsw > 1: magnetically accelerated
         *                         jet Lorentz factor
         * plotsw       switch to choose plotting or not
         * infosw	switch to print useful info on screen
         */
        double mbh, r_g, eddlum, jetrat, pspec, zsh, r0, hratio, h0, eltemp, plfrac, dist, fsc, zmax, equip, beta_pl, gamfac, k_equip, zacc, velsw, sigsh, heat, brk;
        double rin, rout, tin;
        double tbb2,bbf2, bbf1;
        double mxsw;
        int plotsw;
	int infosw;

        mbh       	= param[0];
        r_g             = gconst*mbh*msun/(cee*cee);
        eddlum          = 1.25e38 * mbh;
	inclin          = param[1]*pi/180;
	dist            = param[2]*kpc;        
        jetrat          = eddlum*param[3];
	r0              = param[4]*r_g;
	hratio		= param[5];
	zsh             = param[6]*r_g;
	zacc 		= param[7]*r_g;
	zmax            = param[8];
	eltemp          = param[9] * emerg/kboltz; //Peak electron gamma, here converted into Te (in K)
	pspec           = param[10];
	heat 		= param[11];
        brk		= param[12];
        fsc             = param[13]; 
      	gamfac          = param[14];  
        beta_pl         = param[15];	
	sigsh		= param[16];
	tin 		= param[17]  / kboltz_kev; //Inner temperature of the disk is in keV	
	rin             = param[18] * r_g;
	rout            = param[19]*r_g;
	tbb2		= param[20];
	bbf2		= param[21];
        plfrac          = param[22];
        mxsw            = param[23];
	velsw		= param[24];
	plotsw          = param[25];	
	infosw		= param[26];
	
	h0 = hratio*r0;
	equip = 1./beta_pl; //input paramter is inverse plasma beta; here converted to old equip for calculations; this is more convenient than rewriting all the code from equip to beta_pl
	
	if (infosw == 1) {
		cout << "Nozzle height: " << h0/r_g << " r_g" << endl;
	}

	/**
          * Setting the speed of sound in nozzle <!-- FIXME -->
          */
        gad4_3  = 4./3.;
        betas0  = sqrt(gad4_3*(gad4_3 - 1.)/(gad4_3 + 1.));        
        gam0    = 1./sqrt(1.-(betas0*betas0));
        
        rvel0   = gam0*betas0;
        vel0    = cee*rvel0;

	/**
	  * Larger Comptonization region and no multiple Compton if treating a blazar
	  */

	if (velsw >= gam0 || param[1] <= 12){
		zfrac = 100*zfrac;
		njet = 1;
		Niter = 1;
	}

	/**
	  * Older radiation grids if not treating an aligned source
	  */
	
	if (velsw <= 1 || param[1] > 15) {		
		nsyn    	= 100;
        	ncom    	= 60;
        	snumin  	= log10(0.5*ear[0]/hkev);
        	snumax  	= log10(ear[ne-1]/hkev);
        	cnumin  	= 13;
        	cnumax  	= 23;

        	snuinc  	= (snumax-snumin)/nsyn;
        	cnuinc  	= (cnumax-cnumin)/ncom;
		//zfrac = zfrac/1000;
		//njet = 2;
		//Niter = 500;
		
	}
	
	//This is Riley's grid
	zmin = int(log10(r0))*1.0; // CHANGED TO NO LONGER DEPEND UPON ZINC
    	nz = 70; // 19/07/17, new grid with fixed number of zones
    	double *zed = new double[nz](); // z is now an array that can be called again!
    	int zedcnt; // A counter such that k is not used to determine zed beyond zcut
    	double zinc; // increment is an even log step between zcut and zmax, with nz zones (grid only becomes logarithmic beyond zcut!)
    	int zone_zcut = 0; // Index of zone given by z = zcut, where comptonisation stops


        if (nz > nzdum){
        	cerr << "nz is greater than array size! Please change zmax." << endl;
             	exit(1);
        }

 /**
         * Table of velocity profile solution
         * from Euler equation based on
         * Falcke et al. 1996
         * ***********************************
         * gbx[120]
         * gby[120]
         *
         * Used to be 120 array, now is 54 because
         * it crashes splines routines with the 0s in gbx
         *
         */

    
	static double gbx_vel2[54] = {1.00000, 1.00001, 1.00005, 1.00023, 1.00101, 1.00456, 1.02053, 1.09237,
            1.41567, 2.87053, 6.26251, 14.3691, 34.0825, 82.8831, 205.572, 518.283, 1325.11, 3429.81, 8975.1,
            23719.5, 63255, 170099, 460962, 1.25824e+06, 3.4578e+06, 9.5633e+06, 2.66095e+07, 7.44655e+07,
            2.09528e+08, 5.92636e+08, 1.6846e+09, 4.81148e+09, 1.38055e+10, 3.9787e+10, 1.15153e+11, 3.34651e+11, 9.76408e+11, 2.85981e+12, 8.40728e+12, 2.4805e+13, 7.34419e+13, 2.18185e+14, 6.50341e+14, 1.94471e+15, 5.83352e+15, 1.75523e+16, 5.29701e+16, 1.60321e+17, 4.86616e+17, 1.48111e+18, 4.52032e+18, 1.38326e+19, 4.24394e+19, 1.0e+20};
        
        
        
        static double gby_vel2[54] = {0.485071, 1.05031, 1.05032, 1.05039, 1.05067, 1.05193, 1.05751, 1.08105,
            1.16389, 1.35278, 1.52406, 1.68077, 1.82495, 1.95888, 2.08429, 2.20255, 2.3147, 2.42158, 2.52386,
            2.62205, 2.71662, 2.80793, 2.89628, 2.98195, 3.06516, 3.14611, 3.22497, 3.30189, 3.37702, 3.45046,
            3.52232, 3.59270, 3.66169, 3.72937, 3.79580, 3.86106, 3.92520, 3.98827, 4.05033, 4.11143, 4.17160,
            4.23089, 4.28933, 4.34696, 4.40381, 4.45992, 4.51530, 4.56999, 4.62402, 4.67740, 4.73015, 4.78230,
            4.83388, 4.87281};
        
        
        /* Tabulated corrected velocity profile */
        bool isIsoTh = 0;
       // int velsw = 1; // velsw = 1 full 1D quasi-isothermal Bernoulli Eq.;
	if ((velsw > 1) & (velsw <= gam0)){
		cout << "The model is running with BLJet profile" << endl;
		cout << "Velsw now represents the final Lorentz factor of the jet after acceleration" << endl;
		cout << "Currently velsw < gam0! Fix your input parameters!" << endl;
		exit(1);
	// velsw > 1 magnetically accelerated Bernoulli Eq.; 
	}
	if (velsw == 0) {
		cout << "The model is running with adiabatic agnjet profile, careful!" << endl;
	}
	
       	int sizegb = 55;
        double zeta = 1; // Uj = n mc2
        double gbx[sizegb-1],gby[sizegb-1];
	if (velsw > 1){
		if (mxsw == 0){ //with BLJet profile, if purely non thermal particles then never re-shock
			zsh = pow(10.,zmax) + (r0+h0);
		}
		double step = (log10(pow(10.,zmax)/r_g + 1. -log10(pow(10.,zmin)/r_g)))/(sizegb-3.);
		
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
        /*
         * velsw = 0 ADIABATIC AGNJET;
         * velsw = 1 AGNJET;
         * velsw > 1 BLJET (this needs to be greater than g0 = 1.1);
         */
    
        /**
         * Files for plotting
         *
         */
        ofstream fluxplotFile, complotFile, presynFile, postsynFile, bbplotFile, zonesFile;

        if(plotsw == 1){
		fluxplotFile.open("outputs/total.dat",ios::trunc);
		if(!fluxplotFile.is_open()){
			cerr << "*** Error: can't open outputs/total.dat" << endl;
			exit(1);
		}
		else{
			fluxplotFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Total Emission" << endl;
		}
		fluxplotFile.close();

		complotFile.open("outputs/com.dat",ios::trunc);
		if(!complotFile.is_open()){
			cerr << "*** Error: can't open outputs/com.dat" << endl;
			exit(1);
		}
		else{
			complotFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Compton Emission" << endl;
		}
		complotFile.close();

		presynFile.open("outputs/presyn.dat",ios::trunc);
		if(!presynFile.is_open()){
			cerr << "*** Error: can't open outputs/presyn.dat" << endl;
			exit(1);
		}
		else{
			presynFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Presyn Emission" << endl;
		}
		presynFile.close();

		postsynFile.open("outputs/postsyn.dat",ios::trunc);
		if(!postsynFile.is_open()){
			cerr << "*** Error: can't open outputs/postsyn.dat" << endl;
			exit(1);
		}
		else{
			postsynFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Postsyn Emission" << endl;
		}
		postsynFile.close();

		bbplotFile.open("outputs/bb.dat",ios::trunc);
		if(!bbplotFile.is_open()){
			cerr << "*** Error: can't open outputs/bb.dat" << endl;
			exit(1);
		}
		else{
			bbplotFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Black-body Emission" << endl;
		}
		bbplotFile.close();

		if(infosw == 1){

			zonesFile.open("outputs/zones.dat",ios::trunc);
			if(!zonesFile.is_open()){
				cerr << "*** Error: can't open zones.dat" << endl;
				exit(1);
			}
			else{
				zonesFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Total Emission" << endl;
			}
			zonesFile.close();
		}
	}

        /*********************************
         * STARTING RUNNING xrbjet model *
         *********************************/
        //cout << "Running BH jet model..." << endl;

        /**
         * Define energy array from ear array, and then convert into frequencies
         * to be compatible with old routines.
         *
         * ear, ebin in [KeV]
         *
         */
        for(i=0; i<(ne-1); i++){
        	ebin[i] = ear[i] + (ear[i+1]-ear[i])/2;
                nutot[i]= log10(ebin[i]/hkev);
	}

        /**
         * Define radiation grids
         *
         */
        //Synchrotron
        for(i=0; i<nsyn; i++){
                nubb[i] = snumin + (i+0.5)*snuinc;
                nurad[i]= pow(10., nubb[i]);
                energ[i]= herg*nurad[i];
                ephot[i]   = log10(energ[i]);
//                cout << energ[i] << endl;
        }
        //Compton
        for(i=0; i<ncom; i++){
                ephxr[i]= pow(10, cnumin+(i+0.5)*cnuinc);
	}

        /* New electron energy/dist arrays for printing to output file */
        double *tdyntot = new double[nz]();
        double *tsyntot = new double[nz]();
        double *ebrarr = new double[nz]();
    
        /**
         * Initialization of the synchrotron interpolation
         *
         */
        gsl_interp_accel *acc_syn       = gsl_interp_accel_alloc();
        gsl_spline *spline_syn          = gsl_spline_alloc(gsl_interp_cspline, 47);     //XXX HERE change interpolation function -- 47 == size of arg[] and var[]
        gsl_spline_init(spline_syn, arg, var, 47);

	/*
         * Checking if disk exists
         *
         */
        outfac  = rout/rin;
        if(outfac < 1){
                disksw  = 0;
                if (infosw ==1) {
                	cout << "rin > rout, no disk component!" << endl;
                }
        }

        /**
         * Setting disk temperature/checking if it's super Eddington
         * Setting normalizations of disk and second blackbody
         */
         
        if (tin > 0) { 
        	thrlum  = 4*pi*sbconst*pow(tin,4)*pow(rin,2)*(1-1/outfac);
        	eddrat  = thrlum/eddlum;
        	if(eddrat > 1){
                	cout << "Disk flux super-Eddington, please readjust rin/tin. eddrat larger than 1: " << eddrat << endl;                
        	}       
        }
        
        else {
        	eddrat = - tin * kboltz_kev;
        	thrlum  = eddrat*eddlum;
		tin = pow(thrlum/(4*pi*sbconst*pow(rin,2)*(1-1/outfac)),0.25);	
        }	
        
        bbf1 = sqrt(bbf2*jetrat/(sbconst*pi*pow(tbb2,4)*pow(rout,2)));
        reff = bbf1*rout;
        reff2 = reff*reff;

        /* Setting thmfrac as 1-plfrac */
        thmfrac = 1.-plfrac;

        //[M.Nowak] H/R~L/L_edd
        hbb      = rin*eddrat;

        /* distance in the jet above which Comptonization is neglected*/
        zcut    = zfrac*r0;	

        /**
         * Calculating internal energy in electron's distribution (Maxwellian-J\¨uttner)
         *
         * need mj_aven.cc <!-- FIXME -->
         */
        emin    = 2.23*kboltz*eltemp;
        gamin   = emin/emerg;
        emax    = gamfac*gamin*emerg;
	ebreak = emax + 1.;
	
        uplim   = pow(10, -54.32 + 9.416*log10(eltemp) - 0.376*log10(eltemp)*log10(eltemp));
        endncom = enden(eltemp, uplim);
        endnsmj = endncom*emerg;
	
	/*
	* For the magnetically accelerated jet we need to compute the inital particle number differently, because U_e + U_b = U_p doesn't hold
	* In this case, for pspec close to 2 and it is possible for numerical issues to arise, so in this case we compute the average Lorentz 
	* factor as if pspec = 2. The difference is negligible
	*/

	if (velsw <= 1){
        	nprot0  = (jetrat/4.)/(vel0*pmgm*cee*cee*pi*r0*r0);
	} else {
		if (mxsw == 1){
			avgen = endncom;
		} else if (mxsw == 0){
			if (pspec >= 1.99 && pspec <= 2.01){
				avgen = log(emax/emin)*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,1.-pspec))*emerg);
			} else {
				avgen = (pow(emax,2.-pspec)-pow(emin,2.-pspec))*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,1.-pspec)*(2.-pspec))*emerg);
			}
		} else {
			if (pspec >= 1.99 && pspec <= 2.01){
				avgen = mxsw*endncom+(1.-mxsw)*log(emax/emin)*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,1.-pspec))*emerg);
			} else {
				avgen = mxsw*endncom+(1.-mxsw)*(pow(emax,2.-pspec)-pow(emin,2.-pspec))*(1.-pspec)/((pow(emax,1.-pspec)-pow(emin,1.-pspec)*(2.-pspec))*emerg);
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
		eta = (sig0*pmgm)/(avgen*emgm*(2.*equip-sig0*gad4_3));
		numcorr = 1.+eta*(1.+equip)*(avgen*emgm)/(pmgm);
		nprot0  = jetrat/(2.*numcorr*vel0*pi*r0*r0*pmgm*cee*cee); 
		cout << "Changed equip to: " << equip << endl;
		cout << "Check equip value! Pair content now is: " << eta << endl;
	}
	if (nprot0 < 0){
		cout << "Proton number density: " << nprot0 << endl;
		cout << "Too many pairs, the jet has a negative proton number density!" << endl;		
	}
	
        /**
         * Calculating magnetic energy from equipartition assumptions
         * ***********************************************************
         * @param mxsw          maxwellian switch [jet parameter]
         * @param nenpsw        new parameter to allow different partitions
         * @param equip         equipartition factor [jet parameter]
         * @param pspec         power-law index [jet parameter]
         * 
         * @return cnorm        normalization factor of the power-law distribution <!-- FIXME -->
         * @return ntot0        initial number of electron
         * @return b_en         magnetic energy in nozzle
         *
         */
        equipartition(velsw, mxsw, nenpsw, eta, equip, pspec, nprot0, emin, emax, endnsmj, cnorm, ntot0, b_en);
        b0 = sqrt(8.*pi*b_en);

        /**
         * Frequencies in Hz, fluxes in mJy
         * Setting gamax0 to unity so that gamax/gamax0 gives
         * the correct Mach number factor for energy shift
         *
         */
        gamax0  = 1.;

    
        /* generate maxwellian/power law distirbutions here based on mxsw = 1 or mxsw = 0 */
        if(mxsw == 1) {
        	endens  = endnsmj*ntot0;
            	betat   = emerg/(kboltz*eltemp);
           	enorm   = ntot0*betat/(emerg*k2_fnc(betat));
           	elmin	= emerg;
          	elmax	= 50.*kboltz*eltemp;
          	einc	= log10(elmax/elmin)/nelec;
            
            for(i=0; i<nelec; i++){
                rdlgen[i] = log10(elmin)+(i+0.5)*einc;
                rden[i]	= pow(10.,rdlgen[i]);
                game	= rden[i]/emerg;
                if(game < 1.e-3){
                	bete	= sqrt(game*game - 1.)/game;
                }
                else{
                	bete	= 1.-1./(2.*game*game)-1./(8.*game*game*game*game);
                }
                reled	= enorm*game*game*bete*exp(-game*betat);
                rdend[i]= reled*rden[i];
                rdedn[i]= log10(reled);
            }
        }
	else if(mxsw == 0) {
       		if(pspec >= 1.99 && pspec <= 2.01){ //values too close to 2 for the slope give numerical errors, so we treat them as if they were 2
                	endens	= cnorm*log(emax/emin);
            	}
            	else{
            		endens	= cnorm*(pow(emax,2.-pspec) - pow(emin,2.-pspec))/(2.-pspec);
            	}
            	einc	= log10(emax/emin)/nelec;
            	for(i=0; i<nelec; i++){
                	rdlgen[i] = log10(emin) + (i+0.5)*einc;
                	rden[i]	= pow(10,rdlgen[i]);
                	reled	= cnorm*pow(rden[i],-pspec);
                	rdend[i]= reled*rden[i];
                	rdedn[i]= log10(reled);
            	}
        }
        
        /* This case represents a mixed thermal/non-thermal injected particle distribution (i.e. mxsw != 0 && mxsw != 1). We start by setting up a thermal distribution as we would when
         mxsw = 1, with the energy array extending from elmin to elmax, where elmin is game = 1, elmax is 50kT. We assign these electrons to variable edtrm. The energy array is then
         interpolated to fit the new grid, elmin - emax (emax is set by gamfac, it is the maximum energy of the power-law distirbution), using the gsl spline. Then we calculate the
         power-law distribution, and append both edtrm and pltrm to the full distribution, as well as setting up the arrays thmcomp and plcomp for future use. */
        
        else {
            
       		endens  = endnsmj*ntot0;
       		betat   = emerg/(kboltz*eltemp);
        	enorm   = mxsw * ntot0 * betat/(emerg*k2_fnc(betat));
            	elmin	= emerg;
            	elmax	= 50.*kboltz*eltemp;
            	einc	= log10(emax/elmin)/nelec;
            
            	for(i=0; i<nelec; i++){
                	rdlgen[i] = log10(elmin)+(i+0.5)*einc;
                	rden[i]	= pow(10.,rdlgen[i]);
                	if(rden[i] <= elmax) {
                    		game	= rden[i]/emerg;
                    	if(game < 1.e-3){
                        	bete	= sqrt(game*game - 1.)/game;
                    	}
                    	else{
                    		bete	= 1.-1./(2.*game*game)-1./(8.*game*game*game*game);
                    	}
                    	edtrm	= enorm*game*game*bete*exp(-game*betat);
                }
                else {
                	edtrm = 0.;
		}
                if(rden[i] >= emin) {
               		pltrm = cnorm * pow(rden[i], -pspec);
                }
                else{
			pltrm = 0.;
                }
                
                thmbase[i] = edtrm;
                plbase[i] = pltrm;
                
                reled = edtrm + pltrm;
                rdend[i] = reled * rden[i];
		rdedn[i] = log10(reled);
                }
	}
	
        /**
         * if running with nenpsw=0, print out the jetrat (jet power) you would have if running with nenpsw=1. This can be used to calculate what new jetrat value you need when using an old parameter file jetwrap
         */
        
	if(nenpsw == 0) {
		printf("old jetrat=%e\n",jetrat/eddlum);
		jetrat = (endens + b_en) * (4.*vel0*pi*r0*r0);
		printf("Renormalized jetrat using old values=%e\n", jetrat/eddlum);
            
            /**
             *Check new normalizations
             **/
          	printf("new ntot0, nprot0 =%e,%e\n",ntot0,nprot0);
		printf("Ue,Ub=%e,%e\n", endens,b_en);
		printf("Check 2: Up,Ue+Ub,Ub/Ue=%e,%e,%e\n",nprot0*pmgm*cee*cee,endens+b_en,b_en/endens);
            	printf("new bfield=%e\n",sqrt(b_en)*8.*pi);
            	printf("endens,new U_b=%e,%e\n",endens,equip*endens);
	}
    
        /**
         * Initializing jet
         *
         */
        gsl_interp_accel *acc_jet       = gsl_interp_accel_alloc();
        gsl_spline *spline_jet          = gsl_spline_alloc(gsl_interp_akima, sizegb-1);
    
        if(!isInitialized){
		gsl_spline_init(spline_jet, gbx, gby, sizegb-1);
                if (velsw > 1){
			gbs0 = gsl_spline_eval(spline_jet, gbx[0]+1., acc_jet);
		} else {
			gbs0 = gsl_spline_eval(spline_jet, pow(10,0.0), acc_jet);
		}
        }
	
        for(k=0; k<nz; k++){
		for (i=0;i<ncom;i++){
                	for (m=0;m<njet;m++){
				comspc_it[(m*nz+k)*ncom+i] = 0.0;
                	}
		}
        }
        
	   
	if (infosw ==1) {
		cout << "Physical quantities: " << endl;
		cout << "Maxwell-Juttner peak in the nozzle: " << 2.23*eltemp << endl;
		cout << "Inner disk temperature in K: " << tin << "; accretion rate:" << eddrat << endl;
		cout << "Pair content: " << eta << endl;
	}
    
    
        //--------------------------------------------------------------------------------------------
    
        /**
         * START THE BIG LOOP OVER Z
         *
         */

	for(k=0; k<nz; k++){
		
		/*Adjusted grid: delz = r in inner regions to avoid resolution issues, logarithim grid in outer regions*/
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
                		zinc = (zmax - log10(zcut))/(nz - zone_zcut);
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
            		}
            		else{
                		zinc = (zmax - log10(zcut/80.))/(nz - zone_zcut);
                		z = pow(10.,log10(zcut/80.) + zinc*(k-zone_zcut)) * (1. + pow(10.,-6));
                		delz = z - z*pow(10.,-zinc);
				if (infosw == 1) {
					cout << "Grid shifted; distance: " << z/r_g << " r_g; zone number " << k  << endl;
					if (z >= zcut) {
						cout << "Out of the Comptonization region; distance: " << z/r_g << endl;
					}	
				}
				
				
                	}
		}
		
           	zed[k] = z;


                /**
                 * Define properties of jet's segment at z, given initial conditions
                         * ******************************************************************
                 *
                 * @param z		Distance along the jet from black hole
                 * @param b0		Magnetic field in Nozzle
                 * @param ntot0		Electron number density in Nozzle
                 * @param gamax0	Electron Lorentz factor in Nozzle
                 * @param r0		Width of Nozzle
                 * @param h0		Height of Nozzle
                 * @param spline_jet    spline interpolation objet
                 * @param acc_jet       accelerator interpolation object
                 *
                 * @return rvel		Bulk velocity of jet (gamma*beta)
                 * @return bfield	Magnetic filed at z
                 * @return ntot		Density at z
                 * @return gamax	Electron Lorentz factor at z
                 * @return r		Width at z
                 *
                 */
		
		if (velsw <= 1){
			jetpars(velsw, zeta, zmin, z, zsh, r_g, b0, ntot0, gamax0, r0, h0, spline_jet, acc_jet, rvel, drvel, bfield, ntot, gamax, r);
		} else {
			bljetpars(mxsw, velsw, z, zacc, r_g, eta, nprot0, gamax0, r0, h0, endnsmj, pspec, cnorm, emin, emax, ebreak, spline_jet, acc_jet, rvel, bfield, ntot, gamax, r, sigsh); 
		}	
	

		if (infosw == 1) {
			cout << "Z: " << z/r_g << " B: " << bfield << " n: " << ntot << " gmin: " << emin/emerg << " gbreak: " << gebreak << " gmax: " << gemax << endl;
			cout << "velocity: " << rvel << " Doppler factor: " << 1./(sqrt(1.+rvel*rvel)*(1. - beta*cos(inclin))) << " radius: " << r/r_g << " opening angle: " << 57.3*2.*atan(r/z) << endl << endl;	

		}

                if(isVerbose){
                	cout << left << "Jetpars done... pars are: " << setw(15) << z/r_g << setw(15) << rvel << setw(15) << bfield << setw(15) << ntot << setw(15) << gamax << setw(15) << r/r_g << endl;
                }

                /* area and volume of jet's segment */
                area    = pi*r*r;
                vol     = delz*area;

                /* Express gamma and beta in jet's segment <!-- FIXME better definition --> */
                gamv2   = 1.+rvel*rvel;
                gamv      = sqrt(gamv2);
                if (gamv2 > 1.e5){
                        beta	= 1. - 1./(2.*gamv2) - 1./(8.*gamv2*gamv2);
                }
                else{
                        beta	= sqrt(gamv2-1.)/gamv;
                }
                
              /**
                * theff is the angle to the irradiation blackbody. Then calculate relevant
                * effective temperature (see R&L 153 or 156, iv)
                *
                */
               
                theff	= atan(reff/(z-hbb/2.));
                tbbeff	= tbb2*gamv*(1.-beta*cos(theff));

                /* Setting Doppler factor */

		for(l=0; l< njet; l++){
			dopfac[l*nz+k]	= 1./(gamv*(1. - beta*cos(inclin)*pow(-1.,l)));
		}

                /* Express gshift and nw, used in shifting down energy distribution */
                gshift  = gamax/gamax0;
                nw      = 0;

                /**
                 * Depending where you are in the jet, different actions need to be done
                 *
                 */
                /* Here we define the process of shifting down the distributions when we have a mixed particle distribution injected into the base of the jet.
                In this case we must shift down the maxwellian component only, tag on the power law as before (with an updated break energy), and renormalise everything */
            
                if(mxsw != 1 && z < zsh) {
                    
			/*if (log10(zsh) < zmax){
                        	cout << "zsh=" << log10(zsh) << ", zmax= " << zmax << endl;
                        	cerr << "WARNING! Set zsh greater than zmax!" << endl;
                        	exit(1);
                    	} else {*/
//                        cout << "here there is the mxsw: is this working?!" << endl;
                        	ub = bfield * bfield/(8. * pi);
                        	syncon	= 4./3. * sigtom * ub/(emgm*emgm*cee*cee*cee);
                        if(z < h0) {
                            	tdyn = h0/(beta*cee);
                            	tdyntot[k] = tdyn;
                            	tsyntot[k] = (1./syncon);
                            	ebrarr[k] = tsyntot[k]/tdyntot[k];
                            	ebreak = ebrarr[k];

                        }
                        else {
                            	tdyn    = r / (beta * cee);
                            	tdyntot[k] = tdyn + tdyntot[k-1];
                            	tsyntot[k] = (1./syncon) + tsyntot[k-1] ;
                            	delebreak = tsyntot[k]/tdyntot[k];// - tdyntot[k-1]/tsyntot[k-1];
                           
                            	ebrarr[k] += delebreak;
                            	ebreak = ebrarr[k];
			}
                        gebreak = ebreak/emerg;

			//cout << "tdyn: " << tdyn << " tsyntot: " << tsyntot[k] <<  " gebreak: " << gebreak << endl;
                     	
                        if (!isIsoTh) {
                            	mjteff	= eltemp*gshift;
                        } else {
                            	mjteff = eltemp;
                        }
                        
                        emin	= 2.23*kboltz*mjteff; // changed from 2. to 2.23 ....this should be correct, right?
                        mjpeak	= pow(10.,9.514 - 1.028*log10(mjteff));
                        
                        if(ebreak <= emin) {				
			    	cnorm = (1.-mxsw)*ntot*(pspec)/(pow(emin, -pspec) - pow(emax, -pspec));  // if break energy less than minimum energy, set normalisation for cooled distribution
                        }
                        else if(ebreak >= emax){
			   	cnorm = (1.-mxsw)*ntot*(1.-pspec)/(pow(emax, 1.-pspec) - pow(emin, 1.-pspec)); // if break greater than max energy, set normalisation for uncooled distribution
                        }
                        else{
                            	cnorm = (1.-mxsw)*ntot/((pow(ebreak, 1.-pspec)/(1.-pspec)) - (pow(emin, 1.-pspec)/(1.-pspec)) + (pow(ebreak, 1.-pspec)/pspec) - (ebreak*pow(emax, -pspec)/pspec));
                            // If break lies within energy limits, set normalisation for broken lower law
                        }
                        
                        if(ebreak <= emin) {
                            	enorm = cnorm*pow(emin,-(pspec+1.))*betat/(k2_fnc(betat)*mjpeak); // if break energy < min energy, normalisation at boundary between thermal/power law distributions is product of thermal and cooled power law
                        }
                        else{
                            	enorm = cnorm*pow(emin, -pspec)*betat/(k2_fnc(betat)*mjpeak); // if break energy > min energy, normalisation at boundary between thermal/power law distributions is product of thermal and uncooled power law
                        }
          
                        for(i=0; i<nelec; i++){
                            
				etemp[i] = rden[i];
                            	if(etemp[i] >= emin){
                                	if(ebreak <= emin){
                                    		pltrm = cnorm*pow(etemp[i],-(pspec+1)); // if break energy < min energy, power law distribution completely cooled
                                    
                                }
                                else if(ebreak >= emax){
                                    	pltrm = cnorm*pow(etemp[i],-pspec); // if break energy > emax, power law distribution completely uncooled
                                }
                                else{
                                    	if(etemp[i] < ebreak){
                                        	pltrm = cnorm*pow(etemp[i],-pspec); // if break energy inside energy limits, broken power law distribution with break at ebreak
                                    	}
                                    	else{
                                        	pltrm = cnorm*ebreak*pow(etemp[i],-(pspec+1));
                                    		}
                                	}
                      		}
                            	else{
                                	game = etemp[i]/emerg;
                                	if(game <= 1e3){
                                    		bete = sqrt(game*game-1.)/game;
                                	}
                                	else{
                                    		bete = 1.-1./(2.*game*game)-1.*(8*game*game*game*game);
                                	}
                                	pltrm = enorm * game*game*bete * exp(-game*betat); // was game/betat, which I think was wrong, changed 6 may 2015
                            	}
                            	plcomp[i] = log10(pltrm);
                            
                        }
                        
                        nw = 0;
                        
                        
                        for(i=0; i<nelec; i++) {
                            	lelec[i] = rdlgen[i];
                            	elen[i] = pow(10.,lelec[i]);
                            	etemp[i] = log10(gshift * elen[i]);
                            	if(gshift * elen[i]/emerg <= 1){
                                	nw = nw + 1;
                            	}
                            	if(thmbase[i] > 0){
                                	thmcomp[i] = log10((ntot/ntot0)*thmbase[i]);
                            	}
                            	else{
                                	thmcomp[i] = 0.;
                            	}
                            
                        }
                        
                        sum1 = 0.;
                        for(i=0; i<nelec-1; i++) {
                        	if(thmcomp[i] == 0. || thmcomp[i+1] == 0.){
                                	tmp1 = 0;
                            	}
                            	else{
                                	tmp1 = (pow(10.,etemp[i+1]) - pow(10.,etemp[i])) * 0.5 * (pow(10.,thmcomp[i]) + pow(10.,thmcomp[i+1]));
                            	}
                            	sum1 += tmp1;
                        }
                        
                        renorm = log10((ntot*mxsw)/sum1);
                        
                        
                        for(i=0; i<nelec; i++) {
                            	if(thmcomp[i] == 0.){
                                	thmcomp[i] = 0.;
                            	}
                            	else{
                                	thmcomp[i] += renorm;
                            	}
                            
                            	if(etemp[i] > emin && plcomp[i] > (enorm * game*game*bete * exp(-game*betat))){
                                	plcomp[i] = thmcomp[i];
                            	}
                            	if(thmcomp[i] == 0 || plcomp[i] == 0){
                                	if(thmcomp[i] == 0.){
                                    		dtemp[i] = log10(pow(10.,plcomp[i]));
                                	}
                                	else{
                                    		dtemp[i] = log10(pow(10.,thmcomp[i]));
                                	}
                            	}
                            
                            	else{
                                	dtemp[i] = log10(pow(10.,thmcomp[i]) + pow(10.,plcomp[i]));
                            	}
                        }
                        
                        
                        /* Count the lost electrons and renormalise accordingly */
                        
                        if(nw>0){
                            	if(nw==1){
                                	totlos  = (pow(10.,etemp[1])-pow(10.,etemp[0]))*pow(10.,dtemp[0]);
                            	}
                            	else if(nw > 1){        //<!-- FIXME why not simply else? -->
                                	totlos  = 0;
                               		totlos  = 0;
                                	for(i=0; i < nw-1; i++){
                                    		tmlos = (pow(10.,etemp[i+1])-pow(10.,etemp[i]))*0.5*(pow(10.,dtemp[i])+pow(10.,dtemp[i+1]));
                                    		totlos += tmlos;
                                	}
                            	}
                            	oldnum = ntot;
                            	ntot -= totlos;
                            	renorm = log10(ntot/oldnum);
                            	for(i=0; i<nelec; i++){
                                	dtemp[i] += renorm;
                            	}
                            
                            	gsl_interp_accel *acc_shift  = gsl_interp_accel_alloc();
                            	gsl_spline *spline_shift     = gsl_spline_alloc(gsl_interp_cspline, nelec);
                            	gsl_spline_init(spline_shift, etemp, dtemp, nelec);
                          
                            	einc    = (etemp[nelec-1]-etemp[nw])/nelec;
                            	for(i=0; i<nelec; i++){
                                	lelec[i]= etemp[nw]+(i+0.5)*einc;
                                	ledens[i] = gsl_spline_eval(spline_shift, lelec[i], acc_shift);
                                	elen[i] = pow(10., lelec[i]);
                                	eled[i] = elen[i]*pow(10., dtemp[i]);
                            	}
				
                            	gsl_spline_free(spline_shift), gsl_interp_accel_free(acc_shift);
                        }//END nw > 0
                        
                        else{
                            	for(i=0; i<nelec; i++){
                                	lelec[i]= etemp[i];
                                	elen[i] = pow(10.,lelec[i]);
                                	eled[i] = elen[i]*pow(10.,dtemp[i]);
                                	ledens[i]= dtemp[i];
                                
                            	}
                        }//END else nw < 0
                        
                        sum1 = 0.;
                        for(i=0; i<nelec-1; i++) {
                            	tmp1 = (pow(10.,lelec[i+1]) - pow(10.,lelec[i])) * 0.5 * (pow(10.,ledens[i]) + pow(10.,ledens[i+1]));
                            	sum1 += tmp1;
                        }
                        renorm	= log10(ntot/sum1);
                        for(i=0; i<nelec; i++){
                            	ledens[i] += renorm;
                        	}
                  	//}
                }
            
                if(!isShock || z<=zsh){
            	
		/******************
                * NOZZLE + SHOCK *
                ******************/
                /**
                * Take read in distribution, Shift the energies down by mach**{-1/3},
                * keep track of lost ones by renormalizing first and the losing
                * the first bins... and renormalizing
                */
                	for(i=0; i<nelec; i++){
                       		lelec[i]= rdlgen[i];
                       		elen[i]	= pow(10.,lelec[i]);
                       		ledens[i] = rdedn[i];
                       		ledens[i] = log10(ntot/ntot0*pow(10.,ledens[i]));
                	}
                	for(i=0; i< nelec; i++){
                		etemp[i]= log10(gshift*elen[i]);
                        	if(gshift*elen[i]/emerg <= 1.){
                        		nw = nw+1;
                        	}
                        	dtemp[i]= ledens[i];
                    		}
			sum1 = 0.;
                	for(i=0; i < nelec-1; i++){
                		tmp1 = (pow(10.,etemp[i+1])-pow(10.,etemp[i])) * 0.5 * (pow(10.,dtemp[i])+pow(10.,dtemp[i+1]));
                        	sum1 += tmp1;
                	}
                	renorm = log10(ntot/sum1);
                	for(i=0; i<nelec; i++){
                		dtemp[i] += renorm;
               		}
                	/**
                	* spline_dbl and resize the array into the old lelec, ledens
                	* if you lose particles at low energy
                	*/
                	if(nw>0){
                		if(nw==1){
                        		totlos = (pow(10.,etemp[1])-pow(10.,etemp[0]))*pow(10.,dtemp[0]);
                            		}
                            		else if(nw > 1){        //<!-- FIXME why not simply else? -->
                                		totlos  = 0;
                                		totlos  = 0;
                                		for(i=0; i < nw-1; i++){
                                        		tmlos = (pow(10.,etemp[i+1])-pow(10.,etemp[i]))*0.5*(pow(10.,dtemp[i])+pow(10.,dtemp[i+1]));
                                                	totlos += tmlos;
                                    		}
                            		}
                            		oldnum = ntot;
                            		ntot -= totlos;
                            		renorm  = log10(ntot/oldnum);
                            		for(i=0; i<nelec; i++){
                                		dtemp[i] += renorm;
                            		}
				
                            		gsl_interp_accel *acc_shift  = gsl_interp_accel_alloc();
                            		gsl_spline *spline_shift     = gsl_spline_alloc(gsl_interp_cspline, nelec);
                            		gsl_spline_init(spline_shift, etemp, dtemp, nelec);
				
                            		einc    = (etemp[nelec-1]-etemp[nw])/nelec;
                            		for(i=0; i<nelec; i++){
                                    		lelec[i]= etemp[nw]+(i+0.5)*einc;
                                    		ledens[i] = gsl_spline_eval(spline_shift, lelec[i], acc_shift);
                                    		elen[i] = pow(10., lelec[i]);
                                    		eled[i] = elen[i]*pow(10., dtemp[i]);
                            		}
				
				gsl_spline_free(spline_shift), gsl_interp_accel_free(acc_shift);
			}//END nw > 0
			else{
                		for(i=0; i<nelec; i++){
                        		lelec[i]= etemp[i];
                                	elen[i] = pow(10.,lelec[i]);
                                	eled[i] = elen[i]*pow(10.,dtemp[i]);
                                	ledens[i]= dtemp[i];
				}
				}//END else nw >0
                    
                		sum1 = 0.; // ADDED 22/03/16, to renormalise the electron distribution after shifting the energy array. This then gives the ratio ntot/sum1=1
               			for(i=0; i<nelec-1; i++) {
                			tmp1 = (pow(10.,lelec[i+1]) - pow(10.,lelec[i])) * 0.5 * (pow(10.,ledens[i]) + pow(10.,ledens[i+1]));
                			sum1 += tmp1;
                		}
                		renorm = log10(ntot/sum1);
               	 		for(i=0; i<nelec; i++){
                			ledens[i] += renorm;
                		}
//                    cout << "Tot los particles =" << totlos << endl;
			
                        // <!-- FIXME insert shift distribution of protons here -->
//                    cout << "Noshock?!" << z << " " << zsh << endl;
		if(z>zsh){

                	/*********
                        * SHOCK *
                        *********/
                        /* Here at the shock, in contrast to the previous version of code, we need to calculate the normalisations of a thermal + broken power-law electron distribution function, subject to whether the break energy lies within the limits of the power law distirbution (so, emin, emax, and ebreak). The thermal distirbution is read in from the electron distribution as it is at the shock. The minimum energy is set to the peak temperature of the MJ distribution (2.23*k*mjteff), and the maximum energy is calculated by equating the acceleration/cooling rates. The break energy is calculated by equating the synchrotron cooling time with the escape time of the electrons in the zone. Then the normalisation conditions for the power are set according to the 3 cases: (i) ebreak < emin, (ii) ebreak > emax, (iii) emin < ebreak < emax. The power law distribution itself is then set for whatever case is satisfied, and this gives both thermal and power law distributions at the shock, given by thmshock and plshock. The total distribution is the log of the sum. The arrays thmshock and plshock can then be read in at every jet segment for increasing z over the loop. */
                        
                        
                        // cout << "I am inside shock region!/n" << endl;
//                            cout << "zsh:" << zsh << ", zmax:" << zmax << endl;
                            //N.B.: to be run ONLY ONCE
			if(isVerbose){
                        	cout << "THIS SHOULD ONLY APPEAR ONCE FOR z =" << z << " (k = " << k << ")" << endl;
                        }
//                            cout << "THIS SHOULD ONLY APPEAR ONCE FOR z =" << z << " (k = " << k << ")" << endl;
                        gshock  = gshift;
                        //cout << "SHOCK!!!!" << endl;  

			if (isHeated == false) {
			//if (velsw != 0 && velsw != 1 && isHeated == false) {
				isHeated = true;
                                eltemp = heat*eltemp; //NEW: this accounts for shock heating, if running with BLJet profile
				if (infosw ==1) {
					cout << "Electron temperature after heating: " << 2.23*eltemp << endl  << endl;
				}
				if (2.23*eltemp > 1.e12 && inclin >= 10) { // 12/09/12: this keeps the code from crashing with high densities and high temperatures, ie for non-blazar sources
					eltemp = 5.e11;
					cout << "Maxed out shock heating at 5x10^11!" << endl; 
				}				
			}		
	  
                        /**
                        * To get idea of SSC rate, use flux from previous component
                        * at distance delz but it will be overestimate.
                        */
                        if(isVerbose){
                        	cout << left << setw(15) << "ACCELPARS:" << setw(15) << r << setw(15) << z << setw(15) << dopfac[0*nz+k] << setw(15) << dopfac[1*nz+k] << setw(15) << bfield << setw(15) << ntot << endl;
                        }

                        ub = bfield*bfield/(8.*pi);

                        sum1    = 0.;
                        for(m=0; m < njet; m++){
                               	for(i=0; i < nsyn; i++){
                                	phofrq[i]= pow(10.,nusyn[(m*nz+(k-1))*nsyn+i])/dopfac[m*nz +(k-1)];
                                    	//3 powers of dopfac because one was inside synint for angle aberration
                                    	phoint[i]= mjy * pow(10., synabs[(m*nz+(k-1))*nsyn+i])/pow(dopfac[m*nz +(k-1)],3);
                                    	//Handle differently: not flux but assuming that the synchrotron
                                    	//was produced in this segment. That's overestimating but hopefully not
                                    	//too much because r, delz also increased.
                                    	phoint[i]*= 4. * dist*dist/(r*r*cee);
			}
                        //Rough integration to get approx energy density
                        for(i=0; i< nsyn-1; i++){
				tmp1 = (phofrq[i+1]-phofrq[i])*0.5*(phoint[i]+phoint[i+1]);
                        	sum1 += tmp1;
                        	}
			}
                        ucom    = sum1;

                        if(isVerbose){
                        	cout << left << setw(15) << "ucom:" << setw(15) << ucom << endl;
                        }

                        /**
                        * Add in multicolor disk contribution
                        * plus component from irradiation, w/ effective area pi*reff**2
                        */
			if(disksw==1){
				bbdil = 2.*reff2/(2.*reff2+pow(z-hbb/2.,2))-reff2/(reff2+pow(z-hbb/2.,2));
				if(bbsw==1){
                                	nubot = 1.5e-5*kboltz*tin/herg;
                                        nutop = 8.*kboltz*tin/herg;
                                        if(bbdisk(log(nutop), gamv, tin, rin, rout, z, hbb)*nubot/(bbdisk(log(nubot), gamv, tin, rin, rout, z, hbb)*nutop) > 1.e30){
                                        	cerr << "nutop not high enough" << endl;
                                                exit(1);
                                     	}
                                        tst	= bbdisk(log(nutop), gamv, tin, rin, rout, z, hbb);
                                        if(tst <= 1.e-20){
                                        	while(tst < 1.e-20){
                                                	nutop *= 0.75;
                                                        tst = bbdisk(log(nutop), gamv, tin, rin, rout, z, hbb);
                                                }
                                    	}
                                        bbintegrals(log(nubot), log(nutop), gamv, tin, rin, rout, z, hbb, bbint);
				    	uphdil = bbint/cee;
				    	uphdil2	= bbdil*aconst*tbbeff*tbbeff*tbbeff*tbbeff/4.;
   					if(z > hbb/2.){
                                               	ucom = ucom + uphdil + uphdil2;
                                         }
                                         else{
                                               	ucom = ucom + uphdil;
                                         }                                           	
		                }//END if(bbsw==1)
			}//END if(disksw==1)
			else{
                               	ucom = 0.;
                        }//END else(disksw==1)

                        if(isVerbose){
                        	cout << left << setw(15) << "ucom:" << setw(15) << ucom << setw(15) << uphdil << setw(15) << uphdil2 << endl;
                        }

                        /**
                        * Find the maximum energy a particle can reach by solving:
                        * t^{-1}_acc = t^{-1}_syn + t^{-1}_compton + t^{-1}_esc
                        *
                        * This is equivalent to solve a equation du second degre en E:
                        * E**2 + escon/(syncon+comcon)*E - accon/(syncon+comcon)
                        *
                        */
                        accon	= 3./4. * fsc *cee *charg * bfield; // removed factors of cee here, since now fsc can be in sensible units!!
                        syncon	= 4./3. * sigtom * ub/(emgm*emgm*cee*cee*cee);
                        comcon	= syncon * ucom/ub;
                        escom	= beta * cee/z;
                        tdyn    = (brk * r) / (beta * cee);
                        qutrmb	= escom/(syncon+comcon);
                        qutrmc	= accon/(syncon+comcon);
                        emax	= (-qutrmb+sqrt(qutrmb*qutrmb+4.*qutrmc))/2.;
			ebreak  = 1. / (syncon * tdyn);                                     /* Calculation of synchrotron cooling break energy - do this by setting tsyn = tdyn, where electrons cooling time is equal to the time it takes them to escape the current jet segment. re-arrange this equation for the electron energy and you have the break energy */
                        gemax	= emax/emerg;
                        gebreak = ebreak/emerg;
                        //cout << "tdyn: " << tdyn << " syncon: " << 1./syncon << " emax: " << gemax << " ebreak: " << gebreak << endl;    
                        if(gemax >= 100){
                        	bemax	= 1. - 1./(2.*gemax*gemax) - 1./(8.*gemax*gemax*gemax*gemax);
                        }
                        else{
                        	bemax	= sqrt(gemax*gemax - 1)/gemax;
                        }
                            
                        /**
                        * From this emax can generate powerlaw between min and emax.
                        * Use emin normalized to scaled thermal peak at shock, and after shock
                        * then scale the two down together (see p.158, ntbk 6). This is a bit below
                        * the peak for lower temperatures
                        */
				
                        mjteff = eltemp*gshift;
			emin = 2.23*kboltz*mjteff; // changed from 2. to 2.23 ....this should be correct, right?
                        /**
                        * Create powerlaw of normalization plfrac*ntot - 3 options based on where the cooling break energy lies
                        */
                            
			if(ebreak <= emin) {
                        	cnorm = plfrac*ntot*(pspec)/(pow(emin, -pspec) - pow(emax, -pspec));  // if break energy less than minimum energy, set normalisation for cooled distribution
                        }
                        else if(ebreak >= emax){
                        	cnorm = plfrac*ntot*(1.-pspec)/(pow(emax, 1.-pspec) - pow(emin, 1.-pspec)); // if break greater than max energy, set normalisation for uncooled distribution
                        }
                        else{
                        	cnorm = plfrac*ntot/((pow(ebreak, 1.-pspec)/(1.-pspec)) - (pow(emin, 1.-pspec)/(1.-pspec)) + (pow(ebreak, 1.-pspec)/pspec) - (ebreak*pow(emax, -pspec)/pspec));
                                //   cout << left << setw(20) << "The break frequency is at approx.:   " << setw(15) << gebreak*gebreak*(charg*bfield/(emgm*cee)) << "Hz" << endl;
                                
                                // If break lies within energy limits, set normalisation for broken lower law
                        }
                          
                        /**
                        * If normalization of PL at emin greater than thermal then below emin,
                        * fix maxwellian which will require renormalizing again
                        */
                            
                        /**
                        * Make energy array between lelec[0] and log10(emax), add power-law
                        * and particle distribution together in appropriate ratios, then write
                        * into array and after shock, read it in and shift appropriately according
                        * to shock.
                        */
                        /**
                        * In theory this should be correctly normalized, but it's not.
                        * renorm to ntot at end to make more exact
                        */
				
                        gsl_interp_accel *acc_shock     = gsl_interp_accel_alloc();
                        gsl_spline *spline_shock        = gsl_spline_alloc(gsl_interp_cspline, nelec);
                        gsl_spline_init(spline_shock, lelec, ledens, nelec);
                            
                        einc = (log10(emax)-lelec[0])/nelec;
                         
                        //Below is full form, but extra factor of emerg cancels out
                        betat = emerg/(kboltz*mjteff);
                        mjpeak = pow(10.,9.514 - 1.028*log10(mjteff));
                            
                        if(ebreak <= emin) {
                        	enorm = cnorm*pow(emin,-(pspec+1.))*betat/(k2_fnc(betat)*mjpeak); // if break energy < min energy, normalisation at boundary between thermal/power law distributions is product of thermal and cooled power law
                        }
                        else{
                        	enorm = cnorm*pow(emin, -pspec)*betat/(k2_fnc(betat)*mjpeak); // if break energy > min energy, normalisation at boundary between thermal/power law distributions is product of thermal and uncooled power law
			}
                         
                        for(i=0; i<nelec; i++){
                        	etemp[i]= pow(10, lelec[0]+(i+0.5)*einc);
                                if(log10(etemp[i]) < lelec[nelec-1]){
                                    	edtrm = gsl_spline_eval(spline_shock, log10(etemp[i]), acc_shock);
                                    	edtrm = thmfrac*pow(10, edtrm);
                                }
                                else{
                                	edtrm = 0;
                                }
                                
                                thmshock[i] = edtrm; // thermal component assigned to an array at shock, so we can keep track of it
                                if(etemp[i] >= emin){
                                    	if(ebreak <= emin){
                                        	pltrm = cnorm*pow(etemp[i],-(pspec+1)); // if break energy < min energy, power law distribution completely cooled
                                	}
                                    	else if(ebreak >= emax){
                                        	pltrm = cnorm*pow(etemp[i],-pspec); // if break energy > max energy, power law distribution completely uncooled
                                    	}
                                    	else{  // if break energy within energy limits, set broken power law, power law index +1 after break energy
                                        	if(etemp[i] < ebreak){
                                            		pltrm = cnorm*pow(etemp[i],-pspec);					
                                        	}
                                       		else{
                                            		pltrm = cnorm*ebreak*pow(etemp[i],-(pspec+1));
                                        	}
                                    	}
                                }
				else{
                                	game = etemp[i]/emerg;
                               		if(game <= 1e3){
                                        	bete = sqrt(game*game-1.)/game;
                                    	}
                                    	else{
                                        	bete = 1. - 1./(2.*game*game) - 1.*(8*game*game*game*game);
                                    	}
                                    	pltrm = enorm * game*game*bete * exp(-game*betat); // was game/betat, which I think was wrong, changed 6 may 2015
                                }
                               // if(etemp[i] < 2.*emin/heat && pltrm > edtrm){ 	// This avoids a discontinuity between the power-law and thermal distributions, but in reality that discontinuity
                                	//pltrm = edtrm;				// does not impact the spectrum, so it's fine to ignore it. If it turns out to be important, it may crash the code
                               // }						// for high (>20) values of heat; in that case this condition can also be avoided by changing 2.*emin > 2.*emin/heat
                                
                                plshock[i] = pltrm; // power law component assigned to an array at shock, so can keep track of it
                                dtemp[i]= log10(edtrm+pltrm); // total distribution sum of edtrm and pltrm (thermal + power law), remembering to take logs
                                lelec[i]= log10(etemp[i]);
                                elen[i]	= etemp[i];
                                       
			} // end of loop for nelec
			gsl_spline_free(spline_shock), gsl_interp_accel_free(acc_shock); // free the gsl_spline
                            
			sum1    = 0.;
                        for(i=0; i < nelec-1; i++){
                               	tmp1 = (elen[i+1] - elen[i])*0.5*(pow(10, dtemp[i])+pow(10, dtemp[i+1])); // integrate distribution, and then renormalise ntot
                               	sum1 += tmp1;
                        }
                        renorm = ntot/sum1;
                            
                        if(isVerbose){
                        	cout << left << setw(15) << "renorm:" << setw(15) << ntot/sum1 << endl;
                        }
                            
                        shocknd	= ntot;
                        shockr	= r;
                        shockgb	= rvel;
                            
                        for(i=0; i<nelec; i++){
                        	ledens[i]= log10(renorm*pow(10,dtemp[i]));
				eled[i]	= elen[i]*pow(10,ledens[i]);
                                shelen[i]= lelec[i];
                                shden[i]= dtemp[i];
                                plshock[i] = renorm*plshock[i]; // renormalise the power law component at shock, to be read in after shock
                                thmshock[i] = renorm*thmshock[i]; // renormalise the thermal component at shock, to be read in after shock
			}
                            
                        sum1    = 0.;
                        sum2    = 0.;
                        for(i=0; i<nelec-1; i++){
                        	tmp1	= (elen[i+1] - elen[i])*0.5*(eled[i]+eled[i+1]);
                                tmp2	= (elen[i+1] - elen[i])*0.5*(pow(10, ledens[i]) + pow(10, ledens[i+1]));
                                sum1	+= tmp1;
                                sum2	+= tmp2;
			}
                            
                        if(isVerbose){
                        	cout << "EQUIP at Shock: " << bfield*bfield/(8.*pi*sum1) << endl;
				cout << "NTOT after Shock: " << sum2 << " and " << ntot << endl;
			}
                            
			// <!-- FIXME insert proton acceleration here -->
                            
			isShock = true;
                             
			}//END if z>zsh
		}//END nozzle + shock
                else if (z >= zsh && isShock){
                 
                /****************
                * BEYOND SHOCK *
                ****************/
                    
                /* Here we want to read in the distribution as it was at the shock, and then read in the thermal/power law parts of the distribution separately (thmcomp[i] and plcomp[i]). Then we essentially repaeat the process completed at the shock, calculating the new min/max/break energies of the power law distribution according to the now gshifted thermal distribution. In this way we re-calculate the power law distribution as if it were fully replenished from the thermal distribution  */
                          
                    
			//Can read stored distribution directly now
                    	for(i=0; i < nelec; i++){
                        	etemp[i]= pow(10.,shelen[i]);
                        	dtemp[i]= shden[i];
                        	if(thmshock[i] > 0){ // after elmax (initial max energy of thermal distribution), thermal term = 0, thus cannot take logs....only take logs of values > 0, set all others = 0.
                            	thmcomp[i] = log10(thmshock[i]);
                        	}
                        	else{
                            	thmcomp[i] = 0;
                        	}
                        	plcomp[i] = log10(plshock[i]); // load in power law component from shock calculation
                    	}
                    
                    
                    	//Shift number density down from shock, not form ntot0 anymore
                    	ntot	= shocknd*pow(shockr/r,2)*(shockgb/rvel);
               
			//Now shift down in energies and lose the lost guys.
                    	//NOTE: the shift is from the shock = gshift/gshock
                    
                    
                    	/****************************************************************************************************************************/
                    	/*** Repeat the process completed at shock for every jet segment (determine ub/ucom, max energy, break energy, and shift) ***/
                    	/****************************************************************************************************************************/
                                       
                    	if(isVerbose){
                        	cout << left << setw(15) << "ACCELPARS:" << setw(15) << r << setw(15) << z << setw(15) << dopfac[0*nz+k] << setw(15) << dopfac[1*nz+k] << setw(15) << bfield << setw(15) << ntot << endl;
                    	}
                    
                    	ub      = bfield*bfield/(8.*pi);
                    
                    	sum1    = 0.;
                    	for(m=0; m < njet; m++){
                        	for(i=0; i < nsyn; i++){
                            		phofrq[i]= pow(10.,nusyn[(m*nz+(k-1))*nsyn+i])/dopfac[m*nz +(k-1)];
                            		//3 powers of dopfac because one was inside synint for angle aberration
                            		phoint[i]= mjy * pow(10., synabs[(m*nz+(k-1))*nsyn+i])/pow(dopfac[m*nz +(k-1)],3);
                            	//Handle differently: not flux but assuming that the synchrotron
                            	//was produced in this segment. That's overestimating but hopefully not
                            	//too much because r, delz also increased.
                            	phoint[i]*= 4. * dist*dist/(r*r*cee);
                        	}
                        	//Rough integration to get approx energy density
                        	for(i=0; i< nsyn-1; i++){
                            		tmp1	= (phofrq[i+1]-phofrq[i])*0.5*(phoint[i]+phoint[i+1]);
                            		sum1	+= tmp1;
                        	}
                    	}
                    	ucom    = sum1;

                    	if(isVerbose){
                        	cout << left << setw(15) << "ucom:" << setw(15) << ucom << endl;
                    	}
                    
                    	/**
                     	* Add in multicolor disk contribution
                     	* (old versions: plus component from irradiation, w/ effective area pi*reff**2)
                     	*/
			if(disksw==1){
				bbdil = 2.*reff2/(2.*reff2+pow(z-hbb/2.,2))-reff2/(reff2+pow(z-hbb/2.,2));
				if(bbsw==1){
                                	nubot = 1.5e-5*kboltz*tin/herg;
                                        nutop = 8.*kboltz*tin/herg;
                                        if(bbdisk(log(nutop), gamv, tin, rin, rout, z, hbb)*nubot/(bbdisk(log(nubot), gamv, tin, rin, rout, z, hbb)*nutop) > 1.e30){
                                        	cerr << "nutop not high enough" << endl;
                                                exit(1);
                                     	}
                                        tst	= bbdisk(log(nutop), gamv, tin, rin, rout, z, hbb);
                                        if(tst <= 1.e-20){
                                        	while(tst < 1.e-20){
                                                	nutop *= 0.75;
                                                        tst = bbdisk(log(nutop), gamv, tin, rin, rout, z, hbb);
                                                }
                                    	}
                                        bbintegrals(log(nubot), log(nutop), gamv, tin, rin, rout, z, hbb, bbint);
				    	uphdil = bbint/cee;
				    	uphdil2	= bbdil*aconst*tbbeff*tbbeff*tbbeff*tbbeff/4.;
   					if(z > hbb/2.){
                                               	ucom = ucom + uphdil + uphdil2;
                                         }
                                         else{
                                               	ucom = ucom + uphdil;
                                         }                                           	
		                }//END if(bbsw==1)
			}//END if(disksw==1)
			else{
                               	ucom = 0.;
                        }//END else(disksw==1)
                    
                    	if(isVerbose){
                    		cout << left << setw(15) << "ucom:" << setw(15) << ucom << setw(15) << uphdil << setw(15) << uphdil2 << endl;
                    	}
                    
                    	/**
                     	* Find the maximum energy a particle can reach by solving:
                     	* t^{-1}_acc = t^{-1}_syn + t^{-1}_compton + t^{-1}_esc
                     	*
                     	* This is equivalent to solve a equation du second degre en E:
                     	* E**2 + escon/(syncon+comcon)*E - accon/(syncon+comcon)
                     	*
                     	*/
                    	accon	= 3./4. * fsc * cee* charg * bfield;
                    	syncon	= 4./3. * sigtom * ub/(emgm*emgm*cee*cee*cee);
                    	comcon	= syncon * ucom/ub;
                    	escom	= beta * cee/z;
                    	tdyn    = (brk * r) / (beta * cee);
                    	qutrmb	= escom/(syncon+comcon);
                    	qutrmc	= accon/(syncon+comcon);
                    	emax	= (-qutrmb+sqrt(qutrmb*qutrmb+4.*qutrmc))/2.;
                    	ebreak  = 1. / (syncon * tdyn);                                     /* Calculation of synchrotron cooling break energy - do this by setting tsyn = tdyn, where electrons cooling time is equal to the time it takes them to escape the current jet segment. re-arrange this equation for the electron energy and you have the break energy */
                    	double tacc;
                    	tacc = syncon + comcon + escom;
                    
//                    cout << "t_acc: " << 1./tacc << " and kTe/mc2: " << kboltz*eltemp/(emgm*cee*cee) << " and kTe: " << eltemp*kboltz << endl;
                    
                    	gemax	= emax/emerg;
                    	gebreak = ebreak/emerg;
                    	//cout << "tdyn: " << tdyn << " syncon: " << 1./syncon << " emax: " << gemax << " ebreak: " << gebreak << endl;
                    	if(gemax >= 100){
                        	bemax	= 1. - 1./(2.*gemax*gemax) - 1./(8.*gemax*gemax*gemax*gemax);
                    	}
                    	else{
                       		bemax	= sqrt(gemax*gemax - 1)/gemax;
                    	}
                    
                    	/**
                     	* From this emax can generate powerlaw between min and emax.
                     	* Use emin normalized to scaled thermal peak at shock, and after shock
                     	* then scale the two down together (see p.158, ntbk 6). This is a bit below
                     	* the peak for lower temperatures
                     	*/
                    	mjteff	= eltemp*gshift;
                    	emin	= 2.23*kboltz*mjteff; // changed from 2. to 2.23 ....this should be correct, right?
                    	//cout << "gamin: " << emin/emerg << " gabreak: " << gebreak << " gamax: " << gemax << endl; 
                    
                    	/* Calculate power law component according to mjteff and energy limits in this jet segment (including the energy break) */
                    
                    	if(ebreak <= emin) {
                        	cnorm = plfrac*ntot*(pspec)/(pow(emin, -pspec) - pow(emax, -pspec)); // if break energy less than minimum energy, set normalisation for cooled distribution
                        	//cout << "break energy less than min energy\n";
                    	}
                    	else if(ebreak >= emax){
                        	cnorm	= plfrac*ntot*(1.-pspec)/(pow(emax, 1.-pspec) - pow(emin, 1.-pspec)); // if break greater than max energy, set normalisation for uncooled distribution
                        	// cout << "break energy greater than max energy\n";
                    	}
                    	else{
                        	cnorm = plfrac*ntot/((pow(ebreak, 1.-pspec)/(1.-pspec)) - (pow(emin, 1.-pspec)/(1.-pspec)) + (pow(ebreak, 1.-pspec)/pspec) - (ebreak*pow(emax, -pspec)/pspec));
                        	//   cout << left << setw(20) << "The break frequency is at approx.:   " << setw(15) << gebreak*gebreak*(charg*bfield/(emgm*cee)) << "Hz" << endl;
                        	// if break energy > min energy, normalisation at boundary between thermal/power law distributions is product of thermal and uncooled power law
                        	// cout << "Break present in energy limits\n";
                    	}
                    
                    	/**
                     	* If normalization of PL at emain greater than thermal then below emin,
                     	* fix maxwellian which will require renormalizing again
                     	*/
                    
                    	/**
                     	* We do not reset the energy array here, despite re-calculation of the energy limits, because emax only increases, and the electron densities at higher energies are increasingly more negligible. Add power-law
                     	* and particle distribution together in appropriate ratios, then write
                     	* into array and after shock, read it in and shift appropriately according
                     	* to shock.
                     	*/
                    	/**
                     	* In theory this should be correctly normalized, but it's not.
                     	* renorm to ntot at end to make more exact
                     	*/
                    	//Below is full form, but extra factor of emerg cancels out
                    	betat	= emerg/(kboltz*mjteff);
                    	mjpeak	= pow(10.,9.514 - 1.028*log10(mjteff));
                    
                    	if(ebreak <= emin) {
                        	enorm	= cnorm*pow(emin,-(pspec+1))*betat/(k2_fnc(betat)*mjpeak); // if break energy < min energy, normalisation at boundary between thermal/power law distributions is product of thermal and cooled power law
                    	}
                    	else{
                        	enorm = cnorm*pow(emin, -pspec)*betat/(k2_fnc(betat)*mjpeak); // if break energy > min energy, normalisation at boundary between thermal/power law distributions is product of thermal and uncooled power law
		       }
                    
                    	for(i=0; i<nelec; i++){
                        
                        	if(etemp[i] >= emin){
                            		if(ebreak <= emin){
                                		pltrm	= cnorm*pow(etemp[i],-(pspec+1)); // if break energy < min energy, power law distribution completely cooled
                            		}
                            		else if(ebreak >= emax){
                                		pltrm = cnorm*pow(etemp[i],-pspec); // if break energy > emax, power law distribution completely uncooled
                            		}
                            		else{
                                		if(etemp[i] < ebreak){
                                    			pltrm = cnorm*pow(etemp[i],-pspec); // if break energy inside energy limits, broken power law distribution with break at ebreak
                                		}
                                		else{
                                    			pltrm = cnorm*ebreak*pow(etemp[i],-(pspec+1));
                                		}
                            		}
                        	}
                        	else{
                            		game	= etemp[i]/emerg;
                            		if(game <= 1e3){
                                		bete	= sqrt(game*game-1.)/game;
                            		}
                            		else{
                                		bete	= 1. - 1./(2.*game*game) - 1.*(8*game*game*game*game);
                            		}
                            			pltrm	= enorm * game*game*bete * exp(-game*betat); // was game/betat, which I think was wrong, changed 6 may 2015
                        	}
                                                
                        	plcomp[i]= log10(pltrm);
                    	}
                    
                    /* Now shift down the thermal component (thmcomp), determine the power law component, sum them, and renormalise everything accordingly  - this differs from previous version, in that before we were shifting down the full distribution, wheras now we shift down the thermal component only, and recalculate the power law component according to the thermal peak of the adiabatically cooled thermal distribution. In this way the break can be re-calculated at each step in a simpler way.  */
                    
                    	nw	= 0;
                    	for(i=0; i<nelec; i++){
                        	etemp[i] = log10(etemp[i]);
                        	etemp[i]= log10(gshift/gshock*pow(10,etemp[i]));
                        	if(pow(10,etemp[i])/emerg <= 1.){
                            		nw++;
                        	}
                    	}
                    
                    	/* calculate newly shifted down thermal component */
                    	sum1	= 0;
                    	for(i=0; i<nelec-1; i++){
                        
                        	tmp1	= (pow(10,etemp[i+1]) - pow(10,etemp[i]))*0.5*(pow(10,thmcomp[i]) + pow(10,thmcomp[i+1]));
                        	sum1	+= tmp1;
                        
                    	}
                    	renorm	= log10(ntot*thmfrac/sum1);
                    	for(i=0; i<nelec; i++){
                        	if(thmcomp[i]>0){ // only read in the thermal component (and take log value) if it is non-zero. Else thmcomp[i] = 0.
                            		thmcomp[i]+= renorm;
                        	}
                        	else{
                            		thmcomp[i]=0;
                        	}
                        	if(pow(10.,etemp[i])*(gshock/gshift) < 2.*emin && plcomp[i]>thmcomp[i]){ // to smooth out distirbution, for a power law component > thermal at e < 2*emin, power law component = thermal component.
                            		plcomp[i]=thmcomp[i];
                        	}
                        	if(thmcomp[i]>0){ // again, thermal component should only be included in total sum if it is non-zero
                            		dtemp[i] = log10(pow(10.,thmcomp[i])+pow(10.,plcomp[i]));
                        	}
                        	else{
                            		dtemp[i] = log10(pow(10.,plcomp[i]));
                        	}
			}
                    
                    
                    	//spline_dbl and resize the array into the old lelec, ledens if you lose particles at low energy
                    	if(nw>0){
                        	if(nw==1){
                            		totlos	= (pow(10,etemp[1])-pow(10,etemp[0]))*pow(10,dtemp[0]);
                        	}
                        	else if(nw>1){
                            		totlos	= 0;
                            		for(i=0; i<nw-1; i++){
                                		tmlos	= (pow(10,etemp[i+1])-pow(10,etemp[i]))*0.5*(pow(10,dtemp[i])+pow(10,dtemp[i+1]));
                                		totlos	+= tmlos;
                                }
                        }//END else if(nw>1)
                        
                        oldnum  = ntot;
                        ntot    -= totlos;
                        renorm  = log10(ntot/oldnum);
                        for(i=0; i<nelec; i++){
                            	dtemp[i]+= renorm;
                        }
                        
                        gsl_interp_accel *acc_shift  = gsl_interp_accel_alloc();
                        gsl_spline *spline_shift     = gsl_spline_alloc(gsl_interp_cspline, nelec);
                        gsl_spline_init(spline_shift, etemp, dtemp, nelec);
			
                        einc    = (etemp[nelec-1]-etemp[nw])/nelec;
			
                        for(i=0; i<nelec; i++){
                            	lelec[i]= etemp[nw]+(i+0.5)*einc;
                            	ledens[i] = gsl_spline_eval(spline_shift, lelec[i], acc_shift);
                            	elen[i] = pow(10., lelec[i]);
                            	eled[i] = elen[i]*pow(10., dtemp[i]);
                        }
			
                        gsl_spline_free(spline_shift), gsl_interp_accel_free(acc_shift);
                    	}//END if(nw>0)
                    	else{
                    		for(i=0; i<nelec; i++){
                            		lelec[i]= etemp[i];					
                            		elen[i] = pow(10.,lelec[i]);
                            		eled[i] = elen[i]*pow(10.,dtemp[i]);
                            		ledens[i]= dtemp[i];
                            	}
                    	}//END else nw > 0
                    
                    
                    	sum1 = 0.; // ADDED 04/04/16, renormalise electron distribution in each segment, such that ratio ntot/sum1=1 post-normalisation.
                    	for(i=0; i<nelec-1; i++) {
                        	tmp1 = (pow(10.,lelec[i+1]) - pow(10.,lelec[i])) * 0.5 * (pow(10.,ledens[i]) + pow(10.,ledens[i+1]));
                        	sum1 += tmp1;
                    	}
                    	renorm	= log10(ntot/sum1);
                    	for(i=0; i<nelec; i++){
                        	ledens[i] += renorm;
                    	}

                    	if(isVerbose){
                    		/* Double check equipartition */
                        	sum1    = 0.;
                        	for(i=0; i<nelec-1; i++){
                            		tmp1 = (elen[i+1]-elen[i])*0.5*(eled[i]+eled[i+1]);
                            		sum1 +=tmp1;
                        	}
                        	ub=bfield*bfield/(8.*pi);
                        	cout << "CHECK equipartition: " << ub/sum1 << endl;
                    	}
                    
                }//END elseif (z>=zsh && isShock)
                else{
                	cerr << "***Error no correct isShock and zsh values";
			cout << " z: " << z/r_g << " k : " << k << " zsh: " << zsh/r_g << endl;
                    	exit(1);
                }
		
                /***************************
                 * NOZZLE + SHOCK + BEYOND *
                 *       RADIATION         *
                 ***************************/
                //Back to common actions to all jet
                
                //<!-- FIXME insert here hadronic routines -->
                sum1    = 0.;
                sum2    = 0.;
                for(i=0; i<nelec-1; i++){
                	tmp1	= (elen[i+1]-elen[i])*0.5*(eled[i]+eled[i+1]);
                    	tmp2	= (elen[i+1]-elen[i])*0.5*(pow(10,ledens[i+1])+pow(10,ledens[i]));
                    	sum1	+=tmp1;
                    	sum2	+=tmp2;
                }
                avelec	= sum1/sum2;
            
//              cout << "just before running out of particles: Ntot= " << ntot << endl;
                if(ntot <= 0){
                	for(i=0; i<nsyn; i++){
                        	for(m=0; m<njet; m++){
                            		nusyn[(m*nz+k)*nsyn+i]= log10(nurad[i]*dopfac[m*nz+k]);
                            		synabs[(m*nz+k)*nsyn+i]= -100;
                        	}
                   	}
                    	cerr << "ERROR: ***Ran out of particles (not relativistic anymore!)" << endl;
                    	exit(1);
                }
		
                gsl_interp_accel *acc_eldis1    = gsl_interp_accel_alloc();
                gsl_spline *spline_eldis1       = gsl_spline_alloc(gsl_interp_cspline, nelec);
                gsl_spline_init(spline_eldis1, lelec, ledens, nelec);
		
                elenmn   = elen[0];
                elenmx   = elen[nelec-1];
                for(i=0; i<nelec; i++){
                    	tstrm[i]= pow(10,ledens[i])/(elen[i]*elen[i]);
                }
                for(i=0; i<nelec-1; i++){
                    	drtrm[i]= (tstrm[i+1]-tstrm[i])/(elen[i+1]-elen[i]);
                }
                drtrm[nelec-1]	= drtrm[nelec-2];
		
                gsl_interp_accel *acc_derivs    = gsl_interp_accel_alloc();
                gsl_spline *spline_derivs       = gsl_spline_alloc(gsl_interp_cspline, nelec);
                gsl_spline_init(spline_derivs, elen, drtrm, nelec);
		
                /**
                 * 4.\times \pi in synabs is for assumed isotropic source.
                 * synint returns units of ergs/cm^2/s/Hz for absorbed
                 * convert to Jansky=1.e-23 erg/cm^2/s/Hz
                 */
                for(i=0; i<nsyn; i++){
                    	for(m=0;m<njet; m++){
                        	nusyn[(m*nz+k)*nsyn+i]= log10(nurad[i]*dopfac[m*nz+k]);
                    	}
                }

                isBreaknjetnsyn = false;
                for(i=0; i<nsyn; i++){
                        /**
                         * First perform integrations for synchrotron emissivity
                         * and absorption at current frequency.
                         */
			synintegrals(nurad[i], bfield, elenmn, elenmx, spline_syn, acc_syn, spline_eldis1, acc_eldis1, spline_derivs, acc_derivs, esum, asum);

                        for(m=0; m<njet; m++){
                                synint(nurad[i], r, inclin, dopfac[m*nz+k], delz, bfield, dist, esum, asum, absd, phoden);

                                nphot[i]= phoden;
//                      	cout << nurad[i] << " " << delz << " " << esum << " " << asum << " " << absd << " " << phoden << endl;
                            
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
            
//            	cout << "I have commented zcut here!\n";
                if(z>zcut){
                    	gsl_spline_free(spline_eldis1), gsl_interp_accel_free(acc_eldis1);
                    	gsl_spline_free(spline_derivs), gsl_interp_accel_free(acc_derivs);
                    	continue; //start new iteration of BIG LOOP
                }
                /**
                 * spline_dbl photon density (erg/s/Hz) for passing to Compton routines,
                 * divide by herg*cee*energy(erg)*pi*r**2 (erg**2*s**2*cm**3*Hz*s**-1)
                 * to get #/cm**3/erg
                 * -Add in the diluted blackbody from the disk
                 *  -Make the max energy based on tin
                 */
                for(i=0; i<nsyn; i++){
//                  	cout << nurad[i] << " " << nphot[i] << endl;
                    	if(nphot[i] == 0){
                        	nphot[i]= nphot[i-1]-4.;
                    	}
                    	else{
                        	nphot[i]= log10(nphot[i]);
                    	}
                    	jetu[i]	= pow(10,nphot[i])/(cee*herg*energ[i]*pi*r*r);
//                    		cout << i << " " << nphot[i]<< " " << jetu[i] << " " << cee*herg*energ[i]*pi*r*r << endl;
                }

                if(isVerbose){
                        /**
                         * Check SSC energy density to compare to ucom1 at shock
                         */
                        sum1	= 0;
                        for(i=0; i<nsyn-1; i++){
                                tmp1	= (pow(10, nubb[i+1])-pow(10,nubb[i]))*0.5*(pow(10,nphot[i+1])+pow(10,nphot[i]));
                                sum1	+=tmp1;
                        }
                        toten	= sum1/(cee*pi*r*r);
                        cout << "ENDEN FOR SSC: " << toten << endl;
                }
                
                if(disksw==1){
                        bbnumax	= log10(1e2*tin*kboltz/herg);
                        if(bbnumax > snumax){
				cerr << "Need to change BB+syn for IC" << endl;
				cerr << "BB going to higher freq than syn" << endl;
				exit(1);
			}
			bbdil	= 2.*reff2/(2.*reff2+(z-hbb/2.)*(z-hbb/2.))-reff2/(reff2+(z-hbb/2.)*(z-hbb/2.));
			
                        /**
			 * Stop Compton loop when too far past photon peak
			 */
			photmx	= 1.e-20;
			peaksw	= 0;
			bbcon	= 0;

                        for(i=0; i<nsyn; i++){
                                if(bbsw==1){
                                        if (peaksw!=1){
                                                if(nubb[i] <= bbnumax){
                                                        if(photmx > 1.e-20 && bbcon/photmx <= 1.e-13){
                                                                if(isVerbose){
                                                                        cout << left << setw(15) << "gone too low" << setw(15) << i << setw(15) << photmx << setw(15) << bbcon << endl;
                                                                }
                                                            	bbcon	= 0.;
                                                            	peaksw	= 1.;
                                                        }
                                                        else{
                                                            	bbcon	= bbdisk(log(nurad[i]), gamv, tin, rin, rout, z, hbb)/(cee*herg*energ[i]*nurad[i]);
                                                            	photmx	= max(photmx,bbcon);
                                                        }//end else(photmx)
                                                }//END if(nubb[i] <= bbnumax)
                                                else{
                                                	bbcon   = 0.;
                                                }//END nubb
                                        	bbirad	= bbdil*2.*herg*pow(nurad[i],3)/(cee*cee*(exp(energ[i]/(kboltz*tbbeff))-1.)*cee*herg*energ[i]);        
					}//END if(peak!=1)

					if(z > hbb/2){
                                               	phodis[i]= log10(jetu[i]+bbcon+bbirad);
                                               	if (jetu[i]+bbcon == 0) phodis[i]= -500;
                                               	disku[i]= bbcon + bbirad;
                                        }
                                        else{
                                              	phodis[i]= log10(jetu[i]+bbcon);
                                                if (jetu[i]+bbcon == 0) phodis[i]= -500;
                                                disku[i]= bbcon;
                                        }
                                }//END if(bbsw==1)
                                else{
                                    	phodis[i]= log10(jetu[i]);
                                    	if (jetu[i] == 0) phodis[i]= -500;
                                }//end bbsw
                        }//END for-loop nsyn
                }//END if(disksw==1)
                else{
                    	for(i=0; i<nsyn; i++){
                        	phodis[i]= log10(jetu[i]);
                        	if (jetu[i] == 0) phodis[i]= -500;
                    	}
		}
                
        	/**
         	* Check energy densities
         	*/
        	sum1	= 0;
        	sum2	= 0;
    
//        	cout << "\n\nFirst Source Term for Compton: \n";
        	for(i=0; i<nsyn-1; i++){
                	tmp1	= (pow(10,phodis[i+1]+ephot[i+1])+pow(10,phodis[i]+ephot[i]))*0.5*(pow(10,ephot[i+1])-pow(10,ephot[i]));
                	tmp2	= (pow(10,phodis[i+1])+pow(10,phodis[i]))*0.5*(pow(10,ephot[i+1])-pow(10,ephot[i]));
                	sum1	+= tmp1;
                	sum2	+= tmp2;
//                	cout <<  i << " " << ephot[i] << " " <<  phodis[i] << endl;
        	}

		avphot	= sum1/sum2;
//        	cout << "pieces: avphot: " << avphot << ", sum1: " << sum1 << ", sum2: " << sum2 << ", ub: " << b0*b0/(4*pi) << endl;
//        	cout << "First source comtpon ucom/ub: " << sum1/(b0*b0/(4*pi)) << endl;

        	ephmax	= pow(10,ephot[nsyn-1]);
        	ephmin	= pow(10,ephot[0]);

        	if(isVerbose){
        		cout << "Av photon: " << avphot << " and max/min energies: " << ephmax << "/" << ephmin << endl;
        	}

        	gsl_interp_accel *acc_com2      = gsl_interp_accel_alloc();
        	gsl_spline *spline_com2         = gsl_spline_alloc(gsl_interp_cspline, nsyn);
        	gsl_spline_init(spline_com2, ephot, phodis, nsyn);
	
        /**
		 * Call integrate Compton function over electrons, returns #/cm**3/s/erg from #/cm**3/erg
		 */
        	double sumcom = 0.0, sumcom_prev = 0.0, tmpcom = 0.0, com_prev = 0.0;
        	double sumpp = 0.0, tmppp = 0.0;
        	double yParR,yParNR,yPar,tau;
        	int intYpar;
        	tau = r*ntot*sigtom;
		if (kboltz*eltemp/emerg < 1.0) {
            		yParNR = pow(4.0*kboltz*eltemp/emerg,1)*max(tau*tau,tau);
            		yPar = yParNR;
        	} else {
            		yParR = pow(4.0*kboltz*eltemp/emerg,2)*max(tau*tau,tau);
        	    	yPar = yParR;
        	}
        	intYpar = int(yPar+0.5);
//            	cout << int(yPar) << " " << yPar << endl;
//        	cout << "Starting Compton-Y: " << yPar << endl;
            
        	/* Cleaning up the array for new source term in between slices */
        	for(i=0; i<ncom; i++){
            		comi[i] = 1e-100;
        	}
    
//            	Niter = 0;
        	for (int it=0; it<Niter;it++){
            	/* Cleaning up stuff before filling them */
            		for(i=0; i<ncom; i++){
                		nucom_src[i] = -100.;
                		comspc_src[i]= -100.;
//                		if (k==1) cout << nucom_src[i] << " " << comspc_src[i] << endl;
            		}
            
            		/* New Source Terms for Higher Orders */
            		if (it > 0) {
//              		cout << "\n\n" << it+1 << "-th source for Compton\n\n";
                		for(i=0; i<ncom; i++){
                    			if (comi[i]==0) comi[i] = 1e-100;
                    			nucom_src[i] = log10(ephxr[i]*herg);
                    			comspc_src[i]= log10(comi[i])-log10(herg*herg*cee*ephxr[i]*pi*r*r);
                    
                    			/* The energy limits must be changed in order to avoid problems with interpolation */
                    			ephmin = ephxr[0]*herg;
                    			ephmax = ephxr[ncom-1]*herg;
//                    	if (k==1) cout << nucom_src[i] << " " << comspc_src[i] << endl;
                		}
            		} else {
                	/* This is a dummy filling for the it=0 because it's not needed */
                		for(i=0; i<ncom; i++){
                    			nucom_src[i] = log10(ephxr[i]*herg);
                    			comspc_src[i]= -100.;
                		}
            		}
			
			gsl_interp_accel *acc_comiter      = gsl_interp_accel_alloc();
	            	gsl_spline *spline_comiter         = gsl_spline_alloc(gsl_interp_akima, ncom);
	            	gsl_spline_init(spline_comiter, nucom_src, comspc_src, ncom);

// 	           	cout << "\n\n";
            
//      	      	cout << "This is the Compton spectrum for it= " << it << endl;
        	    	sumcom_prev = sumcom;
        	    	sumcom = 0.0;
        	    	for(i=0; i<ncom; i++){
        	        	//ephxr is in Hz remember!!
        	        	eph     = ephxr[i]*herg;

        	        	for(m=0; m<njet; m++){
        	            		nucom[(m*nz+k)*ncom+i]= log10(eph*dopfac[m*nz+k]/herg);
				}
        	        	//avoid out of bounds error
        	        	if(i >= 2 && comspc[(0*nz+k)*ncom+(i-1)] <= -20. && (comspc[(0*nz+k)*ncom+(i-1)]-comspc[(0*nz+k)*ncom+(i-2)]) <= 0. ){
        	            		for(m=0; m<njet; m++){
        	                		comspc[(m*nz+k)*ncom+i]= comspc[(m*nz+k)*ncom+(i-1)]-4.;
        	            		}
        	            	continue;//Out of i(ncom) for loop for a new iteration
        	        	}
        	        	/* Is this like hnu*gamma**2? the max energy gain for an average photon? */
        	        	//typmax	= avphot*(avelec/emerg)*(avelec/emerg);
	
        	        	blim	= max(log(elenmn/emerg),log(eph/emerg))+1.e-15;
                
        	        	if(blim >= log(elenmx/emerg)){
        	            		com	= 1e-100;
//      	              		cout << "I am in here for i= " << i << "\n";
        	        	}
        	        	else{
        	            		if (it==0) {	
						com = comintegral(blim, log(elenmx/emerg), eph, ephmin, ephmax, spline_eldis1, acc_eldis1, spline_com2, acc_com2);
					} else {
						com_prev = com;
        	                		com = comintegral(blim, log(elenmx/emerg), eph, ephmin, ephmax, spline_eldis1, acc_eldis1, spline_comiter, acc_comiter);; 
        	            		}
//      	                  	cout << "Interpolation error at i= " << i << endl;
        	        	}
        	        	comi[i] = com*ephxr[i]*herg*herg*vol;
        	        	if (com == 0) comi[i] = 1e-100;
//      	          	cout << i << " " << log10(eph) << " " << log10(comi[i]) << endl;
                
                		/**
                 		* Convert to ergs/cm**2/s/Hz and the mJy
                 		* Write in Hz and mJy
                 		*/
                		for(m=0; m<njet; m++){
                    			comspc[(m*nz+k)*ncom+i]=com*ephxr[i]*herg*herg*vol*dopfac[m*nz+k]*dopfac[m*nz+k]/(mjy*4.*pi*dist*dist);
                    			comspc_it[(m*nz+k)*ncom+i] = comspc_it[(m*nz+k)*ncom+i] +  comspc[(m*nz+k)*ncom+i];
//                    			if (m==0 && (k==0 || k==5)) cout << i << " " << log10(eph) << " " << comspc_it[(m*nz+k)*ncom+i] << endl;
                    			comspc[(m*nz+k)*ncom+i]=comspc_it[(m*nz+k)*ncom+i];
					if(comspc[(m*nz+k)*ncom+i] == 0){
                        			comspc[(m*nz+k)*ncom+i]=comspc[(m*nz+k)*ncom+(i-1)]-2.;
                    			}
                    			else{
                        			comspc[(m*nz+k)*ncom+i]=log10(comspc[(m*nz+k)*ncom+i]);
                    			}
                		}
	                	/* Integral to compare contributions of Compton orders */
	                	tmpcom	= (com_prev+com)*0.5*(ephxr[i+1]*herg-ephxr[i]*herg);
	                	sumcom	+= tmpcom;
                
	                	/* Integral over enegies to evaluate the photons 
	                 	* which are contributing to pair production.
	                 	*
	                 	*/
	               		sumpp = 0.0;
	                	for(int ii=0; ii<ncom; ii++) {
                    			if (eph*ephxr[ii]*herg < 2.0*emerg*emerg) {
                        			if (ii == ncom-1){
                            				tmppp = (com_prev+com)*ephxr[ii]*herg;
                        			} else {
                            				tmppp = (com_prev+com)*0.5*(ephxr[ii+1]*herg-ephxr[ii]*herg);
                        			}
                        		sumpp += tmppp;
//					cout << i << " " << tmppp << " " << ephxr[ii+1]*herg << " " << ephxr[ii]*herg << " "  << (ephxr[ii+1]*herg-ephxr[ii]*herg) << endl;
                    			}
                		}
            		}//END for(ncom)
       
            		if ((it == intYpar) || (it > 8)) {/*
                		if (kboltz*eltemp/emerg < 1.0) {
                    			cout << "Final N-scatt: " << it+1 << ", Compton-Y par NR: " << yPar << endl;
                		} else {
                    			cout << "Final N-scatt: " << it+1 << ", Compton-Y par R: " << yPar << endl;
                		}
                		if (it == 100) {
                    			cout << "Too many orders!!! Compton-Y par is too big: " << yPar << endl;
                		}*/
				gsl_spline_free(spline_comiter), gsl_interp_accel_free(acc_comiter);
				break;
            		}
            		gsl_spline_free(spline_comiter), gsl_interp_accel_free(acc_comiter);
		} // END Niter
        	gsl_spline_free(spline_com2), gsl_interp_accel_free(acc_com2);
//        	gsl_spline_free(spline_comiter), gsl_interp_accel_free(acc_comiter);
        	gsl_spline_free(spline_eldis1), gsl_interp_accel_free(acc_eldis1);
        	gsl_spline_free(spline_derivs), gsl_interp_accel_free(acc_derivs);
	}//END of BIG LOOP OVER Z
//    cout << "THE end of the loop over z:" << pow(10,zmax) << " " << z << " "<< k << endl;

    /**
	 * spline_dbl arrays for adding to big total array for z<zmax, and add BB flux, nutot in Hz
	 */
	
	for(m=0; m<njet; m++){
		for(k=0; k<nz; k++){
			z	= log10(zed[k]);
			for(i=0; i<nsyn; i++){
				snu[i]	= nusyn[(m*nz+k)*nsyn+i];
				sdump[i]= synabs[(m*nz+k)*nsyn+i];			
			}
          	
            	gsl_interp_accel *acc_snu       = gsl_interp_accel_alloc();
            	gsl_spline *spline_snu          = gsl_spline_alloc(gsl_interp_cspline, nsyn);
            	gsl_spline_init(spline_snu, snu, sdump, nsyn);
       
      		gsl_interp_accel *acc_cnu       = gsl_interp_accel_alloc();
            	gsl_spline *spline_cnu          = gsl_spline_alloc(gsl_interp_akima, ncom);


		if(z<= log10(zcut)){
			for(i=0; i<ncom; i++){
				cnu[i]	= nucom[(m*nz+k)*ncom+i];
				cdump[i]= comspc[(m*nz+k)*ncom+i];				
			}
			
	                gsl_spline_init(spline_cnu, cnu, cdump, ncom);
		}
			
			for(i=0; i<(ne-1); i++){
				if(nutot[i] >= snu[0] && nutot[i] <= snu[nsyn-1]){
		        	        sflx    = gsl_spline_eval(spline_snu, nutot[i], acc_snu);
				}
				else{
					sflx	= -200.;
				}	
				if(z> log10(zcut)){
					cflx	= -200.;
				}
				else{
					if(nutot[i] >= cnu[0] && nutot[i] <= cnu[ncom-1]){
						cflx    = gsl_spline_eval(spline_cnu, nutot[i], acc_cnu);
					} 
					else{
						cflx	= -200.;
					}
				}

				if(plotsw == 1){
					complot[i]= complot[i] + pow(10,cflx);
					if(z < log10(zsh)){
						presyn[i]= presyn[i] + pow(10,sflx);
					}
					else{
						postsyn[i]= postsyn[i] + pow(10,sflx);
					}
				}
				//cflx_array[i] = cflx_array[i] + pow(10,cflx);
				fplot[i]= fplot[i] + pow(10,sflx) + pow(10,cflx);

				if (infosw == 1) {
					zonesFile.open("outputs/zones.dat", ios::app);
						if(!zonesFile.is_open()){
							cerr << "*** Error: can't open zones.dat" << endl;
						exit(1);
						} else{
							zonesFile << left << setw(20) << nutot[i] << setw(20) << sflx << setw(20) << cflx << endl;
						}
					zonesFile.close();
				}
			}//end ne for loop

               		gsl_spline_free(spline_snu), gsl_interp_accel_free(acc_snu);
                	gsl_spline_free(spline_cnu), gsl_interp_accel_free(acc_cnu);
		
		}//end z for loop
	}//end m[njet] for loop
	
	bbnrm	= (bbf1*rout)*(bbf1*rout) * cos(inclin)/(4.*dist*dist);
	blim	= log(rin);
	ulim	= log(rout);

        for(i=0; i<(ne-1); i++){
		frq     = pow(10,nutot[i]);
		if(bbsw	== 1 && disksw == 1){
			bbearthint(blim, ulim, frq, tin, rin, dist, inclin, bbflx);
			bbflx	= bbflx/mjy;
			bbtrm	= bbnrm*2.*herg*frq*frq*frq/(cee*cee*mjy*(exp(herg*frq/(kboltz*tbb2))-1.));
			fplot[i]= fplot[i] + bbflx + bbtrm;
			if(plotsw == 1){
				bbplot[i]= bbflx + bbtrm;
			}
		}

		if(plotsw == 1){
			if(fplot[i] == 0){
				fluxplotFile.open("outputs/total.dat", ios::app);
				if(!fluxplotFile.is_open()){
					cerr << "*** Error: can't open outputs/total.dat" << endl;
					exit(1);
				}
				else{
					fluxplotFile << left << setw(20) << nutot[i] << setw(20) << -200. << endl;
				}
				fluxplotFile.close();
			}
			else{
				fluxplotFile.open("outputs/total.dat", ios::app);
				if(!fluxplotFile.is_open()){
					cerr << "*** Error: can't open outputs/total.dat" << endl;
					exit(1);
				}
				else{
					fluxplotFile << left << setw(20) << nutot[i] << setw(20) << log10(fplot[i]) << endl;
				}
				fluxplotFile.close();
			}

			if(complot[i] == 0){
				complotFile.open("outputs/com.dat", ios::app);
				if(!complotFile.is_open()){
					cerr << "*** Error: can't open outputs/com.dat" << endl;
					exit(1);
				}
				else{
					complotFile << left << setw(20) << nutot[i] << setw(20) << -200. << endl;
				}
				complotFile.close();
			}
			else{
				complotFile.open("outputs/com.dat", ios::app);
				if(!complotFile.is_open()){
					cerr << "*** Error: can't open outputs/com.dat" << endl;
					exit(1);
				}
				else{
 					complotFile << left << setw(20) << nutot[i] << setw(20) << log10(complot[i]) << endl;
				}
				complotFile.close();
			}

			if(presyn[i] == 0){
				presynFile.open("outputs/presyn.dat", ios::app);
				if(!presynFile.is_open()){
					cerr << "*** Error: can't open outputs/presyn.dat" << endl;
					exit(1);
				}
				else{
					presynFile << left << setw(20) << nutot[i] << setw(20) << -200. << endl;
				}
				presynFile.close();
			}
			else{
				presynFile.open("outputs/presyn.dat", ios::app);
				if(!presynFile.is_open()){
					cerr << "*** Error: can't open outputs/presyn.dat" << endl;
					exit(1);
				}
				else{
					presynFile << left << setw(20) << nutot[i] << setw(20) << log10(presyn[i]) << endl;
				}
				presynFile.close();
			}

			if(postsyn[i] == 0){
				postsynFile.open("outputs/postsyn.dat", ios::app);
				if(!postsynFile.is_open()){
					cerr << "*** Error: can't open outputs/postsyn.dat" << endl;
					exit(1);
				}
				else{
					postsynFile << left << setw(20) << nutot[i] << setw(20) << -200. << endl;
				}
				postsynFile.close();
			}
			else{
				postsynFile.open("outputs/postsyn.dat", ios::app);
				if(!postsynFile.is_open()){
					cerr << "*** Error: can't open outputs/postsyn.dat" << endl;
					exit(1);
				}
				else{
					postsynFile << left << setw(20) << nutot[i] << setw(20) << log10(postsyn[i]) << endl;
				}
				postsynFile.close();
			}

			if(bbplot[i] == 0){
				bbplotFile.open("outputs/bb.dat", ios::app);
				if(!bbplotFile.is_open()){
					cerr << "*** Error: can't open outputs/bb.dat" << endl;
					exit(1);
				}
				else{
					bbplotFile << left << setw(20) << nutot[i] << setw(20) << -200. << endl;
				}
				bbplotFile.close();
			}
			else{
				bbplotFile.open("outputs/bb.dat", ios::app);
				if(!bbplotFile.is_open()){
					cerr << "*** Error: can't open outputs/bb.dat" << endl;
					exit(1);
				}
				else{
					bbplotFile << left << setw(20) << nutot[i] << setw(20) << log10(bbplot[i]) << endl;
				}
				bbplotFile.close();
			}
		}//end if plotsw

		if(fplot[i] == 0){
			photeng[i]= nutot[i];
			phot_spect[i]= -200.;
		}
		else{
			photeng[i]= nutot[i];
			phot_spect[i]= log10(fplot[i]);
		}
	}//end for loop ne

       
        /***************
         * Free Memory *
         ***************/
//        free(par2);
        delete[] ebin, delete[] nutot, delete[] nubb, delete[] nurad, delete[] energ, delete[] ephot, delete[] ephxr;
        delete[] rdlgen, delete[] rden, delete[] rdend, delete[] rdedn;
        delete[] dopfac, delete[] lelec, delete[] elen, delete[] ledens, delete[] etemp, delete[] dtemp, delete[] eled;
        delete[] phofrq, delete[] nusyn, delete[] nucom, delete[] comspc, delete[] comspc_src, delete[] nucom_src,  delete[] comi, delete[] comspc_it, delete[] phoint, delete[] synabs, delete[] shelen, delete[] shden;
        delete[] thmshock, delete[] thmcomp, delete[] thmbase,  delete[] plshock, delete[]plcomp, delete[] plbase; delete[] tdyntot, delete[] tsyntot, delete[] ebrarr; 
	delete[] tstrm, delete[] drtrm;
        delete[] phodis, delete[] nphot, delete[] jetu, delete[] disku;
        delete[] snu, delete[] sdump, delete[] cnu, delete[] cdump;
        delete[] complot, delete[] presyn, delete[] postsyn, delete[] fplot, delete[] bbplot;
        delete[] total;
	delete[] zed;


        gsl_spline_free(spline_syn), gsl_interp_accel_free(acc_syn);
        gsl_spline_free(spline_jet), gsl_interp_accel_free(acc_jet);
} //end xrbjet


/**
 * updated equipartition function
 * ***********************
 * @param mxsw          maxwellian switch [jet parameter]
 * @param equip         equipartition factor [jet parameter]
 * @param pspec         power-law index [jet parameter]
 * @param velsw		choice of velocity profile (necessary for BLJet equipartition)
 * @param gam0		initial bulk Lorentz factor (necessary for BLJet equipartition)
 * @return cnorm        normalization factor of the power-law distribution <!-- FIXME -->
 * @return ntot0        initial number of electron
 * @return b_en         magnetic energy in nozzle
 *
 */
void equipartition(double velsw, double mxsw, int nenpsw, double eta, double equip, double pspec, double nprot0, double emin, double emax, double endnsmj, double &cnorm, double &ntot0, double &b_en){
	cnorm   = 0;

    	if(nenpsw==0){
        	if(mxsw==1){
			if(equip==1){
				ntot0 = nprot0;
				b_en = ntot0*endnsmj;
			}
			else if(equip==0){
				ntot0 = nprot0;
				b_en = ntot0*(pmgm*cee*cee-endnsmj);
			} else{
				ntot0 = nprot0;
				b_en = equip*ntot0*endnsmj;
			}
		} else{
			if(equip==1){
				ntot0   = nprot0;
				cnorm = ntot0 * (1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
				if(pspec >= 1.99 && pspec <= 2.01){
					b_en = cnorm*log(emax/emin);
				} else{
					b_en = cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
				}
		} else if(equip==0){
                	ntot0 = nprot0;
                	cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
               		if(pspec >= 1.99 && pspec <= 2.01){
				b_en = ntot0*pmgm*cee*cee - cnorm*log(emax/emin);
                	} else{
				b_en = ntot0*pmgm*cee*cee - cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
                	}
		} else{
                	ntot0 = nprot0;
                	cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
                	if(pspec >= 1.99 && pspec <= 2.01){
				b_en = equip*cnorm*log(emax/emin);
                	} else{
				b_en = equip*cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
                	}
            	}
        }
	} else{
		if (velsw <= 1){
			if(mxsw==1){
				ntot0 = nprot0 * pmgm * cee * cee/((1.+equip) * endnsmj); // Now ntot0 calculated given Upr = Np*mp*c*c = Ue + Ub, where Ub/Ue = equip. Therefore ntot0 = Np*mp*c*c/(1+equip)*Ue
				b_en = equip * ntot0 * endnsmj;
			}else if(mxsw == 0){
				if(pspec >= 1.99 && pspec <= 2.01){
					cnorm = nprot0 * pmgm * cee * cee/((1.+equip) * log(emax/emin));
					b_en = equip * cnorm * log(emax/emin);
				} else{
					cnorm = nprot0 * pmgm * cee * cee * (2.-pspec)/((1.+equip) * (pow(emax,(2.-pspec))-pow(emin,(2.-pspec))));
					b_en = equip * cnorm  *  (pow(emax,(2.-pspec))-pow(emin,(2.-pspec)))/(2.-pspec);
			}
					ntot0 = cnorm * (pow(emax,(1.-pspec))-pow(emin,(1.-pspec)))/(1.-pspec);
			} else {
				ntot0 = nprot0 * pmgm * cee * cee/((1.+equip) * endnsmj);
           			b_en = equip * ntot0 * endnsmj;
           			cnorm = (ntot0 * endnsmj * (1.-mxsw) * (2.-pspec))/(pow(emax,(2.-pspec)) - pow(emin,(2.-pspec)));
			}
		} else { 
			ntot0 = nprot0*eta;
			if(mxsw==1){
				b_en = equip*ntot0*endnsmj;
			} else if (mxsw==0){
				cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
				if(pspec >= 1.99 && pspec <= 2.01){				
					b_en = equip*cnorm*log(emax/emin);
				} else {
					b_en = equip*cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
				}
			} else {
				cnorm = (1.-mxsw)*ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
				if(pspec >= 1.99 && pspec <= 2.01){				
					b_en = equip*(mxsw*ntot0*endnsmj+cnorm*log(emax/emin));
				} else {
					b_en = equip*(mxsw*ntot0*endnsmj+cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec));
				}
			}
		}
//}
	}
}

/**
 * Main function of the jet properties
 * 
 * **************************************
 *
 * @param z		distance from origin along jet axis (same units as r0, h0
 * @param b0		magnetic field in nozzle
 * @param n0		density in nozzle
 * @param g0		electron Lorentz factor in nozzle
 * @param r0		width of nozzle
 * @param h0		height of nozzle
 *
 */

void jetpars(double velsw, double zeta, double z0, double z, double zsh, double r_g, double b0, double n0, double g0, double r0, double h0, gsl_spline *spline, gsl_interp_accel *acc, double &gb, double &dgb, double &b, double &n, double &g, double &r){
    
	double y	= 0;
    	double dy	= 0;
    	double mj, gbs0, Gammas, betas0;    
	double x, zff, gbin, bff;
        
        Gammas	= 4./3.;
        betas0	= sqrt(Gammas*(Gammas-1.)/(Gammas+1.));
        gbs0	= betas0/sqrt((1.-betas0)*(1.+betas0));
        x	= log10((max(z-h0,0.)+r0)/r0);
        zff	= pow(10, z0);
        
        if(z < h0){
        	y	= gbs0;
        }
        else{
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
        }
        else{
            	if (velsw== 0) {
                	b = b0*(r0/r)/pow(mj/zeta, 0.5+1./6.);
                	g = g0/pow(mj, 1./3.)/pow(r/r0, 2./3.);
            	} else {
                	b = b0*sqrt(zeta)*(r0/r)/pow(mj, 0.5+1./6.);
                	g = g0/pow(mj, 1./3.);
            	}
        } 
}

void bljetpars(double mxsw, double velsw, double z, double zacc, double r_g, double eta, double n0, double g0, double r0, double h0, double endnsmj, double pspec, double cnorm, double emin, double emax, double ebreak, gsl_spline *spline, gsl_interp_accel *acc, double &gb, double &b, double &n, double &g, double &r, double sigsh){
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
	n = eta*n0*pow(sin(0.15/gam0)/sin(theta),2.)*pow(min(h0/z,1.),2.)/mj;

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

void b_prof(double mxsw, double gfin, double eta, double n, double endnsmj, double pspec, double cnorm, double emin, double emax, double ebreak, double &field, double g, double sigsh){
	double betas0, g0, sig0, sig, w;
	double Gammas = 4./3.;
	betas0 = sqrt(Gammas*(Gammas-1.)/(Gammas+1.));
   	g0 = 1./sqrt((1.-betas0)*(1.+betas0));

	/*
	NOTE: for simplicity the magnetic field is always computed assuming a purely thermal electron distribution. 
        Because the bulk of the enthalpy is in protons and magnetic field, this is acceptable approximation.
        This also means that the pair content of the jet has to be low (10 or so) in order to be fully self consistent.
	IF ZSH < ZACC MAKE SURE THAT THE AVERAGE ELECTRON LORENTZ FACTOR IS LOW
	*/

	w = Gammas*eta*n*endnsmj;
	sig0 = (1.+sigsh)*gfin/g0 - 1.;
	sig = (g0/g)*(1.+sig0)-1.;
	field = sqrt(sig*4.*pi*(n*pmgm*cee*cee+w));
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
	
	/**
	 * This is what the jet sees from BB in stationary frame.
	 * It's the effective temperature as input.
	 *
	 */
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

	/**
	 * Since we are interested only in the photon energy density magnitude and not direction,
	 * F dot n is equivalent to absolute value
	 *
	 */
	return abs(bbfunc);
}


/**
 * bb_func3 is for inside the disk
 *
 */
double bbfnc3(double thet, void *p){
        struct bbfnc3_params *params = (struct bbfnc3_params *)p;
        double gamv     = (params -> gamv);
        double tin      = (params -> tin);
        double nu       = (params -> nu);

	double beta, tineff, fac, bbfunc3;
	
	/**
	 * This is what the jet sees from BB in stationary frame. see pg 256, iv and R&L pg 153.
	 * It's the effective temperature as input.
	 *
	 */
	beta	= sqrt(gamv*gamv-1.)/gamv;
	tineff	= tin/(gamv*(1.+beta*cos(thet)));
	fac	= herg*nu/(kboltz*tineff);

	if(fac < 1.e-3){
		bbfunc3	= sin(thet)*cos(thet)/fac;
	}
	else{
		bbfunc3	= sin(thet)*cos(thet)/(exp(fac)-1.);
	}

	/**
	 * Since we are interested only in the photon energy density magnitude and not direction,
	 * F dot n is equivalent to absolute value
	 *
	 */
	bbfunc3	= abs(bbfunc3);

	return bbfunc3;
}

/**
 * bbdisk gives flux of the top of the disk as a function of the frequency and
 * assumes that hbb is constant, set by inner radius
 *
 */
double bbdisk(double lnu, double gamv, double tin, double rin, double rout, double z, double hbb){
	double nu, blim, ulim, inside, top;

	
	/*if (lnu > 1.e4){											// For some reason somewhere in the code the frequency is not 
		lnu = log10(lnu); 										// converted to log scale, this makes sure it always happens.
		//cout << "Converted frequencies to log scale, something wierd happened with the disk!" << endl;// The spectrum is not affected, but debugging tools you may run
	}	*/												// get annoyed by it. This is a horrible hack and should be fixed.
	nu	= exp(lnu);	
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
 * <!-- FIXME add comments -->
 *
 */
double bbdiskfnc(double lnu, void *p){
        struct bbdiskfnc_params *params = (struct bbdiskfnc_params *)p;
        double rin      = (params -> rin);
        double rout     = (params -> rout);
        double z        = (params -> z);
        double hbb      = (params -> hbb);
        double gamv     = (params -> gamv);
        double tin      = (params -> tin);

        return bbdisk(lnu, gamv, tin, rin, rout, z, hbb);
}

/**
 * bbintegral calculates the integration over log(nubot)-log(nutop) of bbdisk
 * <!-- FIXME add comments -->
 *
 */
void bbintegrals(double lnub, double lnut, double gamv, double tin, double rin, double rout, double z, double hbb, double &bbint){
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
 * This function is for synchrotron emissivity to be integrated over
 * particle distribution in synint
 *
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
 * This function is for self-absorption to be intergrated over the
 * derived function of the particle distribution in synint
 *
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
 * Synintegrals does frequency and electron spectrum dependent calc, single component
 * returns necessary for emissivity and self-absorption for use in synint,
 * which calculates the actual spectrum
 *
 * @param freq	Frequency of emission
 *
 * @return esum	Synchrotron emission kernal at frequency 'freq'
 * @return asum Synchrotron absorption kernal at frequency 'freq'
 *
 */
void synintegrals(double freq, double bfield, double elenmn, double elenmx, gsl_spline *spline_syn, gsl_interp_accel *acc_syn, gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_derivs, gsl_interp_accel *acc_derivs, double &esum, double &asum){ //<!-- FIXME -->
        //esum
        gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(1000);
        double result1, error1;
        gsl_function F1;
        struct synemis_params F1params = {freq, bfield, spline_syn, acc_syn, spline_eldis1, acc_eldis1};
        F1.function     = &synemis;
        F1.params       = &F1params;
	
        gsl_integration_qag(&F1, log(elenmn), log(elenmx), 1e-1, 1e-1, 1000, 1, w1, &result1, &error1);
	esum    = result1;
        gsl_integration_workspace_free(w1);      //<!-- FIXME should I free the workspace here? -->

        //asum
        gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(1000);
        double result2, error2;
        gsl_function F2;
        struct absfnc_params F2params = {freq, bfield, spline_syn, acc_syn, spline_derivs, acc_derivs};
        F2.function     = &absfnc;
        F2.params       = &F2params;
        gsl_integration_qag(&F2, log(elenmn), log(elenmx), 1e-1, 1e-1, 1000, 1, w2, &result2, &error2);
        asum    = result2;
        gsl_integration_workspace_free(w2);
}


/**
 * For selfabsorption, single component
 * returns specific intensity in erg/s/cm^2/st/Hz
 * integrates distributions over electron spectrum
 * sin(pitch) = 2/3 equiv to angle averaging over isotropic
 *
 * @param freq		Frequency of emission
 * @param r		d
 * @param angle		d
 * @param dopfac	d
 * @param h		Height in the jet along z-axis
 * @param bfield        magnetic field magnitude
 * @param dist          distance to the source <!-- FIXME -->
 * @param esum		Synchrotron emission kernal
 * @param asum		Synchrotron absportion kernal
 *
 * @return fluxa	Synchrotron self-absorption emission
 * @return fluxa2	Resulting photon density
 *
 */
void synint(double freq, double r, double angle, double dopfac, double h, double bfield, double dist, double &esum, double &asum, double &fluxa, double &fluxa2){
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

	/*
	 * This first term is for what is seen externally. It includes skin depth
	 * and angle effects. See pg 12 ntbk 4 for consideration.
	 */
	tsyn	= pi/2. * asyn*r/(dopfac*sin(angle));
	if(tsyn >= 1.){
		absfac	= (1.-exp(-tsyn));
	}
	else{
		absfac	= tsyn*(1.+tsyn*(-0.5+tsyn*(1./6.+tsyn*(-1./24.+tsyn/120))));
	}

	fluxa	= (2.*r*h*sin(angle)*dopfac*absfac*epsasyn)/area;

	/*
	 * This second term is the same as above, but assuming what is "seen" locally
	 * by the particles for compton scattering. Pass 'flux', erg/s/Hz to main code
	 */
	tsyn2	= pi/2.*asyn*r;
	if(tsyn2 >= 1.){
		absfac2	= (1.-exp(-tsyn2));
	}
	else{
		absfac2	= tsyn2*(1.+tsyn2*(-0.5+tsyn2*(1./6.+tsyn2*(-1./24.+tsyn2/120))));
	}

	fluxa2	= pi*r*r*absfac2*epsasyn;
}


/**
 * This function is the kernel of eq 2.48 in Blumenthal & Gould(1970)
 *
 */
double comfnc(double ein, void *p){
        struct comfnc_params *params = (struct comfnc_params *)p;
        double game     = (params -> game);
        double e1       = (params -> e1);
        gsl_spline *spline_com1 = (params -> spline_com1);
        gsl_interp_accel *acc_com1 = (params -> acc_com1);

	double com      = 0;
	double einit, utst, biggam, q, elg;
	double phonum	= 0;
	double tm1, tm2, tm3;
	double btst, eg4;

	einit	= exp(ein);
	btst	= einit/(game*emerg);
	eg4	= 4.*einit*game;
	utst	= eg4/(emerg+eg4);

	/**
	 * [JW]
	 * the following costs 3 seconds in a typical run
	 * are these sanity checks really needed?
	 */
	if((e1 < btst && (btst-e1)/btst >= 4.e-8) || (e1 > utst && (e1-utst)/utst >= 4.e-8)){
		com	= 0.;
		return com;
	}
//    cout  << einit << "\t";
	biggam	= eg4/emerg;
	q	= e1/(biggam*(1.-e1));
	elg	= ein*0.434294481903252; // log10(exp(1)) = 0.434294481903252

//    cout  << elg << "\t";
    
    phonum  = gsl_spline_eval(spline_com1, elg, acc_com1);
//    cout << "Phonum: \t";
//   cout << elg << ", " << phonum << "\t\t";

	tm1	= 2.*q*log(q);
	tm2	= (1.+2.*q)*(1.-q);
	tm3	= 0.5*((biggam*q)*(biggam*q))*(1.-q)/(1.+biggam*q);

	com	= (tm1+tm2+tm3)*pow(10.,phonum);
    
//    cout << com << " " << pow(10.,phonum) << endl;
    
	return com;
}

/**
 * This function is for SSC of synchro flux -- numerically integrates over
 * photon distribution and cross-section
 *
 */
double comint(double gam, void *p){
        struct comint_params *params = (struct comint_params *)p;
        double eph      = (params -> eph);
        double ephmin   = (params -> ephmin);
        double ephmax   = (params -> ephmax);
        gsl_spline *spline_eldis1 = (params -> spline_eldis1);
        gsl_interp_accel *acc_eldis1 = (params -> acc_eldis1);
        gsl_spline *spline_com1 = (params -> spline_com1);
        gsl_interp_accel *acc_com1 = (params -> acc_com1);

	double game, econst, gnorm, blim, ulim, value, ellog, e1;
	double den	= 0; 

	game	= exp(gam);
	e1	= eph/(game*emerg);

	econst	= 2.*pi*re0*re0*cee;
	gnorm	= emerg;

	// DON'T CHANGE!  OTHERWISE GOES OUT OF ARRAY BOUNDS
	// added incremental to avoid hitting blim/utst boundary
	blim	= max(log(eph/(4.*game*(game-eph/emerg))), log(ephmin))+1.e-14;
	ulim	= log(min(eph, ephmax));

//    cout << "integral limits in comint: " << blim << " " << ulim << endl;
    
	if(ulim <= blim){
		return 0;
	}

	// DON'T CHANGE THIS EPS!!!
	//value	= qromb2(comfnc, 5.e-6, blim, ulim);
	gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(1000);
   	double result2, error2;
    	gsl_function F2;
   	struct comfnc_params F2params = {game, e1, spline_com1, acc_com1};
    	F2.function     = &comfnc;
    	F2.params       = &F2params;
    	gsl_integration_qag(&F2, blim, ulim, 1e-1, 1e-1, 1000, 1, w2, &result2, &error2);
    	value   = result2;
    	gsl_integration_workspace_free(w2);

    	ellog	= log10(game*emerg);
    	den     = gsl_spline_eval(spline_eldis1, ellog, acc_eldis1);
    	den     = pow(10., den);
//   	 cout << "Interpolation error here? \n";
//   	 cout << den << endl;

      	  //gsl_spline_free(spline_eldis1), gsl_interp_accel_free(acc_eldis1);
      	  //gsl_spline_free(spline_com1), gsl_interp_accel_free(acc_com1);

	return gnorm*econst*den*value/game;
}

double comintegral(double blim, double ulim, double eph, double ephmin, double ephmax, gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_com1, gsl_interp_accel *acc_com1){
        gsl_integration_workspace *w    = gsl_integration_workspace_alloc(1000);
        double result, error;
        gsl_function F1;
        struct comint_params F1params   = {eph, ephmin, ephmax, spline_eldis1, acc_eldis1, spline_com1, acc_com1};
        F1.function     = &comint;
        F1.params       = &F1params;
        gsl_integration_qag(&F1, blim, ulim, 1e-1, 1e-1, 1000, 1, w, &result, &error);
        gsl_integration_workspace_free(w);

        return result;
}


/**
 * bbearth gives the flux at Earth as a function of r(T) and freq at distance dist
 *
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
 *
 */
void bbearthint(double blim, double ulim, double frq, double tin, double rin, double dist, double inclin, double & bbflx){
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

/**
 * Modified Bessel function
 * <!-- FIXME get rid of these
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
		i0	= 1. + t2*(3.5156229+t2*(3.0899424+t2*(1.2067492*pow(t,6)+t2*(0.2659732*pow(t,8)+t2*(0.0360768+t2*0.0045813)))));
		k0	= log(2./x)*i0 - 0.5772156649 + y2*(0.42278420+y2*(0.23069756+y2*(0.03488590+y2*(0.00262698+y2*(0.00010750+y2*0.00000740)))));
	}
	else{
		k0	= exp(-x)/sqrt(x) * (1.25331414+y1*(-0.07832358+y1*(0.02189568+y1*(-0.01062446+y1*(0.00587872+y1*(-0.00251540+y1*0.00053208))))));
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
		i1=x*(0.5e0+t2*(0.87890594e0 + t2*(0.51498869e0 + t2*(0.15084934e0 + t2*(0.02658733e0 + t2*(0.00301532e0 + t2*0.00032411e0))))));
		 k1 = 1.e0/x*(x*log(y)*i1+1.e0 + y2*(0.15443144e0-y2*(0.67278579e0 - y2*(0.18156897e0-y2*(0.01919402e0 - y2*(0.00110404e0-y2*0.00004686e0))))));
	}
	else{
		k1 = exp(-x)/sqrt(x)*(1.25331414e0+y1*(0.23498618e0 - y1*(0.03655620e0+y1*(0.01504268e0 - y1*(0.00780353e0+y1*(0.00325614e0 - y1*0.00068245e0))))));
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
void xrbinterp(double *ear, double *energ, double *phot, double *photar, int ne, int newne){
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
