// =====================================================================================
// 
//       Filename:  agnjet.cc
// 
//    Description:  Main file of the leptonic model of a black hole jet.
//                  Work based on the agnjet.f code from S.B. Markoff.
//					The description fo each function is in agnjet.hh
// 
//        Version:  2.0
//        Created:  09/15/2012 13:36:59
//       Revision:  04/09/2018
//       Compiler:  gcc
// 
//         Author:  Samia Drappeau (sd), drappeau.samia@gmail.com
//	 Rewritten by:  Chiara Ceccobello
//  Maintained by:  Matteo Lucchini, m.lucchini@uva.nl
//        Company:  API - UvA
//
// =====================================================================================

#include "agnjet.hh"

using namespace std;

/**
 * bhjet (ear, ne, params, phot_spect)
 * Main function of the code
 * ************************************
 * @param ear           energy array
 * @param ne            size of ear
 * @param params        model's parameters
 *
 * @return photend      photon's energy array
 * @return phot_spect   photon's spectra
 */
 
void bhjet(double* ear,int ne,double *param,double *photeng,double *phot_spect){

    int i, k, m; 										//for-loop variables
    double sum1, sum2, tmp1, tmp2;						//temporary stuff for numerics

    bool isShock         = false;						//flag for shock
    bool isBreaknjetnsyn = false;						//flag for synchrotron for loop
    bool isVerbose       = false;						//flag for verbose debug mode
    bool isSSC           = false;						//when false: multiple SSC
    bool isPair          = false;						//when false no pair rate is calculated

    //Synchrotron - Compton - Gamma-ray grid [Hz]
    int nsyn = 100;										//number of bins in Synchrotron grid
    int ncom = 60;										//number of bins in Compton grid
    double snumin;										//lower boundary Synchrotron grid
    double snumax;										//upper boundary Synchrotron grid
    double cnumin;										//lower boundary Compton grid
    double cnumax;										//upper boundary Compton grid
    double snuinc;										//increment of Synchrotron grid
    double cnuinc;										//increment of Compton grid

    //Fixed quantities for computations
    int nz 			= 70;								//number of zones in jet
    int zone_zcut 	= 0;								//index of zone where grid switches to logarithmic
    double thrlum; 										//disk thermal luminosity
    double eddrat; 										//accretion rate in Eddington units 
    double thmfrac; 									//1-plfrac (see jet's parameter plfrac down)
    double zcut; 										//distance above which Compton is neglected
    double theff; 										//angle to the secondary blackbody
    double tbbeff;										//effective blackbody temperature
    int zfrac       = 100;								//multiplier to set Comptonization region
    int bbsw    	= 1;								//switch to get black-body radiation
    int disksw  	= 1;								//switch to get disk radiation
    double zmin;										//minimum distance from BH
    int nzdum   	= 100;								//size of jet's segment array
    double bbf1    	= 0;								//black-body fraction (used for second BB)
    double reff    	= 0;								//bff1*rout (used for second BB)
    double reff2   	= 0;								//reff**2		
    int njet    	= 2;								//number of jets produced
    int Niter		= 100;								//number of IC scatterings
    int sizegb 		= 55;								//size of velocity/distance array before interpolation
    double zeta		= 1;								//Uj = n mc2, see Crumley et al. 2017
    double sig0;										//Initial magnetization for magnetic acceleration

	//Variables related to the speed of sound in nozzle
    double gad4_3; 										//relativistic ideal gas adiabatic index
    double betas0; 										//speed of sound in the nozzle
    double gam0;										//initial jet Lorentz factor 1./sqrt(1.-(betas0^2))
    double rvel0;										//gam0*betas0
    double vel0;										//cee*rvel0
    double gamax0; 										//electron Lorentz factor in the nozzle

    //Variables related to particles and energies
    double eta;											//pair content in jet if using magnetic acceleration
    double emin; 										//minimum electron energy
    double emax; 										//maximum electron energy
    double ebreak; 										//break electron energy
    double gamin; 										//minimum electron Lorentz factor
    double uplim; 										//upper limit in thermal energy density calculation
    double endncom; 									//average thermal Lorentz factor
    double endnsmj;										//electron energy at peak of M-J distribution
    double nprot0;										//initial number density of protons 
    double ntot0;										//initial number density of electrons
    double b_en;										//magnetic energy density in nozzle
    double b0;											//magnetic field strength in nozzle
    double cnorm;										//power-law distribution normalization
    double endens;										//initial electron energy density
    double betat;										//emerg/(kboltz*eltemp)
    double enorm;										//ntot0*betat/(emerg*k2_fnc(betat));
    double einc;										//energy incrementation in thermal distribution
    double bete;										//electron beta in distribution
    int nelec		= 200;								//electron array size in each jet segment
    double elenmn;										//minimum electron energy for radiation integration
    double elenmx;										//maximum electron energy for radiation integration
    int nenpsw 		= 1;								//switch for change between normalizations when using
    													//velsw = 0 or 1, leave to 0 if you dont have a good
    													//reason not to do so

    //Variables related to each jet segment
    double z; 											//distance along the jet
    double delz;										//height of segment
    double area;										//jet transverse area 
    double vol;											//volume of segment
    double rvel;										//jet velocity along z			
    double gamax;										//shift parameterin peak MJ distribution
    double r;											//jet radius
    double bfield;										//magnetic field strength
    double ntot;										//electron number density
    double gamv2;										//bulk Lorentz factor squared
    double gamv;										//bulk Lorentz factor
    double beta;										//jet velocity in units of c
    double gshift;										//actual shift, gamax/gamax0
    double gshock;										//peak MJ value at shock
    double ub;											//magnetic energy density 
    double ucom;										//rough estimate of Urad for IC
    int nw;												//counter to track lost particles when cooling
    double uphdil; 										//some measure of photons from disk, CHECK UNITS
    double uphdil2;										//similar to above, CHECK UNITS
    double gemax;										//maximum electron Lorentz factor
    double bemax;										//maximum electron velocity in units of c
    double gebreak;										//break electron Lorentz factor
    double mjteff;										//peak of MJ distribution when cooling
    double pltrm;										//power-law distribution term for a given energy
    double shocknd;										//number density at shock region
    double shockr;										//shock region radius
    double shockgb;										//shock region cooling shift
    double avelec;										//average electron energy, only used for debugging
    double avphot;										//average photon energy, only used for debugging
    double ephmin;										//lower integration boundary for first Compton
    double ephmax;										//higher integration boundary for first Compton
    double hbb;											//disk height as a function of accretion rate

    //Variables related to pair estimates    
    double rate_pa; 									//pair annihilation rate
    double rate_pp;										//pair production rate
    double n_pa;										//number density of pairs annihilated
    double n_pp;										//number density of pairs produced
    double ratio_ne;									//ratio of number density to initial number density

	//Definition of all pointers used in code
	double *gbx 	= new double[sizegb-1]();			//distance array for velocity profile
	double *gby 	= new double[sizegb-1]();			//gamma*beta array for velocity profile
	double *zed 	= new double[nz](); 				//distance array, needed in the radiation calculation
    double *nutot	= new double[NEBIN]();				//frequency array for total spectrum
    double *nubb	= new double[nsyn]();				//frequency array for black body spectrum
    double *nurad	= new double[nsyn]();				//pow(10.,nubb)
    double *energ	= new double[NEBIN]();				//herg*nurad
    double *ephot	= new double[nsyn]();				//log10(energ)
    double *ephxr	= new double[ncom]();				//X-ray photon energy array for Compton 
    double *dopfac	= new double[njet*nzdum]();			//Doppler factor in each jet segment
    double *rdend	= new double[nelec]();				//electron energy density 
    double *rdedn	= new double[nelec]();				//electron energy distribution, log10(reled)
    double *rdlgen	= new double[nelec]();				//log10(electron energy array)  
    double *rden	= new double[nelec]();				//electron energy array
    double *lelec	= new double[nelec]();				//electron energy array
    double *elen	= new double[nelec]();				//electron energy used in each step of z loop
    double *ledens	= new double[nelec]();				//log10(electron number density)
    double *etemp	= new double[nelec]();				//temporary electron energy array
    double *dtemp	= new double[nelec]();				//temporary electron energy density array
    double *eled	= new double[nelec]();				//temporary energy*density array
    double *nusyn	= new double[njet*nzdum*nsyn]();	//frequency of synchrotron photons
    double *nucom	= new double[njet*nzdum*ncom]();	//frequency of Compton photons
    double *comspc	= new double[njet*nzdum*ncom]();	//Compton spectrum pre-final dump
    double *totcom = new double[ncom]();				//total Compton flux for pair estimate
    double *synabs	= new double[njet*nzdum*nsyn]();	//self-absorbed synchrotron emission
    double *shelen	= new double[nelec]();				//elen at shock
    double *shden	= new double[nelec]();				//total electron density array at shock
    double *thmshock = new double[nelec]();				//thermal component array at shock
    double *thmcomp = new double[nelec]();				//thermal component array in each segment
    double *thmbase = new double[nelec]();				//thermal component array at base
    double *plshock = new double[nelec]();				//powerlaw component array at shock
    double *plcomp = new double[nelec]();				//powerlaw component array at in each segmetn
    double *plbase = new double[nelec]();				//powerlaw component array at base
    double *tstrm	= new double[nelec]();				//first term in derivative of electron array
    double *drtrm	= new double[nelec]();				//second term in derivative of electron array
    double *phodis	= new double[nsyn]();				//photon distribution to use for Compton 
    double *nphot	= new double[nsyn]();				//synchrotron seed photons for Compton
    double *jetu	= new double[nsyn]();				//comoving radiation energy density 
    double *disku	= new double[nsyn]();				//disk radiation energy density
    double *complot	= new double[NEBIN]();				//flux array for Compton plot
    double *cflx_array = new double[NEBIN]();			//temporary flux array for Compton plot
    double *presyn	= new double[NEBIN]();				//flux array for preshock synchrotron plot
    double *postsyn	= new double[NEBIN]();				//flux array for postshock synchrotron plot
    double *bbplot	= new double[NEBIN]();				//flux array for blackbody plot
    double *total	= new double[NEBIN]();				//flux array for total plot
    double *fplot	= new double[NEBIN]();				//frequency array for plot   

    //Free parameters
    double mbh;											//black hole mass in solar masses
    double r_g; 										//black hole gravitational radius in cm
    double eddlum; 										//black hole Eddington luminosity
    double inclin;   									//jet viewing angle wrt line of sight
    double dist;										//luminosity distance in kpc
    double jetrat;										//injected jet power in Eddington units
    double r0;											//nozzle radius in r_g
    double hratio;										//nozzle aspect ratio
    double h0;											//nozzle height
    double zsh;											//shock region location in r_g
    double zacc;										//acceleration distance if using magnetic acceleration
    double zmax;										//maximum jet length in r_g
    double eltemp;										//thermal electron peak Lorentz factor in the nozzle   
    double pspec;										//power-law distribubtion slope
	double heat; 										//shock heating; at shock eltem = heat*eltemp
	double brk; 										//break energy of non-thermal particles
    double fsc;											//acceleration efficiency parameter
    double gamfac;										//maximum injected energy of PL particles
    double beta_pl;										//plasma beta at base Ue/Ub
    double equip;										//1/beta_pl
	double sigsh;										//final magnetization at the shock, if velsw > 1
	double tin;											//inner accretiond disk temperature(if>0)/-Mdot(if<0)
	double rin; 										//inner accretion disk radius
    double rout;										//outer accretion disk
    double tbb2;										//secondary blackbody temperature
	double bbf2;										//secondary blackbody normalization	
    double plfrac;										//percentage of particles accelerated into powerlaw
	double mxsw;										//initial particle distribution switch
    double velsw;										//Velocity profile used, see README
	int plotsw;											//plot switch
    int infosw;											//print info to screen switch  

    //New tabulated velocity for 1D quasi-isothermal Bernoulli eq.
    static double gbx_vel2[54] = {1.0, 1.00001, 1.00005, 1.00023, 1.00101, 1.00456, 1.02053, 1.09237,1.41567,
    2.87053, 6.26251, 14.3691, 34.0825, 82.8831, 205.572, 518.283, 1325.11, 3429.81, 8975.1, 23719.5,
    63255.,170099., 460962., 1.25824e+06, 3.4578e+06, 9.5633e+06, 2.66095e+07, 7.44655e+07,2.09528e+08,
    5.92636e+08, 1.6846e+09, 4.81148e+09, 1.38055e+10, 3.9787e+10, 1.15153e+11, 3.34651e+11, 9.76408e+11,
    2.85981e+12, 8.40728e+12, 2.4805e+13, 7.34419e+13, 2.18185e+14, 6.50341e+14, 1.94471e+15, 5.83352e+15,
    1.75523e+16, 5.29701e+16, 1.60321e+17, 4.86616e+17, 1.48111e+18, 4.52032e+18, 1.38326e+19, 4.24394e+19,
    1.0e+20};

    static double gby_vel2[54] = {0.485071, 1.05031, 1.05032, 1.05039, 1.05067, 1.05193, 1.05751, 1.08105,
    1.16389, 1.35278, 1.52406, 1.68077, 1.82495, 1.95888, 2.08429, 2.20255, 2.3147, 2.42158, 2.52386,
    2.62205, 2.71662, 2.80793, 2.89628, 2.98195, 3.06516, 3.14611, 3.22497, 3.30189, 3.37702, 3.45046,
    3.52232, 3.59270, 3.66169, 3.72937, 3.79580, 3.86106, 3.92520, 3.98827, 4.05033, 4.11143, 4.17160,
    4.23089, 4.28933, 4.34696, 4.40381, 4.45992, 4.51530, 4.56999, 4.62402, 4.67740, 4.73015, 4.78230,
    4.83388, 4.87281};
    
    //Synchrotron tables for interpolations     
    static double arg[47]	= {0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.03, 0.05, 0.07, 0.1,
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5,
    8., 8.5, 9., 9.5, 10., 12., 14., 16., 18., 20., 25., 30., 40., 50.};
    
    static double var[47]	= {-1.002e+0, -9.031e-1, -7.696e-1, -6.716e-1, -5.686e-1, -4.461e-1, -3.516e-1,
    -2.125e-1, -1.537e-1, -1.192e-1, -8.725e-2, -4.383e-2, -3.716e-2, -4.528e-2, -5.948e-2, -7.988e-2,
    -1.035e-1, -1.296e-1, -1.586e-1, -1.838e-1, -3.507e-1, -5.214e-1, -6.990e-1, -8.861e-1, -1.073e+0,
    -1.267e+0, -1.470e+0, -1.670e+0, -1.870e+0, -2.073e+0, -2.279e+0, -2.483e+0, -2.686e+0, -2.893e+0,
    -3.097e+0, -3.303e+0, -3.510e+0, -3.717e+0, -4.550e+0, -5.388e+0, -6.230e+0, -7.075e+0, -7.921e+0,
    -1.005e+1, -1.218e+1, -1.646e+1, -2.076e+1};

    //Initialization of the synchrotron interpolation
    gsl_interp_accel *acc_syn       = gsl_interp_accel_alloc();
    gsl_spline *spline_syn          = gsl_spline_alloc(gsl_interp_cspline, 47);
    gsl_spline_init(spline_syn, arg, var, 47);   

	//Assign free parameter values
	mbh		= param[0];
    r_g		= gconst*mbh*msun/(cee*cee);
    eddlum	= 1.25e38 * mbh;
	inclin	= param[1]*pi/180;
	dist	= param[2]*kpc;        
    jetrat	= eddlum*param[3];
	r0		= param[4]*r_g;
	hratio	= param[5];
	zsh		= param[6]*r_g;
	zacc	= param[7]*r_g;
	zmax	= param[8]*r_g;
	eltemp	= param[9]*emerg/kboltz; 
	pspec	= param[10];
	heat	= param[11];
    brk		= param[12];
    fsc		= param[13]; 
    gamfac	= param[14];  
    beta_pl	= param[15];	
	sigsh	= param[16];
	tin		= param[17]/kboltz_kev; 
	rin		= param[18]*r_g;
	rout	= param[19]*r_g;
	tbb2	= param[20];
	bbf2	= param[21];
    plfrac	= param[22];
    mxsw	= param[23];
	velsw	= param[24];
	plotsw	= param[25];	
	infosw	= param[26];

	//Creating energy arrays for synch, compton, total calculations 
    ene_arrays(isSSC,ne,nsyn,ncom,zfrac,mbh,inclin,njet,snumin,snumax,snuinc,cnumin,cnumax,cnuinc,ear,nutot,
    nubb,nurad,energ,ephot,ephxr);
	
	//Initializing disk quantities
	disk_init(infosw,jetrat,rin,rout,eddlum,bbf2,tbb2,tin,disksw,thrlum,eddrat,bbf1,reff,reff2,hbb);		
	
	//Initializing jet quantities	
	jet_init(zfrac,sizegb,mxsw,velsw,jetrat,r_g,r0,hratio,zacc,zmax,beta_pl,equip,h0,zsh,zcut,gad4_3,betas0,
	gam0,rvel0,vel0,zmin,gbx,gby,gbx_vel2,gby_vel2);
	

   	//Setting up velocity profile spline
 	//WARNING: if you are using Velprof rather than the table you need to change spline into akima
    gsl_interp_accel *acc_jet       = gsl_interp_accel_alloc();
    gsl_spline *spline_jet          = gsl_spline_alloc(gsl_interp_akima, sizegb-1);

	gsl_spline_init(spline_jet, gbx, gby, sizegb-1);
	
    //Initiate plot files
    create_plotFiles(plotsw,infosw);

    /*********************************
     * START RUNNING BHjet model *
     *********************************/

	//Initializing particle distribution and proton number density
	particles_init(mxsw,velsw,jetrat,r0,eltemp,pspec,gamfac,sigsh,plfrac,gad4_3,vel0,gam0,thmfrac,emin,gamin,
	ebreak,emax,uplim,endncom,endnsmj,betat,nprot0,eta,equip,sig0,beta_pl);

    //Calculating magnetic energy density from equipartition assumptions
	equipartition(velsw,mxsw,nenpsw,eta,equip,pspec,nprot0,emin,emax,endnsmj,cnorm,ntot0,b_en,b0);

    //Calculating the electron distribution
    electrons(mxsw,nelec,pspec,gamfac,endnsmj,ntot0,cnorm,eltemp,endens,gamax0,rdlgen,rdedn,thmbase,plbase,
    bete);

    //If running with nenpsw=0, print out the jetrat (jet power) you would have if running with nenpsw=1. 
    //This can be used to calculate what new jetrat value you need when using an old parameter file in jetwrap
    //move this in equipartition
    if(nenpsw == 0) {
		printf("old jetrat=%e\n",jetrat/eddlum);
        jetrat = (endens + b_en) * (4.*vel0*pi*r0*r0);
        printf("Renormalized jetrat using old values=%e\n", jetrat/eddlum);
        
        //Check new normalizations
        printf("new ntot0, nprot0 =%e,%e\n",ntot0,nprot0);
        printf("Ue,Ub=%e,%e\n", endens,b_en);
        printf("Check 2: Up,Ue+Ub,Ub/Ue=%e,%e,%e\n",nprot0*pmgm*cee*cee,endens+b_en,b_en/endens);
        printf("new bfield=%e\n",sqrt(b_en)*8.*pi);
        printf("endens,new U_b=%e,%e\n",endens,equip*endens);        
    }
    //Quality of life output of physical parameters if infosw=1
	if (infosw ==1) {
		cout << "Physical quantities: " << endl;
		cout << "Maxwell-Juttner peak in the nozzle: " << 2.23*eltemp << endl;
		cout << "Inner disk temperature in K: " << tin << "; accretion rate:" << eddrat << endl;
		cout << "Pair content: " << eta << endl;
	}  

    //--------------------------------------------------------------------------------------------

    /**
     * START THE BIG LOOP OVER Z
     */

    for(k=0; k<nz; k++){

		//Adjusted grid to avoid resolution issues
		jetgrid(k,infosw,nz,r_g,velsw,h0,r,zmin,zcut,zmax,zone_zcut,delz,z,zed);

        //Define properties of jet's segment at z, given initial conditions
        if (velsw <=1){
        	jetpars(velsw,zeta,zmin,z,zsh,r_g,b0,ntot0,gamax0,r0,h0,spline_jet,acc_jet,rvel,bfield,ntot,gamax,
        	r);
        } else {
        	bljetpars(mxsw,velsw,z,zacc,r_g,eta,nprot0,gamax0,r0,h0,endnsmj,pspec,cnorm,emin,emax,ebreak,
        	spline_jet,acc_jet,rvel,bfield,ntot,gamax,r,sigsh); 
        }    
      
        if (isPair) {
            //Estimate pair annihilation rate for thermal plasma at base
            if (k==0){
                annihilation(z,zcut,r,ntot,eltemp,rate_pa,n_pa);
            } else{
                annihilation(z,zcut,r,ntot,mjteff,rate_pa,n_pa);
            }
        }
        ratio_ne = ntot/ntot0;

		//calculate specific parameters of each zone
		zonepars(njet,k,nz,z,r,delz,rvel,reff,hbb,tbb2,inclin,gamax0,gamax,nw,area,vol,gamv2,gamv,beta,theff,
		tbbeff,gshift,dopfac); 
		
		if (infosw == 1) {
			if (!isShock && mxsw == 1){                
				cout << "Z: " << z/r_g << " B: " << bfield << " n: " << ntot << " gmin: " << emin/emerg 
				<< endl;
			} else if (mxsw < 1) {
				cout << "Z: " << z/r_g << " B: " << bfield << " n: " << ntot << " gmin: " << emin/emerg 
				<< " gbreak: " << ebreak/emerg << " gmax: " << emax/emerg << endl;
			} else {
				cout << "Z: " << z/r_g << " B: " << bfield << " n: " << ntot << " gmin: " << emin/emerg
				<< " gbreak: " << gebreak << " gmax: " << gemax << endl;
			}
			cout << "velocity: " << rvel << " Doppler factor: " << 1./(sqrt(1.+rvel*rvel)*(1.-beta*
			cos(inclin))) << " radius: " << r/r_g << " opening angle: " << 57.3*2.*atan(r/z) << endl << endl;
		}

        //If MIXED electron distribution: in this case we must shift down the maxwellian component only, tag
        //on the power law as before (with an updated break energy), and renormalise everything
        if(mxsw != 1) {      
			// Udating break energy; Very soon this becomes NaN, but it was the same before: CHECK IT!
            update_ebreak(nz,k,z,r,h0,brk,bfield,beta,ebreak,gebreak);

            //Calculating the normalizations of pl and maxwell component
            elnorm(0,mxsw,ntot,pspec,plfrac,eltemp,gshift,emax,ebreak,mjteff,emin,cnorm,enorm,betat);
                               
            //Renormalizing the pl and thermal components, returning other arrays 
            pl_and_th_comp(isShock,thmfrac,nw,nelec,mxsw,rdlgen,thmbase,ebreak,emin,emax,gshift,gshock,mjteff,
            ratio_ne,ntot,cnorm,enorm,pspec,pltrm,etemp,dtemp,plcomp,thmcomp,bete);
                
            //Counting lost electrons and renormalising accordingly
            lost_ele(nw,nelec,etemp,dtemp,ntot,einc,elen,eled,lelec,ledens); 
        } //end of if on mxsw

        if(!isShock){

            /******************
             * NOZZLE + SHOCK *
             ******************/
            //Reading in distribution, shifting the energies down by mach**{-1/3},keep track of lost ones by
            //renormalizing first, then losing the first bins, then renormalizing
            read_and_shift_eldists(nelec,nw,ntot,ntot0,gshift,rdlgen,rdedn,elen,ledens,etemp,dtemp);
            
            //Counting the lost electrons and renormalising accordingly
            lost_ele(nw,nelec,etemp,dtemp,ntot,einc,elen,eled,lelec,ledens);
            
           	if(z>zsh){
           
            	/*********
             	 * SHOCK *
             	 *********/
            	//Here at the shock, in contrast to the previous version of code, we need to calculate the
            	//normalisations of a thermal + broken power-law electron distribution function, subject to
            	//whether the break energy lies within the limits of the power law distirbution (so, emin,
            	//emax, and ebreak). The thermal distirbution is read in from the electron distribution as it
            	//is at the shock. The minimum energy is set to the peak temperature of the MJ distribution
            	//(2.23*k*mjteff), and the maximum energy is calculated by equating the acceleration/cooling
            	//rates. The break energy is calculated by equating the synchrotron cooling time with the
            	//escape time of the electrons in the zone. Then the normalisation conditions for the power
            	//are set according to the 3 cases: 
            	//(i) ebreak < emin, (ii) ebreak > emax, (iii) emin < ebreak < emax. 
            	//The power law distribution itself is then set for whatever case is satisfied, and this
            	//gives both thermal and power law distributions at the shock, given by thmshock and plshock.
            	//The total distribution is the log of the sum. The arrays thmshock and plshock can then be 
            	//read in at every jet segment for increasing z over the loop.
                
				gshock  = gshift;

				//Shock heating
				eltemp = heat*eltemp;
				if (infosw ==1) {
					cout << "Electron temperature after heating: " << 2.23*eltemp << endl  << endl;
				}
				if (2.23*eltemp > 1.e12 && inclin >= 10) {//sometimes high heat crashes the code
					eltemp = 5.e11;
					cout << "Maxed out shock heating at 5x10^11!" << endl; 				
				}					
                //Estimating energy density for first order Compton
                ucompton(nz,njet,nsyn,k,dist,r,nusyn,dopfac,synabs,ucom);
               
               	//Adding in multicolor disk contribution plus component from second BB, w/ effective area 
               	//pi*reff**2                
                bbdisk_comp(disksw,bbsw,z,reff2,hbb,tin,rin,rout,gamv,tbbeff,ucom,uphdil,uphdil2);
               
				//Calculating the energy breaks; if normalization of PL at emin greater than thermal then	
				//below emin, fix maxwellian which will require renormalizing again
                maxene(z,r,ucom,fsc,bfield,beta,delz,brk,emax,gemax,ebreak,gebreak,bemax);            
                einc	= (log10(emax)-lelec[0])/nelec;
            	
            	//Calculate normalizations; if the flag is =1 the normalization is calculated with pspec, not
				//mxsw
                elnorm(1,mxsw,ntot,pspec,plfrac,eltemp,gshift,emax,ebreak,mjteff,emin,cnorm,enorm,betat);
	
                //Normalizing PL/Th components at shock
                norm_at_shock(thmfrac,heat,pspec,cnorm,enorm,betat,ebreak,emax,emin,einc,nelec,lelec,ledens,
                etemp,dtemp,thmshock,plshock,elen,ntot,eled,shden,shelen,bete);
                
                shocknd	= ntot;
                shockr	= r;
                shockgb	= rvel;
    
				if(isVerbose){
					sum1    = 0.;
                	sum2    = 0.;
                	for(i=0; i<nelec-1; i++){
                		tmp1	= (elen[i+1] - elen[i])*0.5*(eled[i]+eled[i+1]);
                    	tmp2	= (elen[i+1] - elen[i])*0.5*(pow(10, ledens[i]) + pow(10, ledens[i+1]));
                    	sum1	+= tmp1;
                    	sum2	+= tmp2;
					}				
                    cout << "THIS SHOULD ONLY APPEAR ONCE FOR z =" << z << " (k = " << k << ")" << endl;
                    cout << left << setw(15) << "ACCELPARS:" << setw(15) << r << setw(15) << z << setw(15) << 
                    dopfac[0*nz+k] << setw(15) << dopfac[1*nz+k] << setw(15) << bfield << setw(15) << ntot << 
                    endl;
                    cout << left << setw(15) << "ucom:" << setw(15) << ucom << endl;
                    cout << left << setw(15) << "ucom:" << setw(15) << ucom << setw(15) << uphdil << setw(15) 
                    << uphdil2 << endl;
                    cout << "EQUIP at Shock: " << bfield*bfield/(8.*pi*sum1) << endl;
                    cout << "NTOT after Shock: " << sum2 << " and " << ntot << endl;
				}
                               
                isShock = true;                
            }//END if z>zsh
    }//END nozzle + shock
    else if (z >= zsh && isShock){
            
            /****************
             * BEYOND SHOCK *
             ****************/            
            //Here we want to read in the distribution as it was at the shock, and then read in the thermal/
            //power law parts of the distribution separately (thmcomp[i] and plcomp[i]). Then we essentially 
            //repaeat the process completed at the shock, calculating the new min/max/break energies of the 
            //power law distribution according to the now gshifted thermal distribution. In this way we 
            //re-calculate the power law distribution as if it were fully replenished from the thermal 
            //distribution        
            
            //Read stored distribution directly now
            for(i=0; i < nelec; i++){
                etemp[i]= pow(10.,shelen[i]);
                dtemp[i]= shden[i];
                if(thmshock[i] > 0){//after elmax (initial max energy of thermal distribution), 
                					//thermal term = 0, thus cannot take logs....only take logs of values > 0, 
                					//set all others = 0.
                    thmcomp[i] = log10(thmshock[i]);
                }
                else{
                    thmcomp[i] = 0;
                }
                plcomp[i] = log10(plshock[i]); //loading in power law component from shock calculation
            }           
            
            //Shifting number density down from shock value, not from ntot0 anymore
            ntot	= shocknd*pow(shockr/r,2)*(shockgb/rvel);
        
            //Repeating the process completed at shock for every jet segment: determine ub/ucom, max energy, 
			//break energy, and shift distribution down
        
            //Estimating energy density for first order Compton
            ucompton(nz,njet,nsyn,k,dist,r,nusyn,dopfac,synabs,ucom);
        
            //Adding in multicolor disk contribution plus component from second BB, w/ effective area 
            //pi*reff**2   
            bbdisk_comp(disksw,bbsw,z,reff2,hbb,tin,rin,rout,gamv,tbbeff,ucom,uphdil,uphdil2);

            //Calculating the energy breaks; if normalization of PL at emin greater than thermal then
            //below emin, fix maxwellian which will require renormalizing again
            maxene(z,r,ucom,fsc,bfield,beta,delz,brk,emax,gemax,ebreak,gebreak,bemax);
            
            //Calculate normalizations without accounting for lost electrons; if the flag is =1 the 
			//normalization is calculated with pspec, not mxsw
            elnorm(1,mxsw,ntot,pspec,plfrac,eltemp,gshift,emax,ebreak,mjteff,emin,cnorm,enorm,betat);
            
            //Returning the pl and thermal components, plus other arrays
            pl_and_th_comp(isShock,thmfrac,nw,nelec,mxsw,rdlgen,thmbase,ebreak,emin,emax,gshift,gshock,mjteff,
            ratio_ne,ntot,cnorm,enorm,pspec,pltrm,etemp,dtemp,plcomp,thmcomp,bete);
             
            //Counting lost electrons and renormalizing accordingly
            lost_ele(nw,nelec,etemp,dtemp,ntot,einc,elen,eled,lelec,ledens);        
    
            if(isVerbose){
                cout << left << setw(15) << "ACCELPARS:" << setw(15) << r << setw(15) << z << setw(15) << 
                dopfac[0*nz+k] << setw(15) << dopfac[1*nz+k] << setw(15) << bfield << setw(15) << ntot << endl;
                cout << left << setw(15) << "ucom:" << setw(15) << ucom << setw(15) << uphdil << setw(15) << 
                uphdil2 << endl;
                
                //Double checking equipartition
                sum1    = 0.;
                for(i=0; i<nelec-1; i++){
                    tmp1	= (elen[i+1]-elen[i])*0.5*(eled[i]+eled[i+1]);
                    sum1	+=tmp1;
                }
                ub=bfield*bfield/(8.*pi);
                cout << "CHECK equipartition: " << ub/sum1 << endl;
            }
            
        }//END elseif (z>=zsh && isShock)
        else{
            cerr << "***Error: incorrect isShock and zsh values" << endl;
            exit(1);
        }
        /***************************
         * NOZZLE + SHOCK + BEYOND *
         *       RADIATION         *
         ***************************/
        //Back to common actions to all jet

        //Checking if the particles in the the jet are still relativistic
        if(ntot <= 0){
            for(i=0; i<nsyn; i++){
                for(m=0; m<njet; m++){
                    nusyn[(m*nz+k)*nsyn+i]= log10(nurad[i]*dopfac[m*nz+k]);
                    synabs[(m*nz+k)*nsyn+i]= -100;
                }
            }
            cerr << "ERROR: ***Ran out of relativistic particles!" << endl;
            exit(1);
        }

        //Defining min/max energy for integrals
        elenmn   = elen[0];
        elenmx   = elen[nelec-1];

        gsl_interp_accel *acc_eldis1    = gsl_interp_accel_alloc();
        gsl_spline *spline_eldis1       = gsl_spline_alloc(gsl_interp_cspline, nelec);
        gsl_spline_init(spline_eldis1, lelec, ledens, nelec);

        //Differentiating electron density for radiation calculations
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
        
        //The factor of 4.\times \pi in synabs is because an isotropic source is assumed.
        //synint returns units of ergs/cm^2/s/Hz for absorbed flux, so we need to convert to Jansky        
        synchrotron(isBreaknjetnsyn,nsyn,njet,nz,k,bfield,delz,inclin,dist,r,elenmn,elenmx,spline_syn,acc_syn,
		spline_eldis1,acc_eldis1,spline_derivs,acc_derivs,nurad,dopfac,nphot,nusyn,synabs);

        //When z>zcut, the Compton is negligible, so we cut the calculations here to speed up the code
        if(z>zcut){
            gsl_spline_free(spline_eldis1), gsl_interp_accel_free(acc_eldis1);
            gsl_spline_free(spline_derivs), gsl_interp_accel_free(acc_derivs);
            continue; //start new iteration of BIG LOOP
        }

        //spline_dbl sets up a photon density (erg/s/Hz) for passing to Compton routines, so we need to divide
        //by herg*cee*energy(erg)*pi*r**2 (erg**2*s**2*cm**3*Hz*s**-1) to get #/cm**3/erg. 
		//In order to add the diluted blackbody from the disk we take the max energy to be based based on tin
        //Both contributions return arrays phodis(#/cm**3/erg) and ephot(erg), which are used by the compton 
        //integrals as seed photons. Here we define the integration boundaries for the first Compton order
		//ephmax and ephmin
        seed_sync_and_disk_phtns(isVerbose,disksw,bbsw,tin,rin,rout,gamv,nsyn,snumax,z,r,reff2,hbb,
        tbbeff,nphot,nurad,nubb,energ,phodis,ephot,ephmax,ephmin);
        
        if(isVerbose){
            //Checking energy densities            
            for(i=0; i<nelec; i++){
                ledens[i] = pow(10,ledens[i]);
            }
            avelec = average_ene(0,nelec,ledens,elen,eled);
            cout << "Av electron: " << avelec << " and max/min energies: " << elenmx << "/" << elenmn << endl;
            avphot = average_ene(1,nsyn,phodis,ephot,phodis);
            cout << "Av photon: " << avphot << " and max/min energies: " << ephmax << "/" << ephmin << endl;
            for(i=0; i<nelec; i++){
                ledens[i] = log10(ledens[i]);
            }
        }

        //Interpolating the first seed photon distribution for Compton scattering
        gsl_interp_accel *acc_com2      = gsl_interp_accel_alloc();
        gsl_spline *spline_com2         = gsl_spline_alloc(gsl_interp_cspline, nsyn);
        gsl_spline_init(spline_com2,ephot,phodis,nsyn);
        
        //Chiara's single/multiple Compton code
        MICompton(isSSC,Niter,ncom,njet,nz,nzdum,k,r,ntot,vol,dist,dopfac,spline_eldis1,acc_eldis1,
        spline_com2,acc_com2,eltemp,elenmn,elenmx,ephmin,ephmax,ephxr,comspc,nucom,totcom);

        if (isPair) {
            //Dumping local IC-photon-number-density-spectrum at jet base to estimate pair production rate
            pproduction(ncom,r,ephxr,totcom,n_pp,rate_pp);
            cout << "log_10 [ rate_pp/ntot (1/s) ] ~ " << rate_pp/log10(ntot) << endl;
            cout << "log_10 [ rate_pa/ntot (1/s) ] ~ " << rate_pa/log10(ntot) << endl;          
            if (rate_pp > rate_pa && abs(rate_pp/log10(ntot)) > 1.0){
                cout << "Pairs are important @ k= " << k << endl;
            }
        }

        gsl_spline_free(spline_com2), gsl_interp_accel_free(acc_com2);
        gsl_spline_free(spline_eldis1), gsl_interp_accel_free(acc_eldis1);
        gsl_spline_free(spline_derivs), gsl_interp_accel_free(acc_derivs);
    }//END of BIG LOOP OVER Z

    //Calculating/interpolating some of the spectral components over the total array of frequencies
    spcomponents(infosw,plotsw,ne,njet,nz,nsyn,ncom,zed,zcut,zsh,ephxr,nusyn,synabs,nucom,comspc,cflx_array,
    nutot,complot,presyn,postsyn,fplot);
    
    for(i=0; i<(ne-1); i++){
        //Dumping all the spectral components into files/arrays
        write_plotFiles(plotsw,bbsw,disksw,tin,rin,rout,dist,inclin,bbf1,tbb2,nutot[i],fplot[i],complot[i],
        presyn[i],postsyn[i],bbplot[i]);

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
    delete[] gbx, delete[] gby, delete[] zed, delete[] nutot, delete[] nubb, delete[] nurad, delete[] energ,
    delete[] ephot, delete[] ephxr; delete[] dopfac, delete[] rdend, delete[] rdedn, delete[] rdlgen; 
    delete[] rden, delete[] lelec, delete[] elen, delete[] ledens, delete[] etemp, delete[] dtemp; 
    delete[] eled, delete[] nusyn, delete[] nucom, delete[] comspc, delete[] totcom, delete[] synabs; 
    delete[] shelen, delete[] shden, delete[] thmshock, delete[] thmcomp, delete[] thmbase, delete[] plshock;
    delete[] plcomp, delete[] plbase, delete[] tstrm, delete[] drtrm, delete[] phodis, delete[] nphot;
    delete[] jetu, delete[] disku, delete[] complot, delete[] cflx_array, delete[] presyn, delete[] postsyn;
    delete[] bbplot, delete[] total, delete[] fplot;   

    gsl_spline_free(spline_syn), gsl_interp_accel_free(acc_syn);
    gsl_spline_free(spline_jet), gsl_interp_accel_free(acc_jet);
} 
