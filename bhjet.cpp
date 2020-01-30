#include <fstream>
#include <stdarg.h>

#include "bhjet.hh"

using namespace std;

void jetmain(double *ear,int ne,double *param,double *photeng,double *photspec) {

	//STEP 1: VARIABLE/OBJECT DEFINITIONS
	//----------------------------------------------------------------------------------------------
	
	bool IsCounterjet = true;					//flag to decide whether to include a counterjet	
	bool IsShock = false;						//flag to check shock heating

	int nz = 80;								//total number of zones
	int nel = 70;
	int syn_res = 10;							//number of bins per decade in synch frequency;
	int com_res = 6;							//number of bins per decade in compton frequency;
	int nsyn,ncom;								//number of bins in synch/compton frequency;
	int Niter = 1;								//number of IC interactions
	int npsw;									//switch to define number of protons calculations
												//0: no protons
												//1: Up = Ue+Ub
												//2: ne = np	
												//pair only jet not implemented yet.
	
	double Mbh;									//black hole mass in solar masses
	double Eddlum;								//black hole Eddington luminosity
	double Rg;									//black hole  gravitational radius in cm	
	double theta;								//viewing angle	
	double dist;								//distance in kpc	
	double redsh;								//source redshift
	double jetrat;								//injected power in Eddington units	
	double zmin;								//jet launching point
	double r0;									//initial jet radius in rg
	double h;									//jet nozzle aspect ratio
	double zacc;								//end of bulk magnetic acceleration region in rg
	double zdiss;								//dissipation/nonthermal particle injection region in rg	
	double zmax;								//distance from bh up to which calculation continues
	double Te;									//temperature in kev, converted to erg	
	double plfrac;								//percentage of nonthermal particles post dissipation region
	double pldist;								//parameter to change plfrac over distance
	double pspec;								//slope of nonthermal distribution
	double heat;								//shock heating paramter
	double betaeff;								//effective expansion velocity used to set adiabatic cooling
	double fsc;									//particle acceleration timescale parameter		
	double pbeta;								//plasma beta in the jet base
	double sigmaf;								//final sigma when using magnetic acceleration
	double Tin;									//disk inner temperature in keV/luminosity in Eddington units
	double Rin;									//disk inner radius in rg
	double Rout;								//disk outer radius in rg
	//double Tcor;								//corona temperature in kev
	//double taucor;							//corona optical depth
	//double Rcor;								//corona effective radius/normalization, in units of Rg
	//double cordens;							//corona lepton number density
	double compar1;								//external inverse Compton parameters; different meanings
	double compar2;								//depending on the value of compsw
	double compar3;	
	double compsw;								//switch to activate different external Compton fields
	double velsw;								//velocity profile parameter
	int infosw;									//switch to print info
	
	double z;									//distance along the jet axis
	double tshift;								//temperature shift from initial value due to ad. cooling	
	double Urad;								//estimate of total radiation energy density in each zone
	double Ubb1,Ubb2;							//estimate of comoving energy density of external photons
	double gmin,gmax;							//minimum/maximum Lorentz factors over which to integrate
	double syn_min,syn_max;						//interval for synchrotron calculation in each zone
	double com_min,com_max;						//interval for inverse Compton calculation in each zone

	double *syn_en;								//sychrotron energy array for jet+counterjet summed
	double *syn_lum;							//synchrotron luminosity array for jet+counterjet summed
	double *com_en;								//compton energy array for jet+counterjet summed
	double *com_lum;							//compton luminosity array for jet+counterjet summed
	double *tot_en = new double[ne];			//energy arrray for sum of all zones and/or components
	double *tot_syn_pre = new double[ne];		//specific synchrotron luminosity arrays for all zones 
	double *tot_syn_post = new double[ne];		//pre/post particle	acceleration
	double *tot_com_pre = new double[ne];		//same as above but for the inverse Compton part
	double *tot_com_post = new double[ne];
	double *tot_lum = new double[ne];			//specific luminosity array for sum of all components	
	
	ofstream Numdensfile;						//ofstream plot file for particle distribution
	ofstream Presyn,Postsyn,Syn_zones;			//ofstream plot files for synchrotron emission
	ofstream Precom,Postcom,Com_zones;			//same as above but for inverse Comtpon
	ofstream Diskfile,Corfile,BBfile;			//same as above but for disk/corona/blackbody
	ofstream Totfile;							//same as above but for total model emission
	
	grid_pars grid;								//structure with grid parameters
	jet_dynpars jet_dyn;						//structure with jet dynamical parameters
	jet_enpars nozzle_ener;						//structure with jet energetic parameters
	zone_pars zone;								//strucutre with parameters of each individual zone
	com_pars agn_com;							//structure with parameters for inverse Compton fields in AGN
	
	//splines for jet acceleration
	gsl_interp_accel *acc_speed = gsl_interp_accel_alloc();
    gsl_spline *spline_speed = gsl_spline_alloc(gsl_interp_akima,54);
    
    //splines for electron distribution
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);

    gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
	gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  

	//STEP 2: PARAMETER/OBJECT INITIALIZATION
	//----------------------------------------------------------------------------------------------
		
	//read in parameter values
	Mbh = param[0];
	Eddlum = 1.25e38*Mbh;				
	Rg = gconst*Mbh*msun/(cee*cee);
	theta = param[1];		
	dist = param[2]*kpc;	
	redsh = param[3];
	jetrat = param[4]*Eddlum;
	r0 = param[5]*Rg;
	h = param[6]*r0;
	zdiss = param[7]*Rg;
	zacc = param[8]*Rg;	
	zmax = param[9]*Rg;
	Te = param[10]*kboltz_kev2erg;
	plfrac = param[11];
	pldist = param[12];
	pspec = param[13];
	heat = param[14];
	betaeff = param[15];
	fsc = param[16];
	pbeta = param[17];
	sigmaf = param[18];	
	Tin = param[19];
	Rin = param[20]*Rg;
	Rout = param[21]*Rg;
	//Tcor = param[20]*kboltz_kev2erg;
	//taucor = param[21];
	//Rcor = param[22]*Rg;
	//cordens = taucor/(Rcor*sigtom);
	compar1 = param[22];
	compar2 = param[23];
	compar3 = param[24];
	compsw = param[25];
	velsw = param[26];
	infosw = param[27];
	
	//initialize total energy/luminosity arrays
    for(int i=0;i<ne;i++){
       	tot_en[i] = (ear[i] + (ear[i+1]-ear[i])/2.)*herg/hkev;	
       	tot_syn_pre[i] = 1.e-50;
       	tot_syn_post[i] = 1.e-50;
       	tot_com_pre[i] = 1.e-50;
       	tot_com_post[i] = 1.e-50;
       	tot_lum[i] = 1.e-50;	
	}
	
	npsw = 1;	
	zmin = 2.*Rg;
	
	//Initialize disk+corona classes/objects
	//TODO sort out implementation of disk temperature/luminosity
	ShSDisk Disk(1,Mbh,Tin,Rin,Rout);
	//Thermal cor_elec(nel,1,Tcor);  
    //ncom = int(log10(1.e21)-log10(1.e15))*com_res*3;
    //Compton Corona(ncom,50,Niter,1.e15,1.e21,emgm);
    
	//Initialize external photon classes	
	BBody BLR(0,0);
	BBody Torus(0,0);
	BBody BlackBody(0,0);
	
	if(infosw>=1){			
		clean_file("Output/Presyn.dat",1);
		clean_file("Output/Postsyn.dat",1);
		clean_file("Output/Precom.dat",1);
		clean_file("Output/Postcom.dat",1);
		clean_file("Output/Disk.dat",1);
		clean_file("Output/BB.dat",1);
		//clean_file("Output/Corona.dat",1);
		clean_file("Output/Total.dat",1);	
	}
	if(infosw>=2){	
		clean_file("Output/Numdens.dat",0);
		clean_file("Output/Cyclosyn_zones.dat",1);
		clean_file("Output/Compton_zones.dat",1);			
	}
	
	//OPTIONAL STEP 3: DISK/CORONA/EXTERNAL PHOTON CALCULATIONS
	//----------------------------------------------------------------------------------------------
	
	//calculate disk+corona spectrum if desired
	//The disk is disabled by setting Rin<Rout, the corona by setting tau <= 0
	if(Rin<Rout){		
		Disk.set_inclination(theta);
		Disk.disk_spectrum();			
		/*if(taucor>0){
			cor_elec.set_p();    	
			cor_elec.set_norm(cordens);
			cor_elec.set_ndens();
				
			gmin = cor_elec.get_gamma()[0];
			gmax = cor_elec.get_gamma()[nel-1];
				   
			gsl_spline_init(spline_eldis,cor_elec.get_gamma(),cor_elec.get_gdens(),nel);
		
			//calculate corona spectrum, scattering only disk photons	
			Corona.set_beaming(theta,0.,1.);
			Corona.set_geometry(1,Rcor,0);
			Corona.set_tau(cordens,Rcor,cor_elec.av_gamma());
			//Multiple scatters only if ypar and tau are large enough
			if(Corona.get_ypar() > 1.e-2 && Corona.get_tau() > 1.e-2){
				Corona.set_niter(15);    
			}
			Corona.shsdisk_seed(Disk.get_energ(),Disk.tin(),Rin,Rout,Disk.hdisk(),0);        
			Corona.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);		
			if (infosw>=3){
				cout << "Disk and corona parameters: " << endl;
				Disk.test();
				Corona.test();
				cout << "Power in the corona: " << pi*pow(Rcor,2.)*cordens*cor_elec.av_gamma()*emerg*cee;
				cout << endl << endl;
			}
			sum_ext(ncom,ne,Corona.get_energ_obs(),Corona.get_nphot_obs(),tot_en,tot_lum);	
		}*/
		sum_ext(50,ne,Disk.get_energ_obs(),Disk.get_nphot_obs(),tot_en,tot_lum);   	
		if(infosw >= 3) {
			Disk.test();
			cout << endl;
		}	
	}	
	
	//Depending on the value of compsw, we either include a) an extra homogeneous black body in every zone or
	//b) the radiation field of the broad line region/torus of a bright AGN
	if (compsw==1){
		BlackBody.set_temp_k(compar1);
		BlackBody.set_lum(compar2);
		Ubb1 = compar3;	
		BlackBody.set_inclination(0.);
		BlackBody.bb_spectrum();
		sum_ext(40,ne,BlackBody.get_energ(),BlackBody.get_nphot(),tot_en,tot_lum);
	}
	else if (compsw==2 && Rin<Rout){
		agn_photons_init(Disk.total_luminosity(),compar1,compar2,agn_com);
			
		BLR.set_temp_kev(agn_com.tblr);
		BLR.set_lum(agn_com.lblr);
		BLR.set_inclination(0.);
		BLR.bb_spectrum();			
		
		Torus.set_temp_k(agn_com.tdt);
		Torus.set_lum(agn_com.ldt);	
		Torus.set_inclination(0.);
		Torus.bb_spectrum();		
		if(infosw>=3){
			cout << "BLR radius in Rg: " << agn_com.rblr/Rg << endl;
			cout << "DT radius in Rg: " << agn_com.rdt/Rg << endl;
		}
		sum_ext(40,ne,Torus.get_energ(),Torus.get_nphot(),tot_en,tot_lum);	
	}	
		
	//STEP 4: JET BASE EQUIPARTITION CALCULATIONS AND SETUP
	//----------------------------------------------------------------------------------------------	
	
    //Dummy particle distribution, needed for average lorentz factor in equipartition function
    //The parameters for the method to set the momentum array are set to dummy values that result in a fully
    //thermal distribution
    Mixed dummy_elec(nel,1,Te,10.,0);
	dummy_elec.set_plfrac(0);
	dummy_elec.set_temp(Te);
	dummy_elec.set_p(0.,1.,0.0,r0,1.e-16);	
	dummy_elec.set_norm(1.);	
	dummy_elec.set_ndens();

    grid.nz = nz;
    grid.cut = 0;
    grid.zcut = 1.e3*Rg;
    
    jet_dyn.min = zmin;
    jet_dyn.max = zmax;
    jet_dyn.h0 = h+zmin;
    jet_dyn.r0 = r0;
    jet_dyn.acc = zacc;
    jet_dyn.beta0 = sqrt(4./3.*(4./3.-1.)/(4./3.+1.));	//set initial jet speed for relativistic fluid, g=4/3 
    jet_dyn.gam0 = 1./sqrt(1.-(pow(jet_dyn.beta0,2.)));	//set corresponding lorentz factor
    jet_dyn.gamf = velsw;
    jet_dyn.Rg = Rg;  
    
    nozzle_ener.pbeta = pbeta;
    nozzle_ener.Nj = jetrat;
   	nozzle_ener.sigf = sigmaf; 
	//set up jet velocity profile depending on choice of adiabatic,isothermal,magnetically dominated jet   
	//note: the adiabatic jet only runs correctly if the final temperature is above ~1kev, which means the
	//initial temperature has to be ~10^4 kev
    if(velsw==0){   
		velprof_ad(spline_speed);
		equipartition(IsCounterjet,npsw,jet_dyn,nozzle_ener);	
    } else if(velsw==1){
      	velprof_iso(spline_speed);
      	equipartition(IsCounterjet,npsw,jet_dyn,nozzle_ener);
    } else {
    	velprof_mag(jet_dyn,spline_speed);
    	equipartition(IsCounterjet,jetrat,jet_dyn,nozzle_ener);
    }	
	
	if(infosw>=3){
		cout << "Jet base parameters: " << endl;
		cout << "Pair content (ne/np): " << nozzle_ener.eta << endl;
		cout << "Initial magnetization: " << nozzle_ener.sig0 << endl;
		cout << "Particle average Lorenz factor: " << dummy_elec.av_gamma() << endl;
	}	
	
	//STEP 5: TOTAL JET CALCULATIONS, LOOPING OVER EACH SEGMENT OF THE JET
	//----------------------------------------------------------------------------------------------

	for(int i=0;i<nz;i++){
		//calculate dynamics/energetics in each zone
		jetgrid(i,grid,jet_dyn,zone.r,zone.delz,z);	
		if(velsw==0){
			adjetpars(z,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);	 
		} else if (velsw==1){
			isojetpars(z,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);
		} else {
			bljetpars(z,betaeff,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);
		}
		zone.delta = 1./(zone.gamma*(1.-zone.beta*cos(theta*pi/180.)));
		
		//Before the dissipation/particle acceleration region follow the temperature profile from the jetpars
		//function. Afterwards do the same, but include an extra dissipation term as a function of distance
		//pow(log10(zdiss)/log10(z),pldist) in order to obtain an inverted radio spectrum. The max() use is
		//to because the particle distrubution code currently crashes for temperatures below 1kev
		if(z<zdiss){
			zone.eltemp = tshift*Te;
		} else {
			zone.eltemp = max(tshift*Te*pow(log10(zdiss)/log10(z),pldist),kboltz_kev2erg);
		}
				
		//if compsw==1 the boosting is the same in all zones, if compsw==2 we boost either blr,torus or both
		//depending on the distance from the BH of the jet region/seed photon production
		if(compsw==1){			
			Urad = pow(zone.gamma,2.)*Ubb1;
		} else if (compsw==2 && Rin<Rout){
			zone_agn_phfields(z,zone,Ubb1,Ubb2,agn_com);
			Urad = agn_com.urad_total;
		} else {
			Urad = 0.;
		}

		//calculate particle distribution in each zone
		if(z<zdiss){		
			Thermal th_lep(nel,1,zone.eltemp);
			th_lep.set_p();
			th_lep.set_norm(zone.lepdens);
			th_lep.set_ndens();
			
			gmin = th_lep.get_gamma()[0];
			gmax = th_lep.get_gamma()[nel-1];
			
			zone.avgammasq = pow(th_lep.av_gamma(),2.);
					
   	 		gsl_spline_init(spline_eldis,th_lep.get_gamma(),th_lep.get_gdens(),nel);
    		gsl_spline_init(spline_deriv,th_lep.get_gamma(),th_lep.get_gdens_diff(),nel); 
  		
			if (infosw >=2){
				plot_write(nel,th_lep.get_p(),th_lep.get_gamma(),th_lep.get_pdens(),th_lep.get_gdens(),
				"Output/Numdens.dat");
			}
		} else {		
			if (IsShock==false){
				Te = heat*Te;
				IsShock = true;
			}					
			Mixed acc_lep(nel,1,zone.eltemp,pspec,1);
			//Mixed acc_lep(nel,1,tshift*Te,pspec,1);
			acc_lep.set_plfrac(plfrac*pow(log10(zdiss)/log10(z),pldist));
			//some combination of plfrac*pow(log10(zdiss)/log10(z),pldist.) or
			//tshift*Te*pow(log10(zdiss)/log10(z),pldist.)  might work for inverting spectra
			//doing both is more consistent and makes physical sense

			//if fsc < 10 it's the acceleration efficiency, else it's the desired maximum lorentz factor
			if (fsc<10.){
				acc_lep.set_p(Urad,zone.bfield,betaeff,zone.r,fsc);
			} else{
				acc_lep.set_p(fsc);
			}	
			
			acc_lep.set_norm(zone.lepdens);	
			acc_lep.set_ndens();
			acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,betaeff);
			//Note: this assumes sycnhrotron+ec cooling only, ssc+disk cooling negligible
		
			gmin = acc_lep.get_gamma()[0];
			gmax = acc_lep.get_gamma()[nel-1];
			
			zone.avgammasq = pow(acc_lep.av_gamma(),2.);
				
			gsl_spline_init(spline_eldis,acc_lep.get_gamma(),acc_lep.get_gdens(),nel);
    		gsl_spline_init(spline_deriv,acc_lep.get_gamma(),acc_lep.get_gdens_diff(),nel); 
    			 	
			if(infosw>=2){
				plot_write(nel,acc_lep.get_p(),acc_lep.get_gamma(),acc_lep.get_pdens(),acc_lep.get_gdens(),
				"Output/Numdens.dat");
			}
		}
		
		//Note: the energy density below assumes only cold protons!
		if (infosw>=4){
			double Up,Ue,Ub;
			Ue = sqrt(zone.avgammasq)*zone.lepdens*emerg;
			Up = (zone.lepdens/nozzle_ener.eta)*pmgm*pow(cee,2.);
			Ub = pow(zone.bfield,2.)/(8.*pi);
			cout << endl <<  "Jetpars; Bfield: " << zone.bfield <<  ", Lepton ndens: " << zone.lepdens 
				 << ", speed: " << zone.gamma << ", delta: " << zone.delta << endl;
			cout << "tshift: " << tshift << ", Temperature in keV: " << zone.eltemp/kboltz_kev2erg << endl;
			cout << "Grid; R: " << zone.r/Rg << ", delz: " << zone.delz/Rg << ", z: " << z/Rg << ", z+delz: " 
			     << (zone.delz+z)/Rg << endl;
			cout << "Equipartition check; Sigma: " << 2.*Ub/Up << " Ue/Ub: " << Ue/Ub << endl;
		}
	
		//calculate emission of each zone		
		//note: the syn_en array is used for the seed photon fields in the IC part, so it needs to include
		//both the black body and disk part. This is why the maximum frequency is taken as the maximum of the
		//two scale frequencies.
		syn_min = 0.1*pow(gmin,2.)*charg*zone.bfield/(2.*pi*emgm*cee);
		if(Rin<Rout){
			syn_max = max(50.*pow(gmax,2.)*charg*zone.bfield/(2.*pi*emgm*cee),20.*Disk.tin()*kboltz/herg);
		}
		else {
			syn_max = 50.*pow(gmax,2.)*charg*zone.bfield/(2.*pi*emgm*cee);
		}
		nsyn = int(log10(syn_max)-log10(syn_min))*syn_res;
		syn_en = new double[nsyn];			
		syn_lum = new double[nsyn];	
		Cyclosyn Syncro(nsyn,syn_min,syn_max,emgm);	
		
		com_min = 0.1*Syncro.nu_syn();
		com_max = ear[ne-1]/hkev;
		ncom = int(log10(com_max)-log10(com_min))*com_res;
		com_en = new double[ncom];			
		com_lum = new double[ncom];		
		Compton InvCompton(ncom,nsyn,Niter,com_min,com_max,emgm);		
		for (int k=0;k<ncom;k++){
			com_lum[k] = 0;
			com_en[k] = InvCompton.get_energ()[k];
		}
		//Note: initializing these two arrays is not stricly necessary for running the code correctly, but it
		//avoids messing up the color scale for the zone plots done by Plots.py
		
		//calculate cyclosynchrotron spectrum
		//Set up the calculation by reading in magnetic field,beaming,volume,counterjet presence	
		Syncro.set_bfield(zone.bfield);
		Syncro.set_beaming(theta,zone.beta,zone.delta);
		Syncro.set_geometry(0,zone.r,zone.delz);
		Syncro.set_counterjet(IsCounterjet);
    	Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);
    	
		if (infosw>=5){
			Syncro.test();
		}
		//Sum counterjet if present, then save emission from the zone in syn_en/syn_lum
		if (IsCounterjet==true){
			sum_counterjet(nsyn,Syncro.get_energ_obs(),Syncro.get_nphot_obs(),syn_en,syn_lum);			
		} else {
			for(int k=0;k<nsyn;k++){
				syn_en[k] = Syncro.get_energ_obs()[k];
				syn_lum[k] = Syncro.get_nphot_obs()[k];
			}
		}
		//Include zone's emission to the pre/post particle acceleration spectrum
		if(z<zdiss){
			sum_zones(nsyn,ne,syn_en,syn_lum,tot_en,tot_syn_pre);
		} else {
			sum_zones(nsyn,ne,syn_en,syn_lum,tot_en,tot_syn_post);
		}			
		//calculate inverse Compton spectrum, if it's expected to be bright enough	
		if (Compton_check(IsShock,i,Mbh,Urad,zone) == true){
		//if(z>zmax){
			//Set up the calculation by reading in/calculating beaming,volume,counterjet presence,tau
			InvCompton.set_beaming(theta,zone.beta,zone.delta);
			InvCompton.set_geometry(0,zone.r,zone.delz);
			InvCompton.set_counterjet(IsCounterjet);	
			InvCompton.set_tau(zone.lepdens,zone.r,Te/emerg);
			//Multiple scatters only if ypar and tau are large enough
			if(InvCompton.get_ypar() > 1.e-2 && InvCompton.get_tau() > 1.e-2){
				InvCompton.set_niter(15);		
			}
			//Cyclosynchrotron photons are always considered in the scattering						
			InvCompton.cyclosyn_seed(Syncro.get_energ(),Syncro.get_nphot());
			//Disk photons are included only if the disk is present
			if(Rin<Rout){
				InvCompton.shsdisk_seed(Syncro.get_energ(),Disk.tin(),Rin,Rout,Disk.hdisk(),z+zone.delz/2.);
			}
			//Black body photons included only if compsw==1
			if(compsw==1){						
				InvCompton.bb_seed(Syncro.get_energ(),Ubb1,zone.delta*BlackBody.temp_kev());
			}
			//AGN photon fields photons are considered only if disk is present and compsw==2
			if(compsw==2 && Rin<Rout){
				InvCompton.bb_seed(Syncro.get_energ(),Ubb1,zone.delta*BLR.temp_kev());
				InvCompton.bb_seed(Syncro.get_energ(),Ubb2,zone.delta*Torus.temp_kev());
			}
			//Calculate the spectrum with whichever fields have been invoked		
			InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
	
			if (infosw>=5){
				InvCompton.test();
			}
			//Sum counterjet if present, then save emission from the zone in com_en/com_lum		
			if (IsCounterjet==true){
				sum_counterjet(ncom,InvCompton.get_energ_obs(),InvCompton.get_nphot_obs(),com_en,com_lum);
			} else {
				for(int k=0;k<ncom;k++){
					com_en[k] = InvCompton.get_energ_obs()[k];
					com_lum[k] = InvCompton.get_nphot_obs()[k];
				}			
			}
			//Include zone's emission to the pre/post particle acceleration spectrum
			if(z<zdiss){
				sum_zones(ncom,ne,com_en,com_lum,tot_en,tot_com_pre);
			} else {
				sum_zones(ncom,ne,com_en,com_lum,tot_en,tot_com_post);
			}
		} else if (infosw>=4){
			cout << "Out of the Comptonization region" << endl;
		}					
		if(infosw>=2){
			plot_write(nsyn,syn_en,syn_lum,"Output/Cyclosyn_zones.dat",dist,redsh);
			plot_write(ncom,com_en,com_lum,"Output/Compton_zones.dat",dist,redsh);
		}			
				
		delete[] syn_en,delete[] syn_lum;
		delete[] com_en,delete[] com_lum;			
	} 

	//FINAL STEP: SUM JET COMPONENTS TO TOTAL OUTPUT, WRITE/CLOSE PLOT FILES, FREE MEMORY
	//----------------------------------------------------------------------------------------------
	//include redshfit stuff here
	for(int k=0;k<ne;k++){
		tot_lum[k] = (tot_lum[k]+tot_syn_pre[k]+tot_syn_post[k]+tot_com_pre[k]+tot_com_post[k]);
		photeng[k] = log10(tot_en[k]/herg);
	}
	output_spectrum(ne,tot_en,tot_lum,photspec,redsh,dist);
		
	if(infosw>=1){
		plot_write(ne,tot_en,tot_syn_pre,"Output/Presyn.dat",dist,redsh);
		plot_write(ne,tot_en,tot_syn_post,"Output/Postsyn.dat",dist,redsh);
		plot_write(ne,tot_en,tot_com_pre,"Output/Precom.dat",dist,redsh);
		plot_write(ne,tot_en,tot_com_post,"Output/Postcom.dat",dist,redsh);
		plot_write(50,Disk.get_energ_obs(),Disk.get_nphot_obs(),"Output/Disk.dat",dist,redsh);
		//plot_write(ncom,Corona.get_energ_obs(),Corona.get_nphot_obs(),"Output/Corona.dat",dist,redsh);
		if(compsw==2){
			plot_write(40,Torus.get_energ(),Torus.get_nphot(),"Output/BB.dat",dist,redsh);
		} else {
			plot_write(40,BlackBody.get_energ(),BlackBody.get_nphot(),"Output/BB.dat",dist,redsh);
		}			
		plot_write(ne,tot_en,tot_lum,"Output/Total.dat",dist,redsh);
	}	
	
	gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
	gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);	
	gsl_spline_free(spline_speed), gsl_interp_accel_free(acc_speed);
	delete[] tot_en,delete[] tot_lum;
	delete[] tot_syn_pre,delete[] tot_syn_post;
	delete[] tot_com_pre,delete[] tot_com_post;
}
