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
    double plfrac_0;							//percentage of nonthermal particles at the dissipation region
    double pldist;								//parameter to change plfrac over distance
    double pspec;								//slope of nonthermal distribution
    double heat;								//shock heating paramter
    double betaeff;								//effective expansion velocity used to set adiabatic cooling
    double fsc;									//particle acceleration timescale parameter		
    double pbeta;								//plasma beta in the jet base
    double sig_acc;								//final sigma when using magnetic acceleration
    double Ldisk;								//disk luminosity in Eddington units
    double Rin;									//disk inner radius in rg
    double Rout;								//disk outer radius in rg
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
    gsl_spline *spline_speed = gsl_spline_alloc(gsl_interp_steffen,54);

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
    zdiss = param[6]*Rg;
    zacc = param[7]*Rg;	
    zmax = param[8]*Rg;
    Te = param[9]*kboltz_kev2erg;
    plfrac_0 = param[10];
    pldist = param[11];
    pspec = param[12];
    heat = param[13];
    betaeff = param[14];
    fsc = param[15];
    pbeta = param[16];
    sig_acc = param[17];	
    Ldisk = param[18];
    Rin = param[19]*Rg;
    Rout = param[20]*Rg;
    compar1 = param[21];
    compar2 = param[22];
    compar3 = param[23];
    compsw = param[24];
    velsw = param[25];
    infosw = param[26];
    
    if (infosw>=1) {
        param_write(param,"Output/Starting_pars.dat");
    }

    //initialize total energy/luminosity arrays
    for(int i=0;i<ne;i++){
       	tot_en[i] = (ear[i] + (ear[i+1]-ear[i])/2.)*herg/hkev;	
       	tot_syn_pre[i] = 1.;
       	tot_syn_post[i] = 1.;
       	tot_com_pre[i] = 1.;
       	tot_com_post[i] = 1.;
       	tot_lum[i] = 1.;	
    }

    npsw = 1;	
    zmin = 2.*Rg;

    //Initialize disk+external photon classes
    ShSDisk Disk;
    BBody BLR;
    BBody Torus;
    BBody BlackBody;

    if(infosw>=1){			
        clean_file("Output/Presyn.dat",2);
        clean_file("Output/Postsyn.dat",2);
        clean_file("Output/Precom.dat",2);
        clean_file("Output/Postcom.dat",2);
        clean_file("Output/Disk.dat",2);
        clean_file("Output/BB.dat",2);
        clean_file("Output/Total.dat",2);	
    }
    if(infosw>=2){	
        clean_file("Output/Numdens.dat",4);
        clean_file("Output/Cyclosyn_zones.dat",2);
        clean_file("Output/Compton_zones.dat",2);			
    }

    //STEP 3: DISK/EXTERNAL PHOTON CALCULATIONS

    //The disk is disabled by setting Rin<Rout; its contribution is summed to the total only if there are no 
    //AGN photon fields that reprocess part of the luminosity, otherwise it is done later
    if(Rin<Rout){
    	Disk.set_mbh(Mbh);
	    Disk.set_rin(Rin);
	    Disk.set_rout(Rout);
	    Disk.set_luminosity(Ldisk);		
	    Disk.set_inclination(theta);
	    Disk.disk_spectrum();
	    if (compsw != 2) {
	        sum_ext(50,ne,Disk.get_energy_obs(),Disk.get_nphot_obs(),tot_en,tot_lum);   	    
	    }			
        if(infosw >= 3) {
            Disk.test();
            cout << endl;
        }	
    }	

    //Depending on the value of compsw, we either include a) an extra homogeneous black body in every zone or
    //b) the radiation field of the broad line region/torus of a bright AGN. For bright AGN, the reprocessed
    //fraction of disk luminosity is removed from the observed disk luminosity
    if (compsw==1){
        BlackBody.set_temp_k(compar1);
        BlackBody.set_lum(compar2);
        Ubb1 = compar3;	
        BlackBody.bb_spectrum();
        sum_ext(40,ne,BlackBody.get_energy_obs(),BlackBody.get_nphot_obs(),tot_en,tot_lum);
    }
    else if (compsw==2 && Rin<Rout){
        agn_photons_init(Disk.total_luminosity(),compar1,compar2,agn_com);

        BLR.set_temp_kev(agn_com.tblr);
        BLR.set_lum(agn_com.lblr);
        BLR.bb_spectrum();			

        Torus.set_temp_k(agn_com.tdt);
        Torus.set_lum(agn_com.ldt);	
        Torus.bb_spectrum();
        
        Disk.cover_disk(compar1+compar2);		
        if(infosw>=3){
            cout << "BLR radius in Rg: " << agn_com.rblr/Rg << endl;
            cout << "DT radius in Rg: " << agn_com.rdt/Rg << endl;
        }
        sum_ext(40,ne,Torus.get_energy_obs(),Torus.get_nphot_obs(),tot_en,tot_lum);
        sum_ext(50,ne,Disk.get_energy_obs(),Disk.get_nphot_obs(),tot_en,tot_lum);   		
    }	
	
    //STEP 4: JET BASE EQUIPARTITION CALCULATIONS AND SETUP
    //----------------------------------------------------------------------------------------------	

    //Dummy particle distribution, needed for average lorentz factor in equipartition function
    //The parameters for the method to set the momentum array are set to dummy values that result in a fully
    //thermal distribution 
    Thermal dummy_elec(nel);
    dummy_elec.set_mass(emgm);
    dummy_elec.set_temp(Te);
    dummy_elec.set_p();	
    dummy_elec.set_norm(1.);	
    dummy_elec.set_ndens();

    grid.nz = nz;
    grid.cut = 0;
    grid.zcut = 1.e3*Rg;

    jet_dyn.min = zmin;
    jet_dyn.max = zmax;
    jet_dyn.h0 = 2.*r0+zmin;
    jet_dyn.r0 = r0;
    jet_dyn.acc = zacc;
    jet_dyn.beta0 = sqrt(4./3.*(4./3.-1.)/(4./3.+1.));	//set initial jet speed for relativistic fluid, g=4/3 
    jet_dyn.gam0 = 1./sqrt(1.-(pow(jet_dyn.beta0,2.)));	//set corresponding lorentz factor
    jet_dyn.gamf = velsw;
    jet_dyn.Rg = Rg;  

    nozzle_ener.pbeta = pbeta;
    nozzle_ener.Nj = jetrat;
    nozzle_ener.sig_acc = sig_acc; 
    nozzle_ener.av_gamma = dummy_elec.av_gamma();
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
        cout << "Jet nozzle ends at: " << jet_dyn.h0/Rg << " Rg" << endl ;
        cout << "Jet nozzle optical depth: " << jet_dyn.r0*nozzle_ener.lepdens*sigtom << endl << endl;
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

        //This is to avoid crashes due to low (sub 1 kev) particle temperatures
        if (z < zdiss) {
            zone.eltemp = max(tshift*Te,kboltz_kev2erg);    
        } else {
            zone.eltemp = max(tshift*Te*pow(log10(zdiss)/log10(z),pldist),kboltz_kev2erg);
        }
        
        
        //This is to evolve the fraction of non thermal particles along the jet, and change the distribution 
        //of non thermal particles appropriately
        if (z < zdiss) {
            zone.nth_frac = 0.;
        } else {
            zone.nth_frac = plfrac_0*pow(log10(zdiss)/log10(z),pldist);            
        }


        //Include the disk for radiative cooling if it's on;
        //if compsw==1 the boosting is the same in all zones, if compsw==2 we boost either blr,torus or both
        //depending on the distance from the BH of the jet region/seed photon production
        if (Rin < Rout) {
            double Rdisk = pow(Rin,2.)+pow(z,2.);
            double delta_disk, theta_disk;
            theta_disk = pi-atan(Rin/z);
            delta_disk = 1./(zone.gamma-zone.beta*cos(theta_disk));
            Urad = pow(delta_disk,2.)*Ldisk*Eddlum/(4.*pi*Rdisk*cee);
        } else {
            Urad = 0.;
        }
        
        if(compsw==1){			
            Urad = Urad + pow(zone.delta,2.)*Ubb1;
        } else if (compsw==2 && Rin<Rout){
            zone_agn_phfields(z,zone,Ubb1,Ubb2,agn_com);
            Urad = Urad + agn_com.urad_total;
        }          

        //calculate particle distribution in each zone
        if(zone.nth_frac == 0.){		
            Thermal th_lep(nel);
            th_lep.set_mass(emgm);
            th_lep.set_temp(zone.eltemp);
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
        } else if (zone.nth_frac < 0.5){ 
            if (IsShock==false){
                Te = heat*Te;
                IsShock = true;
            }					
            Mixed acc_lep(nel);
            acc_lep.set_mass(emgm);
            acc_lep.set_temp(zone.eltemp);
            acc_lep.set_pspec(pspec);
            acc_lep.set_plfrac(zone.nth_frac);
            
            //if fsc < 10 it's the acceleration efficiency, else it's the desired maximum lorentz factor
            if (fsc<10.){
                acc_lep.set_p(Urad,zone.bfield,betaeff,zone.r,fsc);
            } else{
                acc_lep.set_p(fsc);
            }	

            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,betaeff);
            //Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel-1];

            zone.avgammasq = pow(acc_lep.av_gamma(),2.);

            gsl_spline_init(spline_eldis,acc_lep.get_gamma(),acc_lep.get_gdens(),nel);
            gsl_spline_init(spline_deriv,acc_lep.get_gamma(),acc_lep.get_gdens_diff(),nel); 
             	
            if(infosw>=2){
                plot_write(nel,acc_lep.get_p(),acc_lep.get_gamma(),acc_lep.get_pdens(),acc_lep.get_gdens(),
                "Output/Numdens.dat");
            }
        } else if (zone.nth_frac < 1.) {
            if (IsShock==false){
                Te = heat*Te;
                zone.eltemp = max(tshift*Te*pow(log10(zdiss)/log10(z),pldist),kboltz_kev2erg);
                IsShock = true;
            }
            Thermal dummy_elec(nel);
            dummy_elec.set_mass(emgm);
            dummy_elec.set_temp(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pbrk = dummy_elec.av_p();
            
            Bknpower acc_lep(nel);
            acc_lep.set_mass(emgm);
            acc_lep.set_pspec1(-2.);
            acc_lep.set_pspec2(pspec);
                        
            if (fsc<10.){
                acc_lep.set_p(0.1*pbrk,pbrk,Urad,zone.bfield,betaeff,zone.r,fsc);
            } else{
                acc_lep.set_p(0.1*pbrk,pbrk,fsc);
            }	
            
            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,betaeff);
            //Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel-1];

            zone.avgammasq = pow(acc_lep.av_gamma(),2.);

            gsl_spline_init(spline_eldis,acc_lep.get_gamma(),acc_lep.get_gdens(),nel);
            gsl_spline_init(spline_deriv,acc_lep.get_gamma(),acc_lep.get_gdens_diff(),nel); 
             	
            if(infosw>=2){
                plot_write(nel,acc_lep.get_p(),acc_lep.get_gamma(),acc_lep.get_pdens(),acc_lep.get_gdens(),
                "Output/Numdens.dat");
            }
        } else if (zone.nth_frac == 1.) {
            if (IsShock==false){
                Te = heat*Te;
                zone.eltemp = max(tshift*Te*pow(log10(zdiss)/log10(z),pldist),kboltz_kev2erg);
                IsShock = true;
            }
            Thermal dummy_elec(nel);
            dummy_elec.set_mass(emgm);
            dummy_elec.set_temp(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pmin = dummy_elec.av_p();
            
            Powerlaw acc_lep(nel);
            acc_lep.set_mass(emgm);
            acc_lep.set_pspec(pspec);
                        
            if (fsc<10.){
                acc_lep.set_p(pmin,Urad,zone.bfield,betaeff,zone.r,fsc);
            } else{
                acc_lep.set_p(pmin,fsc);
            }	
            
            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,betaeff);
            //Note: this assumes ssc cooling is negligible

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

        //Note: the energy density below assumes only cold protons
        if (infosw>=5){
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
        Cyclosyn Syncro(nsyn);
        Syncro.set_frequency(syn_min,syn_max);	

        com_min = 0.1*Syncro.nu_syn();
        com_max = ear[ne-1]/hkev;
        ncom = int(log10(com_max)-log10(com_min))*com_res;
        com_en = new double[ncom];			
        com_lum = new double[ncom];		
        Compton InvCompton(ncom,nsyn);
        InvCompton.set_frequency(com_min,com_max);	
        
        if(infosw > 1) {
            for (int k=0;k<ncom;k++){
                com_lum[k] = 0;
                com_en[k] = InvCompton.get_energy()[k];
            }
        }	        
        //Note: initializing these two arrays is only done to plot each zone correctly using the colorscale in
        //Plot.py

        //calculate cyclosynchrotron spectrum
        //Set up the calculation by reading in magnetic field,beaming,volume,counterjet presence	
        Syncro.set_bfield(zone.bfield);
        Syncro.set_beaming(theta,zone.beta,zone.delta);
        Syncro.set_geometry("cylinder",zone.r,zone.delz);
        //IsCounterjet = true;
        Syncro.set_counterjet(IsCounterjet);
        Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);

        if (infosw>=5){
            Syncro.test();
        }
        //Sum counterjet if present, then save emission from the zone in syn_en/syn_lum
        if (IsCounterjet==true){
            sum_counterjet(nsyn,Syncro.get_energy_obs(),Syncro.get_nphot_obs(),syn_en,syn_lum);			
        } else {
            for(int k=0;k<nsyn;k++){
                syn_en[k] = Syncro.get_energy_obs()[k];
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
        if (Compton_check(IsShock,i,Mbh,Urad,velsw,zone) == true){
        //if(z>zmax){
            //Set up the calculation by reading in/calculating beaming,volume,counterjet presence,tau
            InvCompton.set_beaming(theta,zone.beta,zone.delta);
            InvCompton.set_geometry("cylinder",zone.r,zone.delz);
            InvCompton.set_counterjet(IsCounterjet);	
            InvCompton.set_tau(zone.lepdens,Te/emerg);
            //Multiple scatters only if ypar and tau are large enough
            if(InvCompton.get_ypar() > 1.e-2 && InvCompton.get_tau() > 5.e-2){
                InvCompton.set_niter(15);		
            }
            //Cyclosynchrotron photons are always considered in the scattering						
            InvCompton.cyclosyn_seed(Syncro.get_energy(),Syncro.get_nphot());
            //Disk photons are included only if the disk is present
            if(Rin<Rout){
                InvCompton.shsdisk_seed(Syncro.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),z+zone.delz/2.);
            }
            //Black body photons included only if compsw==1
            if(compsw==1){						
                InvCompton.bb_seed(Syncro.get_energy(),Ubb1,zone.delta*BlackBody.temp_kev());
            }
            //AGN photon fields photons are considered only if disk is present and compsw==2
            if(compsw==2 && Rin<Rout){
                InvCompton.bb_seed(Syncro.get_energy(),Ubb1,zone.delta*BLR.temp_kev());
                InvCompton.bb_seed(Syncro.get_energy(),Ubb2,zone.delta*Torus.temp_kev());
            }
            //Calculate the spectrum with whichever fields have been invoked		
            InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);

            if (infosw>=5){
                InvCompton.test();
            }
            //Sum counterjet if present, then save emission from the zone in com_en/com_lum		
            if (IsCounterjet==true){
                sum_counterjet(ncom,InvCompton.get_energy_obs(),InvCompton.get_nphot_obs(),com_en,com_lum);
            } else {
                for(int k=0;k<ncom;k++){
    	            com_en[k] = InvCompton.get_energy_obs()[k];
                    com_lum[k] = InvCompton.get_nphot_obs()[k];
	            }			
            }
            //Include zone's emission to the pre/post particle acceleration spectrum
            if(z<zdiss){
                sum_zones(ncom,ne,com_en,com_lum,tot_en,tot_com_pre);
            } else {
                sum_zones(ncom,ne,com_en,com_lum,tot_en,tot_com_post);
            }
        } else if (infosw>=5){
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
        plot_write(50,Disk.get_energy_obs(),Disk.get_nphot_obs(),"Output/Disk.dat",dist,redsh);
        //plot_write(ncom,Corona.get_energy_obs(),Corona.get_nphot_obs(),"Output/Corona.dat",dist,redsh);
        if(compsw==2){
            plot_write(40,Torus.get_energy_obs(),Torus.get_nphot_obs(),"Output/BB.dat",dist,redsh);
        } else {
            plot_write(40,BlackBody.get_energy_obs(),BlackBody.get_nphot_obs(),"Output/BB.dat",dist,redsh);
        }			
        plot_write(ne,tot_en,tot_lum,"Output/Total.dat",dist,redsh);
    }
    if (infosw >=3) {
        cout << "Observed 0.1-5 keV disk luminosity: "
             << integrate_lum(50,0.1*2.41e17,5.*2.41e17,Disk.get_energy_obs(),Disk.get_nphot_obs()) << endl;
        cout << "Observed 0.1-300 keV Inverse Compton luminosity: " 	
             << integrate_lum(ne,0.1*2.41e17,300.*2.41e17,tot_en,tot_com_pre) << endl; 
        cout << "Observed 0.1-300 keV total luminosity: " 
             << integrate_lum(ne,0.1*2.41e17,300.*2.41e17,tot_en,tot_lum) << endl; 
        cout << "Observed 3-7 GHz luminosity: " << integrate_lum(ne,3e9,7e9,tot_en,tot_lum) << endl;
        cout << "X-ray 10-100 keV photon index estimate: " 
             << photon_index(ne,10.*2.41e17,100.*2.41e17,tot_en,tot_lum) << endl;
        cout << "Radio 10-100 GHz spectral index estimate: " 
             << 1.+photon_index(ne,1e10,1e11,tot_en,tot_lum) << endl;
        double compactness = integrate_lum(ne,0.1*2.41e17,300.*2.41e17,tot_en,tot_com_pre)*sigtom
                             /(r0*emerg*cee);      
        cout << "Jet base compactness: " << compactness << endl;
        if (compactness >= 10.*(param[9]/511.)*exp(511./param[9])) {
            cout << "Possible runaway pair production in the jet base!" << endl; 
            cout << "Lower limit on allowed compactness: " << 10.*(param[9]/511.)*exp(511./param[9]) << endl;
            cout << "Note: this is for a slab, a cylinder allows higher l by a factor of a few" << std::endl;
        }
    }	

    gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
    gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);	
    gsl_spline_free(spline_speed), gsl_interp_accel_free(acc_speed);
    delete[] tot_en,delete[] tot_lum;
    delete[] tot_syn_pre,delete[] tot_syn_post;
    delete[] tot_com_pre,delete[] tot_com_post;
}
