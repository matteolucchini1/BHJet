#include <fstream>
#include <stdarg.h>

#include "bhjet.hh"

using namespace std;

void jetmain(double *ear,int ne,double *param,double *photeng,double *photspec) {

    //STEP 1: VARIABLE/OBJECT DEFINITIONS
    //----------------------------------------------------------------------------------------------

    bool IsShock = false;						//flag to set shock heating

    int nz = 100;								//total number of zones
    int nel = 70;
    int syn_res = 10;							//number of bins per decade in synch frequency;
    int com_res = 6;							//number of bins per decade in compton frequency;
    int nsyn,ncom;								//number of bins in synch/compton frequency;
    int npsw = 1;								//switch to define number of protons calculations in agnjet
											    //0: no protons
											    //1: Up = Ue+Ub
											    //2: ne = np	
											    //NOTE pair only jet not implemented self-consistently yet!

    double Mbh;									//black hole mass in solar masses
    double Eddlum;								//black hole Eddington luminosity
    double Rg;									//black hole  gravitational radius in cm	
    double theta;								//viewing angle	
    double dist;								//distance in kpc	
    double redsh;								//source redshift
    double jetrat;								//injected power in Eddington units	
    double zmin;								//jet launching point
    double r_0;									//initial jet radius in rg
    double h;									//jet nozzle aspect ratio
    double z_acc;								//end of bulk magnetic acceleration region in rg
    double z_diss;								//dissipation/nonthermal particle injection region in rg	
    double z_max;								//distance from bh up to which calculation continues
    double t_e;									//temperature in kev, converted to erg	
    double f_nth;							    //percentage of nonthermal particles at the dissipation region
    double f_pl;								//parameter to change plfrac over distance
    double pspec;								//slope of nonthermal distribution
    double f_heat;								//shock heating paramter
    double f_beta;								//effective expansion velocity used to set adiabatic cooling
    double f_sc;								//particle acceleration timescale parameter		
    double p_beta;								//plasma beta in the jet base
    double sig_acc;								//final sigma when using magnetic acceleration
    double l_disk;								//disk luminosity in Eddington units
    double r_in;								//disk inner radius in rg
    double r_out;								//disk outer radius in rg
    double compar1;								//external inverse Compton parameters; different meanings
    double compar2;								//depending on the value of compsw
    double compar3;	
    double compsw;								//switch to activate different external Compton fields
    double velsw;								//velocity profile parameter
    int infosw;									//switch to print info
    int EBLsw;									//switch to activate EBL attenuation

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
    
    //External photon object declarations
    ShSDisk Disk;
    BBody BLR;
    BBody Torus;
    BBody BlackBody;

    //splines for jet acceleration
    gsl_interp_accel *acc_speed = gsl_interp_accel_alloc();
    gsl_spline *spline_speed = gsl_spline_alloc(gsl_interp_steffen,54);

    //splines for electron distribution
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);

    gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
    gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  

    //STEP 2: PARAMETER/FILE INITIALIZATION
    Mbh = param[0];
    Eddlum = 1.25e38*Mbh;				
    Rg = gconst*Mbh*msun/(cee*cee);
    theta = param[1];		
    dist = param[2]*kpc;	
    redsh = param[3];
    jetrat = param[4]*Eddlum;
    r_0 = param[5]*Rg;
    z_diss = param[6]*Rg;
    z_acc = param[7]*Rg;	
    z_max = param[8]*Rg;
    t_e = param[9];
    f_nth = param[10];
    f_pl = param[11];
    pspec = param[12];
    f_heat = param[13];
    f_beta = param[14];
    f_sc = param[15];
    p_beta = param[16];
    sig_acc = param[17];	
    l_disk = param[18];
    r_in = param[19]*Rg;
    r_out = param[20]*Rg;
    compar1 = param[21];
    compar2 = param[22];
    compar3 = param[23];
    compsw = param[24];
    velsw = param[25];
    infosw = param[26];
    EBLsw = param[27];
    zmin = 2.*Rg;
    
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

    if (infosw>=1) {			
        clean_file("Output/Presyn.dat",2);
        clean_file("Output/Postsyn.dat",2);
        clean_file("Output/Precom.dat",2);
        clean_file("Output/Postcom.dat",2);
        clean_file("Output/Disk.dat",2);
        clean_file("Output/BB.dat",2);
        clean_file("Output/Total.dat",2);	
    }
    if (infosw>=2) {	
        clean_file("Output/Numdens.dat",4);
        clean_file("Output/Cyclosyn_zones.dat",2);
        clean_file("Output/Compton_zones.dat",2);			
    }
    if (infosw>=3) {
        clean_file("Output/Spectral_properties.dat",7);    
    }
    if (infosw>=5) {
        clean_file("Output/Profiles.dat",6);
    }
    //STEP 3: DISK/EXTERNAL PHOTON CALCULATIONS
    //The disk is disabled by setting r_in<r_out; its contribution is summed to the total only if there are no 
    //AGN photon fields that reprocess part of the luminosity, otherwise it is done later
    if(r_in<r_out){
    	Disk.set_mbh(Mbh);
	    Disk.set_rin(r_in);
	    Disk.set_rout(r_out);
	    Disk.set_luminosity(abs(l_disk));		
	    Disk.set_inclination(theta);
	    Disk.disk_spectrum();
	    if (compsw != 2 && l_disk > 0) {
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
    else if (compsw==2 && r_in<r_out){
        agn_photons_init(Disk.total_luminosity(),compar1,compar2,agn_com);

        BLR.set_temp_kev(agn_com.tblr);
        BLR.set_lum(agn_com.lblr);
        BLR.bb_spectrum();			

        Torus.set_temp_k(agn_com.tdt);
        Torus.set_lum(agn_com.ldt);	
        Torus.bb_spectrum();
        
        Disk.cover_disk(compar1+compar2);		
        if(infosw>=3){
            cout << "BLR radius in Rg: " << agn_com.rblr/Rg << " and in cm: " << agn_com.rblr<< endl;
            cout << "DT radius in Rg: " << agn_com.rdt/Rg  << " and in cm: " << agn_com.rdt<< endl;
        }
        sum_ext(40,ne,Torus.get_energy_obs(),Torus.get_nphot_obs(),tot_en,tot_lum);
        if (l_disk > 0) {
            sum_ext(50,ne,Disk.get_energy_obs(),Disk.get_nphot_obs(),tot_en,tot_lum);    
        }          		
    }	
	
    //STEP 4: JET BASE EQUIPARTITION CALCULATIONS AND SETUP
    //Dummy particle distribution, needed for average lorentz factor in equipartition function
    //The number density is just set to unity, the normalisation is not needed to calculate the average Lorenz factor
    //of the thermal distribution anyway
    Thermal dummy_elec(nel);
    dummy_elec.set_temp_kev(t_e);
    dummy_elec.set_p();	
    dummy_elec.set_norm(1.);	
    dummy_elec.set_ndens();

    grid.nz = nz;
    grid.cut = 0;
    grid.zcut = 1.e3*Rg;

    jet_dyn.min = zmin;
    jet_dyn.max = z_max;
    jet_dyn.h0 = 2.*r_0+zmin;
    jet_dyn.r0 = r_0;
    jet_dyn.acc = z_acc;
    jet_dyn.beta0 = sqrt(4./3.*(4./3.-1.)/(4./3.+1.));	//set initial jet speed for relativistic fluid, g=4/3 
    jet_dyn.gam0 = 1./sqrt(1.-(pow(jet_dyn.beta0,2.)));	//set corresponding lorentz factor
    jet_dyn.gamf = velsw;
    jet_dyn.Rg = Rg;  

    nozzle_ener.pbeta = p_beta;
    nozzle_ener.Nj = jetrat;
    nozzle_ener.sig_acc = sig_acc; 
    nozzle_ener.av_gamma = dummy_elec.av_gamma();
    //set up jet velocity profile depending on choice of adiabatic,isothermal,magnetically dominated jet   
    //note: the adiabatic jet only runs correctly if the final temperature is above ~1kev, which means the
    //initial temperature has to be ~10^4 kev to avoid numerical issues
    if(velsw==0){   
        velprof_ad(spline_speed);
        equipartition(npsw,jet_dyn,nozzle_ener);	
    } else if(velsw==1){
      	velprof_iso(spline_speed);
      	equipartition(npsw,jet_dyn,nozzle_ener);
    } else {
        velprof_mag(jet_dyn,spline_speed);
        equipartition(jetrat,jet_dyn,nozzle_ener);
    }	

    //check that the pair content is not negative, and also if running bljet that it's not too high
    if(nozzle_ener.eta<1){
        cout << "Unphysical pair content: " << nozzle_ener.eta << " pairs per proton. Check the value of " <<
                "plasma beta!" << endl;
    } else if (velsw>1 && dummy_elec.av_gamma()*nozzle_ener.eta >= 3e2){
           cout << "Pair content or temperature too high for  for bljet! " << endl;
           cout << "Pair content: " << nozzle_ener.eta << " pairs per proton" << endl;
           cout << "Average lepton Lorenz factor: " << dummy_elec.av_gamma() << endl;
           cout << "Check the value of Te and/or plasma beta!" << endl;
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
    for(int i=0;i<nz;i++){
        //calculate dynamics/energetics in each zone
        jetgrid(i,grid,jet_dyn,zone.r,zone.delz,z);	
        if(velsw==0){
            adjetpars(z,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);	 
        } else if (velsw==1){
            isojetpars(z,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);
        } else {
            bljetpars(z,f_beta,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);
        }
        zone.delta = 1./(zone.gamma*(1.-zone.beta*cos(theta*pi/180.)));

        //This is to avoid crashes due to low (sub 1 kev) particle temperatures
        if (z < z_diss) {
            zone.eltemp = max(tshift*t_e,1.);    
        } else {
            zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
        }
                
        //This is to evolve the fraction of non thermal particles along the jet, and change the distribution 
        //of non thermal particles appropriately
        if (z < z_diss) {
            zone.nth_frac = 0.;
        } else {
            zone.nth_frac = f_nth*pow(log10(z_diss)/log10(z),f_pl);           
        }

        //Include the disk for radiative cooling if it's on;
        //if compsw==1 the boosting is the same in all zones, if compsw==2 we boost either blr,torus or both
        //depending on the distance from the BH of the jet region/seed photon production
        if (r_in < r_out) {
            double Rdisk = pow(r_in,2.)+pow(z,2.);
            double delta_disk, theta_disk;
            theta_disk = pi-atan(r_in/z);
            delta_disk = 1./(zone.gamma-zone.beta*cos(theta_disk));
            Urad = pow(delta_disk,2.)*l_disk*Eddlum/(4.*pi*Rdisk*cee);
        } else {
            Urad = 0.;
        }
        
        if(compsw==1){			
            Urad = Urad + pow(zone.delta,2.)*Ubb1;
        } else if (compsw==2 && r_in<r_out){
            zone_agn_phfields(z,zone,Ubb1,Ubb2,agn_com);
            Urad = Urad + agn_com.urad_total;
        }          

        //calculate particle distribution in each zone
        if(zone.nth_frac == 0.){		
            Thermal th_lep(nel);
            th_lep.set_temp_kev(zone.eltemp);
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
                t_e = f_heat*t_e;
                IsShock = true;
                zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
            }					
            Mixed acc_lep(nel);
            acc_lep.set_temp_kev(zone.eltemp);
            acc_lep.set_pspec(pspec);
            acc_lep.set_plfrac(zone.nth_frac);
            
            //if f_sc < 10 it's the acceleration efficiency, else it's the desired maximum lorentz factor
            if (f_sc<10.){
                acc_lep.set_p(Urad,zone.bfield,f_beta,zone.r,f_sc);
            } else{
                acc_lep.set_p(f_sc);
            }	

            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,f_beta);
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
                t_e = f_heat*t_e;
                zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
                IsShock = true;
            }
            Thermal dummy_elec(nel);
            dummy_elec.set_temp_kev(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pbrk = dummy_elec.av_p();
            
            Bknpower acc_lep(nel);
            acc_lep.set_pspec1(-2.);
            acc_lep.set_pspec2(pspec);
                        
            if (f_sc<10.){
                acc_lep.set_p(0.1*pbrk,pbrk,Urad,zone.bfield,f_beta,zone.r,f_sc);
            } else{
                acc_lep.set_p(0.1*pbrk,pbrk,f_sc);
            }	
            
            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,f_beta);
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
                t_e = f_heat*t_e;
                zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
                IsShock = true;
            }
            Thermal dummy_elec(nel);
            dummy_elec.set_temp_kev(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pmin = dummy_elec.av_p();
            
            Powerlaw acc_lep(nel);
            acc_lep.set_pspec(pspec);
                        
            if (f_sc<10.){
                acc_lep.set_p(pmin,Urad,zone.bfield,f_beta,zone.r,f_sc);
            } else{
                acc_lep.set_p(pmin,f_sc);
            }	
            
            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,f_beta);
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
            cout << "tshift: " << tshift << ", Temperature in keV: " << zone.eltemp << endl;
            cout << "Grid; R: " << zone.r/Rg << ", delz: " << zone.delz/Rg << ", z: " << z/Rg << ", z+delz: " 
                 << (zone.delz+z)/Rg << endl;
            cout << "Equipartition check; Sigma: " << 2.*Ub/Up << " Ue/Ub: " << Ue/Ub << endl;
            
            std::ofstream file;
            file.open("Output/Profiles.dat",std::ios::app);	
            file << z/Rg << " " << zone.r/Rg << " " << zone.bfield << " " << zone.lepdens << " " << zone.gamma << " " <<
                    zone.eltemp << " " << std::endl;
            file.close();	
        }
        //calculate emission of each zone		
        //note: the syn_en array is used for the seed photon fields in the IC part, so it needs to include
        //both the black body and disk part. This is why the maximum frequency is taken as the maximum of the
        //two scale frequencies.
        syn_min = 0.1*pow(gmin,2.)*charg*zone.bfield/(2.*pi*emgm*cee);
        if(r_in<r_out){
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
        
        if(infosw>1) {
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
        Syncro.set_counterjet(true);
        Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);
        sum_counterjet(nsyn,Syncro.get_energy_obs(),Syncro.get_nphot_obs(),syn_en,syn_lum);	        
        if (infosw>=4){
            Syncro.test();
        }  
        //Include zone's emission to the pre/post particle acceleration spectrum
        if(z<z_diss){
            sum_zones(nsyn,ne,syn_en,syn_lum,tot_en,tot_syn_pre);
        } else {
            sum_zones(nsyn,ne,syn_en,syn_lum,tot_en,tot_syn_post);
        }		
        	
        //calculate inverse Compton spectrum, if it's expected to be bright enough	
        if (Compton_check(IsShock,i,Mbh,jetrat,Urad,velsw,zone) == true){
        //if(z>z_max){
            //Set up the calculation by reading in/calculating beaming,volume,counterjet presence,tau
            InvCompton.set_beaming(theta,zone.beta,zone.delta);
            InvCompton.set_geometry("cylinder",zone.r,zone.delz);
            InvCompton.set_counterjet(true);	
            InvCompton.set_tau(zone.lepdens,zone.eltemp);
            //Multiple scatters only if ypar and tau are large enough
            if(InvCompton.get_ypar() > 1.e-2 && InvCompton.get_tau() > 5.e-2){
                InvCompton.set_niter(15);		
            }
            //Cyclosynchrotron photons are always considered in the scattering						
            InvCompton.cyclosyn_seed(Syncro.get_energy(),Syncro.get_nphot());
            
            //Disk photons are included only if the disk is present
            if(r_in<r_out){
                InvCompton.shsdisk_seed(Syncro.get_energy(),Disk.tin(),r_in,r_out,Disk.hdisk(),z+zone.delz/2.);
            }
            //Black body photons included only if compsw==1
            if(compsw==1){						
                InvCompton.bb_seed_k(Syncro.get_energy(),Ubb1,zone.delta*BlackBody.temp_k());
            }
            //AGN photon fields photons are considered only if disk is present and compsw==2
            if(compsw==2 && r_in<r_out){
                InvCompton.bb_seed_k(Syncro.get_energy(),Ubb1,zone.delta*BLR.temp_k());
                InvCompton.bb_seed_k(Syncro.get_energy(),Ubb2,zone.delta*Torus.temp_k());
            }
            //Calculate the spectrum with whichever fields have been invoked		
            InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
            sum_counterjet(ncom,InvCompton.get_energy_obs(),InvCompton.get_nphot_obs(),com_en,com_lum);             
            if (infosw>=4){
                InvCompton.test();
            }
            
            //Include zone's emission to the pre/post particle acceleration spectrum
            if(z<z_diss){
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
    for(int k=0;k<ne;k++){
        tot_lum[k] = (tot_lum[k]+tot_syn_pre[k]+tot_syn_post[k]+tot_com_pre[k]+tot_com_post[k]);
        photeng[k] = log10(tot_en[k]/herg);
    }
    
    // Apply EBL attenuation factor for extragalactic sources
    if(redsh > 0. && EBLsw == 1){
        ebl_atten_gil(ne,tot_en,tot_lum,redsh); //correction for total luminosity
    	ebl_atten_gil(ne,tot_en,tot_com_post,redsh); //correction for post Compton luminosity
    }
    output_spectrum(ne,tot_en,tot_lum,photspec,redsh,dist);
	
	//Output to files and print information on terminal if user requires it
    if(infosw>=1){
        plot_write(ne,tot_en,tot_syn_pre,"Output/Presyn.dat",dist,redsh);
        plot_write(ne,tot_en,tot_syn_post,"Output/Postsyn.dat",dist,redsh);
        plot_write(ne,tot_en,tot_com_pre,"Output/Precom.dat",dist,redsh);
        plot_write(ne,tot_en,tot_com_post,"Output/Postcom.dat",dist,redsh);
        plot_write(50,Disk.get_energy_obs(),Disk.get_nphot_obs(),"Output/Disk.dat",dist,redsh);
        if(compsw==2){
            plot_write(40,Torus.get_energy_obs(),Torus.get_nphot_obs(),"Output/BB.dat",dist,redsh);
        } else {
            plot_write(40,BlackBody.get_energy_obs(),BlackBody.get_nphot_obs(),"Output/BB.dat",dist,redsh);
        }			
        plot_write(ne,tot_en,tot_lum,"Output/Total.dat",dist,redsh);
    }
    if (infosw >=3) {
        double disk_lum,IC_lum,Xray_lum,Radio_lum,Xray_index,Radio_index,compactness;
        disk_lum = integrate_lum(50,0.3*2.41e17,5.*2.41e17,Disk.get_energy_obs(),Disk.get_nphot_obs());
        IC_lum = integrate_lum(ne,0.3*2.41e17,300.*2.41e17,tot_en,tot_com_pre);
        Xray_lum = integrate_lum(ne,1.*2.41e17,10.*2.41e17,tot_en,tot_lum);
        Radio_lum = integrate_lum(ne,4e9,6e9,tot_en,tot_lum);
        Xray_index = photon_index(ne,10.*2.41e17,100.*2.41e17,tot_en,tot_lum);
        Radio_index = 1.+photon_index(ne,1e10,1e11,tot_en,tot_lum);
        compactness = integrate_lum(ne,0.1*2.41e17,300.*2.41e17,tot_en,tot_com_pre)*sigtom/(r_0*emerg*cee);   
        cout << "Observed 0.3-5 keV disk luminosity: " << disk_lum << endl;
        cout << "Observed 0.3-300 keV Inverse Compton luminosity: " << IC_lum << endl; 
        cout << "Observed 1-10 keV total luminosity: " << Xray_lum << endl; 
        cout << "Observed 4-6 GHz luminosity: " << Radio_lum << endl;
        cout << "X-ray 10-100 keV photon index estimate: " << Xray_index << endl;
        cout << "Radio 10-100 GHz spectral index estimate: " << Radio_index << endl;                   
        cout << "Jet base compactness: " << compactness << endl << endl;
        std::ofstream file;
        file.open("Output/Spectral_properties.dat",std::ios::app);	
        file << disk_lum << " " << IC_lum << " " << Xray_lum << " " << Radio_lum << " " << Xray_index << " " << Radio_index <<
            " " << compactness << std::endl;
        file.close();
        if (compactness >= 10.*(param[9]/511.)*exp(511./param[9])) {
            cout << "Possible pair production in the jet base!" << endl; 
            cout << "Lower limit on allowed compactness: " << 10.*(param[9]/511.)*exp(511./param[9]) << endl;
            cout << "Note: this is for a slab, a cylinder allows higher l by a factor of ~10" << std::endl;}
    }	

    gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
    gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);	
    gsl_spline_free(spline_speed), gsl_interp_accel_free(acc_speed);
    delete[] tot_en,delete[] tot_lum;
    delete[] tot_syn_pre,delete[] tot_syn_post;
    delete[] tot_com_pre,delete[] tot_com_post;
}
