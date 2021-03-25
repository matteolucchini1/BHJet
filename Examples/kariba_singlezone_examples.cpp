#include "kariba_examples.hh"

int main(){

    int nel = 100;			//array size for particle distributions
    
    //Input parameters read from file:
	clean_file("Output/Singlezone_Syn",1);
	clean_file("Output/Singlezone_SSC",1);
	clean_file("Output/Singlezone_Particles.dat",1);
	
	Mbh = 1e9;
	Eddlum = 1.25e38*Mbh;				
	Rg = gconst*Mbh*msun/(cee*cee);
	Rin = 10.*Rg;
	Rout = 1e4*Rg;
	Ldisk = 1e-4;

    //Note: the radii are just re-normalisation factors	
	double Tau[3] = {2.6,0.76,0.19};
	double Te[2] = {90.,900.};
	double R[6] = {45.*Rg,75.*Rg,110.*Rg,4.*Rg,25.*Rg,130.*Rg};
	double ndens[6];
	
	for (int i=0;i<6;i++) {
        ndens[i] = Tau[i%3]/(sigtom*R[i]);
	}

    ShSDisk Disk;
    Disk.set_mbh(Mbh);
    Disk.set_rin(Rin);
    Disk.set_rout(Rout);
    Disk.set_luminosity(Ldisk);		
    Disk.set_inclination(0.);
    Disk.disk_spectrum();
    plot_write(50,Disk.get_energy_obs(),Disk.get_nphot_obs(),"Output/Disk.dat",1.,0.);
    
    Thermal elec_Tau260Te90(nel);
    elec_Tau260Te90.set_temp_kev(Te[0]);
    elec_Tau260Te90.set_p();	
    elec_Tau260Te90.set_norm(ndens[0]);	
    elec_Tau260Te90.set_ndens();
    
    double gmin = elec_Tau260Te90.get_gamma()[0];
    double gmax = elec_Tau260Te90.get_gamma()[nel-1];
    
	//splines for electron distribution
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);
	gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
	gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  

    gsl_spline_init(spline_eldis,elec_Tau260Te90.get_gamma(),elec_Tau260Te90.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau260Te90.get_gamma(),elec_Tau260Te90.get_gdens_diff(),nel); 

    Compton IC_Tau260Te90(100,50);
    IC_Tau260Te90.set_frequency(1e15,1e22);
    IC_Tau260Te90.set_beaming(0.,0.,1.);
    IC_Tau260Te90.set_geometry("sphere",R[0]);
    IC_Tau260Te90.set_counterjet(false);	
    IC_Tau260Te90.set_tau(ndens[0],Te[0]);
    IC_Tau260Te90.set_niter(20);		
    IC_Tau260Te90.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau260Te90.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(100,IC_Tau260Te90.get_energy_obs(),IC_Tau260Te90.get_nphot_obs(),"Output/IC_Tau260Te90.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau076Te90(nel);
    elec_Tau076Te90.set_temp_kev(Te[0]);
    elec_Tau076Te90.set_p();	
    elec_Tau076Te90.set_norm(ndens[1]);	
    elec_Tau076Te90.set_ndens();
    
    gmin = elec_Tau076Te90.get_gamma()[0];
    gmax = elec_Tau076Te90.get_gamma()[nel-1];
    
	//splines for electron distribution
    gsl_spline_init(spline_eldis,elec_Tau076Te90.get_gamma(),elec_Tau076Te90.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau076Te90.get_gamma(),elec_Tau076Te90.get_gdens_diff(),nel); 

    Compton IC_Tau076Te90(100,50);
    IC_Tau076Te90.set_frequency(1e15,1e22);
    IC_Tau076Te90.set_beaming(0.,0.,1.);
    IC_Tau076Te90.set_geometry("sphere",R[1]);
    IC_Tau076Te90.set_counterjet(false);	
    IC_Tau076Te90.set_tau(ndens[1],Te[0]);
    IC_Tau076Te90.set_niter(20);		
    IC_Tau076Te90.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau076Te90.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(100,IC_Tau076Te90.get_energy_obs(),IC_Tau076Te90.get_nphot_obs(),"Output/IC_Tau076Te90.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau019Te90(nel);
    elec_Tau019Te90.set_temp_kev(Te[0]);
    elec_Tau019Te90.set_p();	
    elec_Tau019Te90.set_norm(ndens[2]);	
    elec_Tau019Te90.set_ndens();
    
    gmin = elec_Tau019Te90.get_gamma()[0];
    gmax = elec_Tau019Te90.get_gamma()[nel-1];
    
	//splines for electron distribution
    gsl_spline_init(spline_eldis,elec_Tau019Te90.get_gamma(),elec_Tau019Te90.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau019Te90.get_gamma(),elec_Tau019Te90.get_gdens_diff(),nel); 

    Compton IC_Tau019Te90(100,50);
    IC_Tau019Te90.set_frequency(1e15,1e22);
    IC_Tau019Te90.set_beaming(0.,0.,1.);
    IC_Tau019Te90.set_geometry("sphere",R[2]);
    IC_Tau019Te90.set_counterjet(false);	
    IC_Tau019Te90.set_tau(ndens[2],Te[0]);
    IC_Tau019Te90.set_niter(20);		
    IC_Tau019Te90.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau019Te90.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(100,IC_Tau019Te90.get_energy_obs(),IC_Tau019Te90.get_nphot_obs(),"Output/IC_Tau019Te90.dat",1.,0.);
    
    //----------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------
    
    Thermal elec_Tau260Te900(nel);
    elec_Tau260Te900.set_temp_kev(Te[1]);
    elec_Tau260Te900.set_p();	
    elec_Tau260Te900.set_norm(ndens[3]);	
    elec_Tau260Te900.set_ndens();
    
    gmin = elec_Tau260Te900.get_gamma()[0];
    gmax = elec_Tau260Te900.get_gamma()[nel-1];
    
    gsl_spline_init(spline_eldis,elec_Tau260Te900.get_gamma(),elec_Tau260Te900.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau260Te900.get_gamma(),elec_Tau260Te900.get_gdens_diff(),nel); 

    Compton IC_Tau260Te900(100,50);
    IC_Tau260Te900.set_frequency(1e15,1e22);
    IC_Tau260Te900.set_beaming(0.,0.,1.);
    IC_Tau260Te900.set_geometry("sphere",R[3]);
    IC_Tau260Te900.set_counterjet(false);	
    IC_Tau260Te900.set_tau(ndens[3],Te[1]);
    IC_Tau260Te900.set_niter(20);		
    IC_Tau260Te900.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau260Te900.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(100,IC_Tau260Te900.get_energy_obs(),IC_Tau260Te900.get_nphot_obs(),"Output/IC_Tau260Te900.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau076Te900(nel);
    elec_Tau076Te900.set_temp_kev(Te[1]);
    elec_Tau076Te900.set_p();	
    elec_Tau076Te900.set_norm(ndens[4]);	
    elec_Tau076Te900.set_ndens();
    
    gmin = elec_Tau076Te900.get_gamma()[0];
    gmax = elec_Tau076Te900.get_gamma()[nel-1];
    
	//splines for electron distribution
    gsl_spline_init(spline_eldis,elec_Tau076Te900.get_gamma(),elec_Tau076Te900.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau076Te900.get_gamma(),elec_Tau076Te900.get_gdens_diff(),nel); 

    Compton IC_Tau076Te900(100,50);
    IC_Tau076Te900.set_frequency(1e15,1e22);
    IC_Tau076Te900.set_beaming(0.,0.,1.);
    IC_Tau076Te900.set_geometry("sphere",R[4]);
    IC_Tau076Te900.set_counterjet(false);	
    IC_Tau076Te900.set_tau(ndens[4],Te[1]);
    IC_Tau076Te900.set_niter(20);		
    IC_Tau076Te900.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau076Te900.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(100,IC_Tau076Te900.get_energy_obs(),IC_Tau076Te900.get_nphot_obs(),"Output/IC_Tau076Te900.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau019Te900(nel);
    elec_Tau019Te900.set_temp_kev(Te[1]);
    elec_Tau019Te900.set_p();	
    elec_Tau019Te900.set_norm(ndens[5]);	
    elec_Tau019Te900.set_ndens();
    
    gmin = elec_Tau019Te900.get_gamma()[0];
    gmax = elec_Tau019Te900.get_gamma()[nel-1];
    
	//splines for electron distribution
    gsl_spline_init(spline_eldis,elec_Tau019Te900.get_gamma(),elec_Tau019Te900.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau019Te900.get_gamma(),elec_Tau019Te900.get_gdens_diff(),nel); 

    Compton IC_Tau019Te900(100,50);
    IC_Tau019Te900.set_frequency(1e15,1e22);
    IC_Tau019Te900.set_beaming(0.,0.,1.);
    IC_Tau019Te900.set_geometry("sphere",R[5]);
    IC_Tau019Te900.set_counterjet(false);	
    IC_Tau019Te900.set_tau(ndens[5],Te[1]);
    IC_Tau019Te900.set_niter(20);		
    IC_Tau019Te900.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau019Te900.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(100,IC_Tau019Te900.get_energy_obs(),IC_Tau019Te900.get_nphot_obs(),"Output/IC_Tau019Te900.dat",1.,0.);
        
    
	system("python3 Coronae.py");		

	gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
	gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);

	return 0;
}
