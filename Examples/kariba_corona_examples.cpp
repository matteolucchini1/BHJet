#include "kariba_examples.hh"

//This example shows how to set up thermal Comptonisation of accretion disk photons in Kariba by a spherical coronae, using the 
//ShSDisk, Thermal and Compton classes for a range of electrons temperatures and optical depths. The output is then compared to 
//CompPS in the python plotting script.

int main(){

    int nel = 100;			//array size for particle distributions
    int nfreq = 100;        //array size for frequency arrays
    
    //Input parameters:
	double Mbh;	
	double Eddlum;
	double Rg;
	double Rin,Rout;
	double Ldisk;

    //These calls remove the output of previous runs from the output files
	clean_file("Output/Disk.dat",1);
	clean_file("Output/IC_Tau260Te90.dat",1);
	clean_file("Output/IC_Tau076Te90.dat",1);
	clean_file("Output/IC_Tau019Te90.dat",1);
	clean_file("Output/IC_Tau260Te900.dat",1);
	clean_file("Output/IC_Tau076Te900.dat",1);
	clean_file("Output/IC_Tau019Te900.dat",1);
	
	//Disk parameters: black hole mass and corresponding Eddington luminosity and gravitational radii, innermost and outermost 
	//radii, bolometric luminosity in Eddington units.
	Mbh = 10.;
	Eddlum = 1.25e38*Mbh;				
	Rg = gconst*Mbh*msun/(cee*cee);
	Rin = 10.*Rg;
	Rout = 1e4*Rg;
	Ldisk = 1e-4;

    //Set up the optical depths, temepratures, number densities, and sizes for all our coronae. Note that because the seed photon
    //energy density in Kariba is calculated in simplified fashion at the center of the emititng region and neglecting GR, in this
    //case the radii of the corona should be thought of less as physical sizes and more as normalisations
	double Tau[3] = {2.6,0.76,0.19};
	double Te[2] = {90.,900.};
	double R[6] = {45.*Rg,75.*Rg,110.*Rg,4.*Rg,25.*Rg,130.*Rg};
	double ndens[6];
	
	//For each combination of temperature and optical depth, set up the number density. The values for the radii were chosen
	//exclusively to make the SED plot look nice.
	for (int i=0;i<6;i++) {
        ndens[i] = Tau[i%3]/(sigtom*R[i]);
	}

    //Set up the disk:note that the constructor does not require an array size (as many other classes, see later), because for a
    //disk it's always taken to be 50. Also note that you ALWAYS should specify the  inner and outer radii, THEN the luminosity. 
    //This is because the innermost radius is necessary to calculate the temperature profile, starting from the luminosity. The
    //mass of the BH is needed because the luminosity is expressed in Eddington units. After all this is done set the inclination,
    //calculate the spectrum, and save it to file
    ShSDisk Disk;
    Disk.set_mbh(Mbh);
    Disk.set_rin(Rin);
    Disk.set_rout(Rout);
    Disk.set_luminosity(Ldisk);		
    Disk.set_inclination(0.);
    Disk.disk_spectrum();
    plot_write(50,Disk.get_energy_obs(),Disk.get_nphot_obs(),"Output/Disk.dat",1.,0.);
    
    //Set up the first electron distribution. Call the constructor, which only requires the size of the arrays for the particles. 
    //As for the disk class, the arrays need to be called in the appropriate order. You can only set the momentum and Lorenz
    //factor arrays with set_p() after specifying a temperature, and you can only set up the number density arrays with 
    //set_ndens() after specifying the temperature and number density. 
    Thermal elec_Tau260Te90(nel);
    elec_Tau260Te90.set_temp_kev(Te[0]);
    elec_Tau260Te90.set_p();	
    elec_Tau260Te90.set_norm(ndens[0]);	
    elec_Tau260Te90.set_ndens();
    
    //Retrive the bounds of the electron distribution and set up the interpolation of the arrays. Both of these are needed as 
    //arguments for the Compton class.
    double gmin = elec_Tau260Te90.get_gamma()[0];
    double gmax = elec_Tau260Te90.get_gamma()[nel-1];
    
	//splines for electron distribution
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);
	gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
	gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  

    gsl_spline_init(spline_eldis,elec_Tau260Te90.get_gamma(),elec_Tau260Te90.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau260Te90.get_gamma(),elec_Tau260Te90.get_gdens_diff(),nel); 
    
    //Set up the inverse Compton calculation. As usual, you need to do some book-keeping before running the code.
    //The constructor this time requires both the size of the output arrays, and of the seed photon arrays, which for a disk is 
    //always set to 50. You then need to specify, in no particular order, the range over which you want the emission to be 
    //computed, what the geometry of the emitting region is and whehter it is beamed. 
    //After that, specify the temperature and number density of the radiating particles. This is required to calculate the 
    //radiative transfer corrections to the spectrum, if multiple scatterings are computed. Here, 20 scatterings are hard coded
    //for illustration purposes; the class constructor by default sets 20 scatterings, and if the optical depth is below 0.05 this
    //is set automatically to 1 by set_tau to save computational time.
    //Finally, specify the seed photon field - including only disk photons here, by passing an array of energies, the innermost
    //temperature of the disk, its inner and outer radii, the scale height of the disk (although it has only a minor effect on the
    //spectra), and the distance from the BH, over the BH vertical axis, of the emitting region.
    //After doing this, take the Lorenz factor interval and GSL interpolations of the particle distribution, compute the spectrum,
    //and save it to file
    Compton IC_Tau260Te90(nfreq,50);
    IC_Tau260Te90.set_frequency(1e15,1e22);
    IC_Tau260Te90.set_beaming(0.,0.,1.);
    IC_Tau260Te90.set_geometry("sphere",R[0]);
    IC_Tau260Te90.set_tau(ndens[0],Te[0]);
    IC_Tau260Te90.set_niter(20);		
    IC_Tau260Te90.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau260Te90.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(nfreq,IC_Tau260Te90.get_energy_obs(),IC_Tau260Te90.get_nphot_obs(),"Output/IC_Tau260Te90.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    //These parts are exactly the same as above. Rather than declare new objects we could also just update all the parameters of
    //of the initial Thermal and Compton objects with the setter methods, which can in principle improve performance by saving up
    //a bit of memory. Obviously that is not necessary for this example.      
    Thermal elec_Tau076Te90(nel);
    elec_Tau076Te90.set_temp_kev(Te[0]);
    elec_Tau076Te90.set_p();	
    elec_Tau076Te90.set_norm(ndens[1]);	
    elec_Tau076Te90.set_ndens();
    
    gmin = elec_Tau076Te90.get_gamma()[0];
    gmax = elec_Tau076Te90.get_gamma()[nel-1];
    
    gsl_spline_init(spline_eldis,elec_Tau076Te90.get_gamma(),elec_Tau076Te90.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau076Te90.get_gamma(),elec_Tau076Te90.get_gdens_diff(),nel); 

    Compton IC_Tau076Te90(nfreq,50);
    IC_Tau076Te90.set_frequency(1e15,1e22);
    IC_Tau076Te90.set_beaming(0.,0.,1.);
    IC_Tau076Te90.set_geometry("sphere",R[1]);	
    IC_Tau076Te90.set_tau(ndens[1],Te[0]);
    IC_Tau076Te90.set_niter(20);		
    IC_Tau076Te90.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau076Te90.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(nfreq,IC_Tau076Te90.get_energy_obs(),IC_Tau076Te90.get_nphot_obs(),"Output/IC_Tau076Te90.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau019Te90(nel);
    elec_Tau019Te90.set_temp_kev(Te[0]);
    elec_Tau019Te90.set_p();	
    elec_Tau019Te90.set_norm(ndens[2]);	
    elec_Tau019Te90.set_ndens();
    
    gmin = elec_Tau019Te90.get_gamma()[0];
    gmax = elec_Tau019Te90.get_gamma()[nel-1];
    
    gsl_spline_init(spline_eldis,elec_Tau019Te90.get_gamma(),elec_Tau019Te90.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau019Te90.get_gamma(),elec_Tau019Te90.get_gdens_diff(),nel); 

    Compton IC_Tau019Te90(nfreq,50);
    IC_Tau019Te90.set_frequency(1e15,1e22);
    IC_Tau019Te90.set_beaming(0.,0.,1.);
    IC_Tau019Te90.set_geometry("sphere",R[2]);
    IC_Tau019Te90.set_tau(ndens[2],Te[0]);
    IC_Tau019Te90.set_niter(20);		
    IC_Tau019Te90.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau019Te90.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(nfreq,IC_Tau019Te90.get_energy_obs(),IC_Tau019Te90.get_nphot_obs(),"Output/IC_Tau019Te90.dat",1.,0.);
    
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

    Compton IC_Tau260Te900(nfreq,50);
    IC_Tau260Te900.set_frequency(1e15,1e22);
    IC_Tau260Te900.set_beaming(0.,0.,1.);
    IC_Tau260Te900.set_geometry("sphere",R[3]);
    IC_Tau260Te900.set_tau(ndens[3],Te[1]);
    IC_Tau260Te900.set_niter(20);		
    IC_Tau260Te900.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau260Te900.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(nfreq,IC_Tau260Te900.get_energy_obs(),IC_Tau260Te900.get_nphot_obs(),"Output/IC_Tau260Te900.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau076Te900(nel);
    elec_Tau076Te900.set_temp_kev(Te[1]);
    elec_Tau076Te900.set_p();	
    elec_Tau076Te900.set_norm(ndens[4]);	
    elec_Tau076Te900.set_ndens();
    
    gmin = elec_Tau076Te900.get_gamma()[0];
    gmax = elec_Tau076Te900.get_gamma()[nel-1];
    
    gsl_spline_init(spline_eldis,elec_Tau076Te900.get_gamma(),elec_Tau076Te900.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau076Te900.get_gamma(),elec_Tau076Te900.get_gdens_diff(),nel); 

    Compton IC_Tau076Te900(nfreq,50);
    IC_Tau076Te900.set_frequency(1e15,1e22);
    IC_Tau076Te900.set_beaming(0.,0.,1.);
    IC_Tau076Te900.set_geometry("sphere",R[4]);	
    IC_Tau076Te900.set_tau(ndens[4],Te[1]);
    IC_Tau076Te900.set_niter(20);		
    IC_Tau076Te900.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau076Te900.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(nfreq,IC_Tau076Te900.get_energy_obs(),IC_Tau076Te900.get_nphot_obs(),"Output/IC_Tau076Te900.dat",1.,0.);
    
    //---------------------------------------------------------------------------------------------------------------
    Thermal elec_Tau019Te900(nel);
    elec_Tau019Te900.set_temp_kev(Te[1]);
    elec_Tau019Te900.set_p();	
    elec_Tau019Te900.set_norm(ndens[5]);	
    elec_Tau019Te900.set_ndens();
    
    gmin = elec_Tau019Te900.get_gamma()[0];
    gmax = elec_Tau019Te900.get_gamma()[nel-1];
    
    gsl_spline_init(spline_eldis,elec_Tau019Te900.get_gamma(),elec_Tau019Te900.get_gdens(),nel);
    gsl_spline_init(spline_deriv,elec_Tau019Te900.get_gamma(),elec_Tau019Te900.get_gdens_diff(),nel); 

    Compton IC_Tau019Te900(nfreq,50);
    IC_Tau019Te900.set_frequency(1e15,1e22);
    IC_Tau019Te900.set_beaming(0.,0.,1.);
    IC_Tau019Te900.set_geometry("sphere",R[5]);	
    IC_Tau019Te900.set_tau(ndens[5],Te[1]);
    IC_Tau019Te900.set_niter(20);		
    IC_Tau019Te900.shsdisk_seed(Disk.get_energy(),Disk.tin(),Rin,Rout,Disk.hdisk(),0.);
    IC_Tau019Te900.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    
    plot_write(nfreq,IC_Tau019Te900.get_energy_obs(),IC_Tau019Te900.get_nphot_obs(),"Output/IC_Tau019Te900.dat",1.,0.);       
    
	system("python3 Coronae.py");		

	gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
	gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);

	return 0;
}
