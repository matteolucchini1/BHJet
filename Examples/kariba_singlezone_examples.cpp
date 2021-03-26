#include "kariba_examples.hh"

//This example shows the simplest possible blazar-type homogoeneous one zone jet model. The parameters are set to replicate Model
//2 from the EHT MWL paper on the 2017 M87 campaign, DOI: <insert>. The classes used here are Cyclosyn, Compton, and Powerlaw.

int main(){

    int nel = 100;			    //array size for particle distributions
    int nfreq = 200;            //array size for frequency arrays
    
    double B,R,n;               //plasma quantities of emitting region: bfield, radius, number density
    double gmin,gmax,p;         //details of particle distribution: minimum/maximum gamma factor, powerlaw slope
    double pmin;                //minimum momentum corresponding to gmin
    double delta,gamma,beta;    //plasma speed/beaming factor
    double theta;               //jet viewing angle
    double Pj;                  //total power carried in the jet 
    double Pe, Ue;              //power/energy density in electrons
    double Pb, Ub;              //power/energy density in bfields
    double Pp, Up;              //power/energy density in cold protons
    double equip;               //standard equipartition factor Ue/Ub
    double Mbh, Eddlum, Rg;     //black hole scale quantities
    double nus_min,nus_max;     //synchrotron frequency range
    double nuc_min,nuc_max;     //SSC frequency range 
    
    //splines for electron distribution. These are needed by the radiation codes
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);

    gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
    gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  
    
    //These calls remove the output of previous runs from the output files
	clean_file("Output/Singlezone_Syn.dat",1);
	clean_file("Output/Singlezone_SSC.dat",1);
	clean_file("Output/Singlezone_Particles.dat",1);
	
	//Set the model parameters
	Mbh = 6.5e9;
	Eddlum = 1.25e38*Mbh;				
	Rg = gconst*Mbh*msun/(cee*cee);
	B = 1.5e-3;
	R = 626.*Rg;
	n = 9.5e-3;
	gmin = 4.1e3;
	pmin = pow(pow(gmin,2.)-1.,1./2.)*emgm*cee;
	gmax = 6.4e7;
	p = 3.03;
	delta = 3.3;
	gamma = 3.;
	beta = 0.942809;
	theta = acos((delta*gamma-1.)/(beta*delta*gamma))*180./pi;
	
	nus_min = 1.e8;
	nus_max = 1.e21;
	nuc_min = 1.e17;
	nuc_max = 1.e28;
	
	//For the power-law class, the constructor only requires the array size for the particle distribution. You then need to set
	//the non-thermal slope of the distribution and specify the minimum momentum and maximum Lorenz factor of the distribution, in
	//no particular order. Then, set the normalisation from the number density and non-thermal slope, and set the distribution
	//array with set_ndens(). 
    Powerlaw Electrons(nel);
    Electrons.set_p(pmin,gmax); 
    Electrons.set_pspec(p);
    Electrons.set_norm(n);	
    Electrons.set_ndens();
    plot_write(nel,Electrons.get_p(),Electrons.get_gamma(),Electrons.get_pdens(),Electrons.get_gdens(),"Output/Singlezone_Particles.dat");
    
    gsl_spline_init(spline_eldis,Electrons.get_gamma(),Electrons.get_gdens(),nel);
    gsl_spline_init(spline_deriv,Electrons.get_gamma(),Electrons.get_gdens_diff(),nel); 
    
    //Set up the cyclo-synchrotron emission, by specifying the size of the arrays for frequency and flux. Then, specify the 
    //frequency range over which to calculate the spectrum, the magnetic field, the amount of beaming and viewing angble, and the
    //geometry of the emitting region> These can be done in no particular order. Finally, compute the spectrum by passing the
    //particle distribution Lorenz factor and interoplation objects. 
    Cyclosyn Syncro(nfreq);
    Syncro.set_frequency(nus_min,nus_max);    
    Syncro.set_bfield(B);
    Syncro.set_beaming(theta,beta,delta);
    Syncro.set_geometry("sphere",R);
    Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);
    plot_write(nfreq,Syncro.get_energy_obs(),Syncro.get_nphot_obs(),"Output/Singlezone_Syn.dat",1.,0.);
    
    //Set up the SSC emission, by specifying the size of both the output arrays, and the seed photon arrays. Then, set in no
    //particular order the amount of beaming, the frequency range where you want the emission to be computed, and the geometry of
    //the emitting region. After that, set the inverse Compton optical depth by passing the number density and temperature/average
    //Lorenz factor of the particles in units of keV. This will check whether multiple scatterings are needed and warn the user
    //if the parameters are outside of the parameters allowed by the physics code.
    //Finally, set up the SSC calculation by calling directly the co-moving photon energy and specify luminosity arrays from
    //Cyclosyn, and calculate the spectrum.
    Compton InvCompton(nfreq,nfreq);
    InvCompton.set_frequency(nuc_min,nuc_max);	    
    InvCompton.set_beaming(theta,beta,delta);
    InvCompton.set_geometry("sphere",R);	
    InvCompton.set_tau(n,Electrons.av_gamma()*511.);			
    InvCompton.cyclosyn_seed(Syncro.get_energy(),Syncro.get_nphot());
    InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    plot_write(nfreq,InvCompton.get_energy_obs(),InvCompton.get_nphot_obs(),"Output/Singlezone_SSC.dat",1.,0.);
    
    //Calculate physically relevant quantities like energy densities and powers, and output to terminal.
    Ue = Electrons.av_gamma()*n*emgm*pow(cee,2.);
    Ub = pow(B,2.)/(8.*pi);
    Up = n*pmgm*pow(cee,2.);
    equip = Ue/Ub;
    Pe = pi*beta*cee*pow(gamma*R,2.)*Ue;
    Pb = pi*beta*cee*pow(gamma*R,2.)*Ub;
    Pp = pi*beta*cee*pow(gamma*R,2.)*Up;
    Pj = Pe+Pb+Pp;
    
    cout << "Physical quantities:" << endl;
    cout << "Average electron Lorenz factor: " << Electrons.av_gamma() << endl;
    cout << "Equipartition Ue/UB: " << equip << endl;
    cout << "Electron power: " << Pe << endl;
	cout << "Magnetic power: " << Pb << endl;
	cout << "Proton power: " << Pp << endl;
	cout << "Total power: " << Pj << " erg s^-1, " << Pj/Eddlum << " Eddington" << endl;
	
	system("python3 Singlezone.py");
	
    gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);

	return 0;
}
