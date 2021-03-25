#include "kariba_examples.hh"

int main(){

    int nel = 100;			    //array size for particle distributions
    int nfreq = 200;            //array size for frequency arrays
    
    double B,R,n;               //plasma quantities of emitting region: bfield, radius, number density
    double gmin,gmax,p;         //details of particle distribution: minimum/maximum gamma factor, powerlaw slope
    double pmin,avgammasq;      //derived quantities from particle distribution
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
    
    //splines for electron distribution
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);

    gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
    gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  
    
    //Clean output files
	clean_file("Output/Singlezone_Syn.dat",1);
	clean_file("Output/Singlezone_SSC.dat",1);
	clean_file("Output/Singlezone_Particles.dat",1);
	
	//These parameters are set to replicate Model 2 from the EHT MWL paper on the 2017 M87 campaign, DOT: <insert>
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
	
    Powerlaw Electrons(nel);
    Electrons.set_pspec(p);
    Electrons.set_p(pmin,gmax); 
    Electrons.set_norm(n);	
    Electrons.set_ndens();
    avgammasq = pow(Electrons.av_gamma(),2.);
    plot_write(nel,Electrons.get_p(),Electrons.get_gamma(),Electrons.get_pdens(),Electrons.get_gdens(),"Output/Singlezone_Particles.dat");
    
    gsl_spline_init(spline_eldis,Electrons.get_gamma(),Electrons.get_gdens(),nel);
    gsl_spline_init(spline_deriv,Electrons.get_gamma(),Electrons.get_gdens_diff(),nel); 
    
    Cyclosyn Syncro(nfreq);
    Syncro.set_frequency(nus_min,nus_max);    
    Syncro.set_bfield(B);
    Syncro.set_beaming(theta,beta,delta);
    Syncro.set_geometry("sphere",R);
    Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);
    plot_write(nfreq,Syncro.get_energy_obs(),Syncro.get_nphot_obs(),"Output/Singlezone_Syn.dat",1.,0.);
    
    Compton InvCompton(nfreq,nfreq);
    InvCompton.set_frequency(nuc_min,nuc_max);	    
    InvCompton.set_beaming(theta,beta,delta);
    InvCompton.set_geometry("sphere",R);	
    InvCompton.set_tau(n,Electrons.av_gamma());
    InvCompton.set_niter(1);					
    InvCompton.cyclosyn_seed(Syncro.get_energy(),Syncro.get_nphot());
    InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    plot_write(nfreq,InvCompton.get_energy_obs(),InvCompton.get_nphot_obs(),"Output/Singlezone_SSC.dat",1.,0.);
    
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
