#ifndef COMPTON_HPP
#define COMPTON_HPP

//Class inverse Compton, inherited from Radiation.hpp

class Compton: public Radiation {
	private:
		int seed_size;					//size of seed photon field
		int Niter;						//number of IC iterations
		double mass;					//Mass of radiating particles
		double tau,ypar;				//optical depth/comtpon Y of emitting region
		
		double *seed_energ;				//array of seed frequencies in Hz
		double *seed_urad; 				//array of seed photon number density in log10(#/erg/cm^3)
		double *iter_urad; 				//array of iterated photon number density in log10(#/erg/cm^3)
				
		gsl_spline *seed_ph;			//interpolation of photon field array seed_urad
		gsl_interp_accel *acc_seed;		//accelerator for above spline
		
		gsl_spline *iter_ph;			//interpolation of photon field for multiple scatters
		gsl_interp_accel *acc_iter;		//accelerator of above spline
		
		gsl_integration_workspace *w2;	//workspace for second numerical integration		
	public:
		~Compton();
		Compton(int s1,int s2,int Niter,double numin,double numax,double m);	
		
		friend double comfnc(double ein,void *p);
		friend double comint(double gam,void *p);
		friend double disk_integral(double alfa,void *p);
		double comintegral(int it,double blim,double ulim,double nu,double numin, double numax,gsl_spline 
						   *eldis,gsl_interp_accel *acc_eldis);
		void compton_spectrum(double gmin,double gmax,gsl_spline *eldis,gsl_interp_accel *acc_eldis);
		//double compton_spectrum_parallel();//same as above but handle integration workspaces differently
		
		void cyclosyn_seed(const double *seed_energ,const double *seed_lum);
		void bb_seed(const double *seed_energ,double Urad,double Tbb);
		void bb_seed(double Urad,double Tbb);	
		void shsdisk_seed(double tin,double rin,double rout,double h,double z);		
		void shsdisk_seed(const double *seed_arr,double tin,double rin,double rout,double h,double z);
		
		void set_tau(double n,double r,double gam);		
		void set_niter(double nu0,double Te);
		void set_niter(int n);
		void seed_freq_array(const double *seed_energ);
		
		const double get_tau(){return tau;};
		const double get_ypar(){return ypar;};
		
		void reset();
		void urad_test();
		void test();
};

#endif
