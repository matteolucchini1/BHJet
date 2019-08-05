#include "Radiation.hpp"
#include "Compton.hpp"

#include <iostream>

Compton::~Compton(){
	delete[] en_phot;
	delete[] num_phot;
	delete[] en_phot_obs;
	delete[] num_phot_obs;
	
	delete[] seed_energ;
	delete[] seed_urad;
	delete[] iter_urad;
		
	gsl_integration_workspace_free (w1);
	gsl_integration_workspace_free (w2);
	
	gsl_spline_free(seed_ph);
	gsl_interp_accel_free(acc_seed);
	
	gsl_spline_free(iter_ph);
	gsl_interp_accel_free(acc_iter);
}

//This constructor initializes the arrays and interpolations. In this case, calculations are done in frequency
//space, not in photon energies.
Compton::Compton(int s1,int s2,int n,double numin,double numax,double m){
	size = s1;
	seed_size = s2;

	en_phot = new double[size];
	num_phot = new double[size];
	iter_urad = new double[size];
	en_phot_obs = new double[2*size];
	num_phot_obs = new double[2*size];
	
	seed_energ = new double[seed_size];
	seed_urad = new double[seed_size]; 		
	
	Niter = n;	
	ypar = 1.;

	counterjet = false;

	seed_ph = gsl_spline_alloc(gsl_interp_steffen,seed_size);
	acc_seed = gsl_interp_accel_alloc();
	
	iter_ph = gsl_spline_alloc(gsl_interp_steffen,size);
	acc_iter = gsl_interp_accel_alloc();			
	
	w1 = gsl_integration_workspace_alloc (100);
	w2 = gsl_integration_workspace_alloc (100);
		
	double nuinc = (log10(numax)-log10(numin))/(size-1);
	
	for(int i=0;i<size;i++){
		en_phot[i] = pow(10.,log10(numin)+i*nuinc)*herg;
		en_phot_obs[i] = en_phot[i];
		en_phot_obs[i+size] = en_phot[i];
		num_phot[i] = 0;
		num_phot_obs[i] = 0;
		num_phot_obs[i+size] = 0;
		iter_urad[i] = 0;
	}	
	
	for(int i=0;i<seed_size;i++){
		seed_energ[i] = 0;
		seed_urad[i] = 0;
	}
	
	mass = m;	
}

//This function is the kernel of eq 2.48 in Blumenthal & Gould(1970), represents the scattered photon spectrum
//for a given electron and includes the Klein-Nishina cross section.
double comfnc(double ein,void *p){
    struct comfnc_params *params = (struct comfnc_params *)p;
    double game     = (params -> game);
    double e1       = (params -> e1);
    gsl_spline *phodis = (params -> phodis);
    gsl_interp_accel *acc_phodis = (params -> acc_phodis);
    
    double einit,utst,biggam,q;
    double phonum;
    double tm1,tm2,tm3;
    double btst,eg4;
  
    einit = exp(ein);
    btst = einit/(game*emerg);
    eg4	= 4.*einit*game;
    utst = eg4/(emerg+eg4);
    
    if((e1 < btst && (btst-e1)/btst >= 4.e-8) || (e1 > utst && (e1-utst)/utst >= 4.e-8)){
        return 0.;
    } else {
    	biggam = eg4/emerg;
    	q = e1/(biggam*(1.-e1));
    	phonum = gsl_spline_eval(phodis,einit,acc_phodis);
		tm1	= 2.*q*log(q);
    	tm2	= (1.+2.*q)*(1.-q);
    	tm3	= 0.5*(pow(biggam*q,2.)*(1.-q))/(1.+biggam*q);
		return (tm1+tm2+tm3)*pow(10.,phonum);
    }
}


//This function is the integral of comfnc above over the total seed photon distribution
double comint(double gam,void *p){
    struct comint_params *params = (struct comint_params *)p;
    double eph      = (params -> eph);
    double ephmin   = (params -> ephmin);
    double ephmax   = (params -> ephmax);
    gsl_spline *eldis = (params -> eldis);
    gsl_interp_accel *acc_eldis = (params -> acc_eldis);
    gsl_spline *phodis = (params -> phodis);
    gsl_interp_accel *acc_phodis = (params -> acc_phodis);
    gsl_integration_workspace *w = (params -> w);
    
    double game,econst,blim,ulim,e1,den;
    double result, error;

    game = exp(gam);

    e1 = eph/(game*emerg);

    econst = 2.*pi*re0*re0*cee;

    blim = log(std::max(eph/(4.*game*(game-eph/emerg)),ephmin));
    ulim = log(std::min(eph,ephmax));
	
    if(ulim <= blim){
        return 0;
    } else {    
    	gsl_function F2;
    	struct comfnc_params F2params = {game,e1,phodis,acc_phodis};
    	F2.function     = &comfnc;
    	F2.params       = &F2params;
    	gsl_integration_qag(&F2,blim,ulim,1e0,1e0,100,2,w,&result,&error);
    	den = gsl_spline_eval(eldis,game,acc_eldis);
		return econst*den*result/game;
    }
}

//This integrates the individual electron spectrum from comint over the total electron distribution
double Compton::comintegral(int it,double blim,double ulim,double enphot,double enphmin,double enphmax,gsl_spline 
							*eldis,gsl_interp_accel *acc_eldis){
    double result, error;
    gsl_function F1;
    struct comint_params F1params = {enphot,enphmin,enphmax,eldis,acc_eldis,seed_ph,acc_seed,w2};
    struct comint_params F1params_it = {enphot,enphmin,enphmax,eldis,acc_eldis,iter_ph,acc_iter,w2};
    F1.function = &comint;
    if (it == 0){
    	F1.params = &F1params;
    } else {
    	F1.params = &F1params_it;
    }

    gsl_integration_qag(&F1,blim,ulim,1e0,1e0,100,2,w1,&result,&error); 
    //NOTE: in some regimes, using a key of 2 in the gsl_integral_qag line instead of 1 makes for smoother 
    //integrals. Not important for the final spectrum, but it makes for better-looking and more accurate plots
    return result;
}

//This calculates the final spectrum in all frequency bins, including multiple scatters. 
//Note: the reason the Doppler boosting is a factor of 2 instead of 3 is because the calculations are done for
//a conical jet in the Lind&BLanford 1985 prescription
void Compton::compton_spectrum(double gmin,double gmax,gsl_spline *eldis,gsl_interp_accel *acc_eldis){
	double blim,ulim,com;
	double dopfac_cj;   	
   	double ephmin, ephmax;
   	
   	ephmin = seed_energ[0];
   	ephmax = seed_energ[seed_size-1];
   
    dopfac_cj = dopfac*(1.-beta*cos(angle))/(1.+beta*cos(angle));

	for (int it=0;it<Niter;it++){
		for(int i=0;i<size;i++){
			blim = log(std::max(gmin,en_phot[i]/emerg));
			ulim = log(gmax);
			if(blim >= ulim){					
       			com	= 1e-100;      			
       		} else {
   				com = comintegral(it,blim,ulim,en_phot[i],ephmin,ephmax,eldis,acc_eldis); 
       		}
			num_phot[i] = num_phot[i]+com*vol*en_phot[i]*herg;
			en_phot_obs[i] = en_phot[i]*dopfac;
			num_phot_obs[i] = num_phot[i]*pow(dopfac,2.);
			if(counterjet == true){ 
    			en_phot_obs[i+size] = en_phot[i]*dopfac_cj;
    			num_phot_obs[i+size] = num_phot[i]*pow(dopfac_cj,2.);
    		} else {
    			en_phot_obs[i+size] = 0;
    			num_phot_obs[i+size] = 0;
    		}
			if(com == 0){
				iter_urad[i] = -50;
			} else {
				iter_urad[i] = log10(com*vol/(pi*pow(r,2.)*cee));
			}
		}
		ephmin = en_phot[0];
		ephmax = en_phot[size-1];
		gsl_spline_init(iter_ph,en_phot,iter_urad,size);
	}
}

//Method to include cyclosynchrotron array from Cyclosyn.hh to the seed field. The two input arrays should be
//get_energ() and get_nphot() methods from Cyclosyn.hh. If the seed_urad array wasn't empty, the contribution 
//is automatically added to the existing one.
//note: the second condition has a <= sign to dodge numerical errors when the seed field energy density is
//extremely low, which can result in negative values/nan for the seed_urad 
void Compton::cyclosyn_seed(const double *seed_arr,const double *seed_lum){	
	for(int i=0;i<seed_size;i++){
		seed_energ[i] = seed_arr[i];		
		if (seed_urad[i] != 0) {
			seed_urad[i] = log10(pow(10.,seed_urad[i])+seed_lum[i]/(cee*herg*seed_energ[i]*pi*r*r));
		} else if (seed_lum[i]/(cee*herg*seed_energ[i]*pi*r*r) <= 0) {
			seed_urad[i] = -100;
		} else {
			seed_urad[i] = log10(seed_lum[i]/(cee*herg*seed_energ[i]*pi*r*r));	
		}
	}
	gsl_spline_init(seed_ph,seed_energ,seed_urad,seed_size);
}

//Method to include black body to seed field for IC; overloaded to initialize seed_freq if desired
//Note: Urad and Tbb need to be passed in the co-moving frame, the function does NOT account for beaming
void Compton::bb_seed(double Urad,double Tbb){
	double ulim,bbfield;

	ulim = 1e2*Tbb*kboltz_kev2erg;
	
	for(int i=0;i<seed_size;i++){
		if(seed_energ[i]<ulim){
			bbfield = (2.*Urad*pow(seed_energ[i]/herg,2.))/(herg*pow(cee,2.)*sbconst*
					  pow(Tbb*kboltz_kev2erg/kboltz,4)*(exp(seed_energ[i]/(Tbb*kboltz_kev2erg))-1.));
		} else {
			bbfield = 1.e-100;
		}
		if (seed_urad[i] != 0) {
			seed_urad[i] = log10(pow(10.,seed_urad[i]) + bbfield);
		} else if (bbfield > 0) {
			seed_urad[i] = log10(bbfield);
		} else {
			seed_urad[i] = -100;
		}
	}
	gsl_spline_init(seed_ph,seed_energ,seed_urad,seed_size);
}

void Compton::bb_seed(const double *seed_arr,double Urad,double Tbb){
	double ulim,bbfield;

	ulim = 1e2*Tbb*kboltz_kev2erg;
	
	seed_freq_array(seed_arr);	
	for(int i=0;i<seed_size;i++){
		if(seed_energ[i]<ulim){
			bbfield = (2.*Urad*pow(seed_energ[i]/herg,2.))/(herg*pow(cee,2.)*sbconst*
					  pow(Tbb*kboltz_kev2erg/kboltz,4)*(exp(seed_energ[i]/(Tbb*kboltz_kev2erg))-1.));
		} else {		
			bbfield = 1.e-100;
		}
		if (seed_urad[i] != 0) {
			seed_urad[i] = log10(pow(10.,seed_urad[i])+bbfield);
		} else if (bbfield > 0) {
			seed_urad[i] = log10(bbfield);
		} else {
			seed_urad[i] = -100;
		}	
	}
	gsl_spline_init(seed_ph,seed_energ,seed_urad,seed_size);
}

//Calculates incident photon field for a Shakura-Sunyaev disk with aspect ratio H, inner temperature Tin,
//inner radius Rin, at a distance z from the disk, taking beaming into account
double disk_integral(double alfa,void *p){
	struct disk_ic_params *params = (struct disk_ic_params *)p;
    double gamma = (params -> gamma);
    double beta = (params -> beta);
    double tin = (params -> tin);
    double rin = (params -> rin);
    double rout = (params -> rout);
    double h = (params -> h);
	double z = (params -> z);
    double nu = (params -> nu); 
    
    double delta = 1./(gamma-beta*cos(alfa));
    double a,b,x,y,r,T_eff,nu_eff,fac,Urad;
    
    a = z - (h*rout*rin)/(2.*(rout-rin));
    b = 1./tan(alfa) + (h*rout)/(2.*(rout-rin));
    
    x = a/b;
    y = x/tan(alfa) + z;

    r = pow(pow(x,2.)+pow(y,2.),1./2.);
	T_eff = delta*tin*pow(rin/r, 0.75);
	nu_eff = delta*nu;
	fac = (herg*nu_eff)/(kboltz*T_eff);
		
	Urad = 4.*pi*pow(nu_eff,2.)/(herg*pow(cee,3.)*(exp(fac)-1.));
	
	return Urad;
}

void Compton::shsdisk_seed(double tin,double rin,double rout,double h,double z){
	double ulim,blim,nulim,Gamma,result,error,diskfield;
	
	Gamma = 1./pow((1.-pow(beta,2.)),1./2.);

	blim = atan(rin/z);
	if(z<h*rout/2.){
		ulim = pi/2.+atan((h*rout/2.-z)/rout);
	} else {
		ulim = atan(rout/(z-h*rout/2.));
	}
	
	nulim = 1e1*tin*kboltz;
	
	for(int i=0;i<seed_size;i++){
		if (seed_energ[i] < nulim){	
			gsl_function F;
        	struct disk_ic_params Fparams = {Gamma,beta,tin,rin,rout,h,z,seed_energ[i]/herg};
        	F.function     = &disk_integral;
        	F.params       = &Fparams;       
        	gsl_integration_qag(&F,blim,ulim,0,1e-7,100,2,w1,&result,&error); 
			diskfield = result;
		} else {
			diskfield = 1.e-100;
		}
		if (seed_urad[i] != 0) {
			seed_urad[i] = log10(pow(10.,seed_urad[i])+diskfield);
		} else if (diskfield > 0) {
			seed_urad[i] = log10(diskfield);	
		} else {
			seed_urad[i] = -100;
		}	
	}
	gsl_spline_init(seed_ph,seed_energ,seed_urad,seed_size);	
}

void Compton::shsdisk_seed(const double *seed_arr,double tin,double rin,double rout,double h,double z){
	double ulim,blim,nulim,Gamma,result,error,diskfield;
	
	Gamma = 1./pow((1.-pow(beta,2.)),1./2.);

	blim = atan(rin/z);
	if(z<h*rout/2.){
		ulim = pi/2.+atan((h*rout/2.-z)/rout);
	} else {
		ulim = atan(rout/(z-h*rout/2.));
	}
	//std::cout<<blim<<" "<<ulim<<std::endl;
	nulim = 1e1*tin*kboltz;
	seed_freq_array(seed_arr);
	
	for(int i=0;i<seed_size;i++){
		if (seed_energ[i] < nulim){	
			gsl_function F;
        	struct disk_ic_params Fparams = {Gamma,beta,tin,rin,rout,h,z,seed_energ[i]/herg};
        	F.function     = &disk_integral;
        	F.params       = &Fparams;       
        	gsl_integration_qag(&F,blim,ulim,0,1e-5,100,2,w1,&result,&error);        
			diskfield = result;
		}else {
			diskfield = 1.e-100;
		}   
		if (seed_urad[i] != 0) {
			seed_urad[i] = log10(pow(10.,seed_urad[i])+diskfield);
		} else if (diskfield > 0) {
			seed_urad[i] = log10(diskfield);
		} else {
			seed_urad[i] = -100;
		}
	}
	gsl_spline_init(seed_ph,seed_energ,seed_urad,seed_size);	
}

//Method to estimate the number of scatterings the electrons go through; scatters are repeated until the 
//photons reach the same typical energy as the electrons
//Input parameters are initial scale frequency of photons in Hz and temperature of scattering electrons in erg
//note: this assumes that a) the final energy of the photons is set only by the electrons and b) that the 
//photons have one scale energy only; it could be incorrect if scattering multiple photon fields (e.g. disk 
//and cyclosynchrotron) 
void Compton::set_niter(double nu0,double Te){
	double x0,xf,k;

	x0 = herg*nu0/emerg;
	xf = Te/emerg;
	k = log(xf/x0);///(16.*pow(xf,2.)+4.*xf);
	Niter = int(k+0.5);
}

void Compton::set_niter(int n){
	Niter = n;
}

//Calculates optical depth for given number density and radius of emitting region, and compton-Y for a given
//electron average Lorentz factor
void Compton::set_tau(double n,double r,double gam){
	tau = n*r*sigtom;
	ypar = std::max(tau*tau,tau)*(pow(4.0*gam,1)+pow(4.0*gam,2));
}

//Method to set up seed frequency array
void Compton::seed_freq_array(const double *seed_arr){
	for(int i=0;i<seed_size;i++){
		seed_energ[i] = seed_arr[i];
	}	
}

//This resets the energy density in case different photon fields want to be calculated separately to see the
//contribution of each
void Compton::reset(){
	for (int i=0;i<seed_size;i++){
		seed_urad[i] = 0;
		seed_energ[i] = 0;
	}
	for (int i=0;i<size;i++){
		num_phot[i] = 0;
		num_phot_obs[i] = 0;
	}
}

void Compton::urad_test(){
	for (int i=0;i<seed_size;i++){
		std::cout<<seed_energ[i]/herg<<" "<<seed_urad[i]<<" "<<iter_urad[i]<<std::endl;
	}
}

void Compton::test(){
	std::cout << "Optical depth: " << tau << " Compton-Y: " << ypar <<  " r: " << r << " z: " << z <<
	 " angle: " << angle << " speed: " << beta << " delta: " << dopfac << std::endl;	
	std::cout << "Number of scatters: " << Niter << std::endl;
	std::cout << "Optically thin slope: " << -log10(tau)/(log10(ypar)-log10(tau)) << std::endl;
	std::cout << "Kompaneets' slope: " << -3./2.+pow(9./4.+4./ypar,1./2.) << std::endl;
}

