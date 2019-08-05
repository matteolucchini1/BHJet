#include "ShSDisk.hpp"

ShSDisk::~ShSDisk(){
	delete[] en_phot;
	delete[] num_phot;
	delete[] en_phot_obs;	
	delete[] num_phot_obs;
	
	gsl_integration_workspace_free (w1);
}

//Constructor for the disk. If Tflag is true, then Tin is calculated as a function of input luminosity in 
//Eddington luminosity, otherwise Tin is taken to be the inner disk temperature in Kev, from which the disk
//luminosity can be calculated
ShSDisk::ShSDisk(bool Tflag,double M,double T,double R1,double R2){
	double emin,emax,einc;
	double eta;			//estimated accretion efficiency as a function of inner radius to calculate luminosity

	size = 50;

	en_phot = new double[size];
	num_phot = new double [size];
	en_phot_obs = new double [size];
	num_phot_obs = new double [size];

	Mbh = M;
	Rg = gconst*Mbh*msun/(cee*cee);
	r = R1;					//note: for disk r is inner radius, z is the outer radius
	z = R2;					//both are input in cm
	eta = 0.04;				//note: this efficiency is tuned so that the estimated and total calculated 
							//luminosity of the disk are the same.

	if (Tflag == true){
		Eddrat = T;
		Tin =  pow(Eddrat*1.25e38*Mbh/(16.*pi*eta*sbconst*pow(r,2.)),0.25);
	} else {
		Tin = T*kboltz_kev2erg/kboltz;
		Eddrat = 16.*pi*eta*sbconst*pow(Tin,4.)*pow(r,2.)/(1.25e38*Mbh);
	}
		
	Hratio = std::max(0.1,Eddrat);

	emin = 0.0005*kboltz*Tin;
	emax = 30.*kboltz*Tin;
	
	einc = (log10(emax)-log10(emin))/(size-1);
	
	for(int i=0;i<size;i++){
		en_phot_obs[i] = pow(10.,log10(emin)+i*einc);
		en_phot[i] = en_phot_obs[i];
		num_phot[i] = 0.;
		num_phot_obs[i] = 0.;
	}
	
	w1 = gsl_integration_workspace_alloc (100);
}	

//return SD spectrum over a given radius, frequency to be integrated over radius
double disk_int(double lr,void *p){
	struct disk_obs_params *params = (struct disk_obs_params *)p;
    double tin      = (params -> tin);
    double rin      = (params -> rin);
    double nu       = (params -> nu);
    
    double r,temp,fac,bb;
    
    r = exp(lr);
    temp = tin*pow(rin/r,0.75);
    fac	= herg*nu/(kboltz*temp);
    
    if(fac<1.e-3){
        bb	= 2.*herg*pow(nu,3.)/(pow(cee,2.)*fac);
    }
    else{
        bb	= 2.*herg*pow(nu,3.)/(pow(cee,2.)*(exp(fac)-1.));
    }
    
    return 2.*pi*pow(r,2.)*bb;
}

void ShSDisk::disk_spectrum(){
	double result, error;

	for(int k=0;k<size;k++){
		gsl_function F1;
    	struct disk_obs_params F1params = {Tin,r,en_phot_obs[k]/herg};
    	F1.function = &disk_int;
    	F1.params = &F1params;
    	gsl_integration_qag(&F1,log(r),log(z),0,1e-2,100,2,w1,&result,&error);
		num_phot[k] = result;
    	num_phot_obs[k] = cos(angle)*result;
	}
}

//Simple integration method to integrate num_phot_obs and get the total luminosity of the disk
double ShSDisk::total_luminosity(){
	double temp = 0.;
	for (int i=0;i<size-1;i++){
		temp = temp+(1./2.)*(en_phot_obs[i+1]/herg-en_phot_obs[i]/herg)*
			   (num_phot_obs[i+1]/cos(angle)+num_phot_obs[i]/cos(angle));
	}
	return temp;
}

void ShSDisk::test(){
	std::cout << "Inner disk radius: " << r << " cm, " << r/Rg << " rg; outer disk radius: " << z <<
			     " cm, " << z/Rg << " rg; disk scale height: " << Hratio << std::endl;	
	std::cout << "Inner disk temperature: " << Tin*kboltz/kboltz_kev2erg << 
			     " kev; Disk luminosity in Eddington units: " << Eddrat << " and erg s^-1: " << 
			     Eddrat*1.25e38*Mbh << std::endl;
}
