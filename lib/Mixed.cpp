#include "Particles.hpp"
#include "Mixed.hpp"

#include <iostream>

//Class constructors to initialize object
Mixed::Mixed(int s,int type,double T,double s1,bool flag){
	size = s;

	p = new double[size];
	ndens = new double[size];	
	gamma = new double[size];
	gdens = new double[size];
	gdens_diff = new double[size];
	
	w1 = gsl_integration_workspace_alloc (100);
	
	if (type==1) {mass = emgm;}
	else  {mass = pmgm;}

	theta = T/(mass*cee*cee);
	thnorm = 1.;										
	
	double emin_th = (1./100.)*(T/kboltz_kev2erg);	//minimum energy in kev, 1/50 lower than peak
	double emax_th = 20.*(T/kboltz_kev2erg); 		//maximum energy in kev, 20 higher than peak
	double gmin_th, gmax_th;
	
	if (mass == emgm){
		gmin_th = emin_th/me_kev+1.;
		gmax_th = emax_th/me_kev+1.;
	} else if (mass == pmgm) {
		std::cout << "Error! Protons currently unsupported!" << std::endl;
	}
	
	pmin_th = pow(pow(gmin_th,2.)-1.,1./2.)*mass*cee;
	pmax_th = pow(pow(gmax_th,2.)-1.,1./2.)*mass*cee;
	
	pspec = s1;
	plnorm = 1.;	
	
	FP_flag = flag;	

	for (int i=0;i<size;i++){
		p[i] = 0;
		ndens[i] = 0;
		gamma[i] = 0;
		gdens[i] = 0;
		gdens_diff[i] = 0;
	}
}


//Methods to set momentum/energy arrays and number density arrays
void Mixed::set_p(double ucom,double bfield,double betaeff,double r,double fsc){	
	pmin_pl = av_th_p();
	pmax_pl = std::max(max_p(ucom,bfield,betaeff,r,fsc),pmax_th);	

	double pinc = (log10(pmax_pl)-log10(pmin_th))/size;
	
	for (int i=0;i<size;i++){
		p[i] = pow(10.,log10(pmin_th)+i*pinc);
		gamma[i] = pow(pow(p[i]/(mass*cee),2.)+1.,1./2.);
	}
}

//Same as above, but assuming a fixed maximum Lorentz factor
void Mixed::set_p(double gmax){
	pmin_pl = av_th_p();
	pmax_pl = pow(pow(gmax,2.)-1.,1./2.)*mass*cee;
	
	double pinc = (log10(pmax_pl)-log10(pmin_th))/size;
	
	for (int i=0;i<size;i++){
		p[i] = pow(10.,log10(pmin_th)+i*pinc);
		gamma[i] = pow(pow(p[i]/(mass*cee),2.)+1.,1./2.);
	}	
}

void Mixed::set_ndens(){
	for (int i=0;i<size;i++){
		if (p[i] <= pmin_pl){
			ndens[i] = thnorm*pow(p[i],2.)*exp(-gamma[i]/theta);
		} else if (p[i] < pmax_th) {
			ndens[i] = thnorm*pow(p[i],2.)*exp(-gamma[i]/theta) + plnorm*pow(p[i],-pspec)*exp(-p[i]/pmax_pl);
		} else {
			ndens[i] = plnorm*pow(p[i],-pspec)*exp(-p[i]/pmax_pl);
		}
	}
	if (FP_flag == false){
		initialize_gdens();
		gdens_differentiate();	
	}
}

//methods to set the temperature, pl fraction, and normalizations. Temperature must be in ergs, no factor kb
void Mixed::set_temp(double T){
	theta = T/(mass*cee*cee);
	
	double emin_th = (1./100.)*(T/kboltz_kev2erg);	//minimum energy in kev, 1/50 lower than peak
	double emax_th = 20.*(T/kboltz_kev2erg); 		//maximum energy in kev, 20 higher than peak
	double gmin_th, gmax_th;
	
	if (mass == emgm){
		gmin_th = emin_th/me_kev+1.;
		gmax_th = emax_th/me_kev+1.;
	} else if (mass == pmgm) {
		std::cout << "Error! Protons currently unsupported!" << std::endl;
	}
	
	pmin_th = pow(pow(gmin_th,2.)-1.,1./2.)*mass*cee;
	pmax_th = pow(pow(gmax_th,2.)-1.,1./2.)*mass*cee;	
}

void Mixed::set_pspec(double s1){
	pspec = s1;
}

void Mixed::set_plfrac(double f){		
	plfrac = f;
}

void Mixed::set_norm(double n){
	thnorm = (1.-plfrac)*n/(pow(mass*cee,3.)*theta*K2(1./theta));
	plnorm = plfrac*n*(1.-pspec)/(pow(pmax_pl,(1.-pspec))-pow(pmin_pl,(1.-pspec)));	
}

//Injection function to be integrated in cooling
double injection_mixed_int(double x,void *p){	
	struct injection_mixed_params * params = (struct injection_mixed_params*)p; 	 
	
	double s = (params->s);
	double t = (params->t); 
	double nth = (params->nth);
	double npl = (params->npl);  
	double m = (params->m);
	double min = (params->min);
	double max = (params->max);	
	double cutoff = (params->cutoff);

	double mom_int = pow(pow(x,2.)-1.,1./2.)*m*cee;	

	if (x<=min) {
		return  nth*pow(mom_int,2.)*exp(-x/t);
	} else if (x < max) {
		return nth*pow(mom_int,2.)*exp(-x/t) + npl*pow(mom_int,-s)*exp(-mom_int/cutoff);
	} else {
		return npl*pow(mom_int,-s)*exp(-mom_int/cutoff);
	} 
}

//Method to solve steady state continuity equation. NOTE: KN cross section not included in IC cooling
void Mixed::cooling_steadystate(double ucom, double n0,double bfield,double r,double betaeff){
	double Urad = pow(bfield,2.)/(8.*pi)+ucom;
	double pdot_ad = betaeff*cee/r;
	double pdot_rad = (4.*sigtom*cee*Urad)/(3.*mass*pow(cee,2.));
	double tinj = r/(cee);

	double gam_min, gam_max;
	gam_min = pow(pow(pmin_pl/(mass*cee),2.)+1.,1./2.);
	gam_max = pow(pow(pmax_th/(mass*cee),2.)+1.,1./2.);		
	
	gsl_function F1;	
	struct injection_mixed_params params = {pspec,theta,thnorm,plnorm,mass,gam_min,gam_max,pmax_pl};
	F1.function = &injection_mixed_int;
	F1.params   = &params;
	
	double integral, error;

	for (int i=0;i<size;i++){
		if (i < size-1) {
			gsl_integration_qag(&F1, gamma[i], gamma[i+1], 1e1, 1e1, 100, 1, w1, &integral, &error);
			ndens[i] = (integral/tinj)/(pdot_ad*p[i]/(mass*cee)+pdot_rad*pow(p[i]/(mass*cee),2.));
		}
		else {
			ndens[size-1] = ndens[size-2]*pow(p[size-1]/p[size-2],-pspec-1);
		}
	}		
	// the last bin is set by arbitrarily assuming cooled distribution; this is necessary because the integral 
	//above is undefined for the last bin

	//The last step requires a renormalization. The reason is that the result of gsl_integration_qag strongly
	//depends on the value of "size". Without doing anything fancy, this can be fixed simply by ensuring that
	//the total integrated number of density equals n0 (which we know), and rescaling the array ndens[i] by
	//the appropriate constant. If this explanation is unclear, comment out the following 4 lines of code and
	//plot ndens[i] for different values of "size".
	double renorm = count_particles()/n0;
		
	for (int i=0;i<size;i++){
		ndens[i] = ndens[i]/renorm;
	}
	
	initialize_gdens();
	gdens_differentiate();	
}

//Method to calculate maximum momentum of non thermal particles based on acceleration and cooling timescales
//The estimate is identical to the old agnjet but in momentum space; see Lucchini et al. 2019 for the math of
//the old version
double Mixed::max_p(double ucom,double bfield,double betaeff,double r,double fsc){
	double Urad, escom, accon, syncon, b, c, gmax;
	Urad = pow(bfield,2.)/(8.*pi)+ucom;
	escom = betaeff*cee/r;
	syncon = (4.*sigtom*Urad)/(3.*mass*cee);
	accon = (3.*fsc*charg*bfield)/(4.*mass*cee);
	
	b = escom/syncon;
	c = accon/syncon;
	
	gmax = (-b+pow(pow(b,2.)+4.*c,1./2.))/2.;

	return pow(pow(gmax,2.)-1.,1./2.)*mass*cee;
}

//Evaluate Bessel function as in old agnjet
double Mixed::K2(double x){
	double res;
	
	if (x < 0.1){
		res = 2./x/x;
	}
	else {
		res = gsl_sf_bessel_Kn(2,x);
	}
	
	return res;
}

//TEST METHODS: THESE ARE ONLY MEANT FOR DEBUGGING

//Methods to calculate number density and average energy in thermal part
double th_num_dens_int(double x,void *p){
	struct th_params * params = (struct th_params*)p; 	 
	
	double t = (params->t); 
	double n = (params->n);  
	double m = (params->m); 	
	 
	double gam_int = pow(pow(x/(m*cee),2.)+1.,1./2.);
	
   	return  n*pow(x,2.)*exp(-gam_int/t);
}

double av_th_p_int(double x,void *p){
	struct th_params * params = (struct th_params *)p; 	 
	
	double t = (params->t); 
	double n = (params->n);  
	double m = (params->m); 	
	 
	double gam_int = pow(pow(x/(m*cee),2.)+1.,1./2.);
	
   	return  n*pow(x,3.)*exp(-gam_int/t);
}

double Mixed::count_th_particles(){
	
	double integral1, error1;
	
	gsl_function F1;	
	struct th_params params = {theta,thnorm,mass};
	F1.function = &th_num_dens_int;
	F1.params   = &params;
	
	gsl_integration_qag(&F1,pmin_th,pmax_th,0,1e-7,100,1,w1,&integral1,&error1);
	
	return integral1;
}

double Mixed::av_th_p(){

	double integral1, error1, integral2;
	
	gsl_function F1;	
	struct th_params params = {theta,thnorm,mass};

	F1.function = av_th_p_int;
	F1.params   = &params;
	
	gsl_integration_qag(&F1,pmin_th,pmax_th,0, 1e-7,100,1,w1,&integral1,&error1);
	
	integral2 = count_th_particles();

	return integral1/integral2;
}

double Mixed::av_th_gamma(){
	double avp;
		
	avp=av_th_p();

	return pow(pow(avp/(mass*cee),2.)+1.,1./2.);
}

//Methods to calculate number density and average energy in non-thermal part
double pl_num_dens_int(double x,void *p){
	struct pl_params * params = (struct pl_params *)p; 	 
	
	double s = (params->s); 
	double n = (params->n);  
	
   	return  n*pow(x,-s);
}

double av_pl_p_int(double x,void *p){
	struct pl_params * params = (struct pl_params *)p; 	 
	
	double s = (params->s); 
	double n = (params->n);  
	
   	return  n*pow(x,-s+1.); 
}

double Mixed::count_pl_particles(){
	
	double integral1, error1;
	
	gsl_function F1;	
	struct pl_params params = {pspec,plnorm};
	F1.function = &pl_num_dens_int;
	F1.params   = &params;
	
	gsl_integration_qag(&F1,pmin_pl,pmax_pl,0,1e-7,100,1,w1,&integral1,&error1);
	
	return integral1;
}

double Mixed::av_pl_p(){
	
	double integral1, error1, integral2;
	
	gsl_function F1;	
	struct pl_params params =  {pspec,plnorm};
	F1.function = &av_pl_p_int;
	F1.params   = &params;
	
	gsl_integration_qag(&F1,pmin_pl,pmax_pl,0,1e-7,100,1,w1,&integral1,&error1);
	
	integral2 = count_pl_particles();

	return integral1/integral2;
}

double Mixed::av_pl_gamma(){
	double avp;
		
	avp=av_pl_p();

	return pow(pow(avp/(mass*cee),2.)+1.,1./2.);
}

//simple method to check quantities.
void Mixed::test(){
	std::cout << "Mixed distribution;" << std::endl;
	std::cout << "Temperature: " << theta*mass*cee*cee/kboltz_kev2erg << " kev" << std::endl;
	std::cout << "Number density: " << count_particles() << std::endl;	
	std::cout << "Thermal monetum limits: " << pmin_th << " " << pmax_th << std::endl;
	std::cout << "Non-thermal momentum limits: " << pmin_pl << " " << pmax_pl << std::endl;
}

