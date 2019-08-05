#include "Particles.hpp"

#include <iostream>

//Class destructor to de-allocate arrays
Particles::~Particles(){
	gsl_integration_workspace_free (w1);

	delete[] p;
	delete[] ndens;
	delete[] gamma;
	delete[] gdens;
	delete[] gdens_diff;
}

//Simple numerical integrals /w trapeze method
double Particles::count_particles(){
	double temp = 0.;
	for (int i=0;i<size-1;i++){
		temp = temp+(1./2.)*(p[i+1]-p[i])*(ndens[i+1]+ndens[i]);
	}
	return temp;
}

double Particles::count_particles_energy(){
	double temp = 0.;
	for (int i=0;i<size-1;i++){
		temp = temp+(1./2.)*(gamma[i+1]-gamma[i])*(gdens[i+1]+gdens[i]);
	}
	return temp;
}

double Particles::av_p(){
	double temp = 0.;
	for (int i=0;i<size-1;i++){
		temp = temp+(1./2.)*(p[i+1]-p[i])*(p[i+1]*ndens[i+1]+p[i]*ndens[i]);
	}
	return temp/count_particles();
}

double Particles::av_gamma(){
	return pow(pow(av_p()/(mass*cee),2.)+1.,1./2.);
}

double Particles::av_psq(){
	double temp = 0.;
	for (int i=0;i<size-1;i++){
		temp = temp+(1./2.)*(p[i+1]-p[i])*(pow(p[i+1],2.)*ndens[i+1]+pow(p[i],2.)*ndens[i]);
	}
	return temp/count_particles();
}

double Particles::av_gammasq(){
	return pow(av_psq()/pow(mass*cee,2.)+1.,1./2.);
}

//Methods to set up energy space number density, as a function of momentum space number density
void Particles::initialize_gdens(){
	for (int i=0;i<size;i++){
		gdens[i] = ndens[i]*gamma[i]*mass*cee/(pow(pow(gamma[i],2.)-1.,1./2.));
	}
}

//Same as above but the other way around
void Particles::initialize_pdens(){
	for (int i=0;i<size;i++){
		ndens[i] = gdens[i]*p[i]/(pow(mass*cee,2.)*pow(pow(p[i]/(mass*cee),2.)+1.,1./2.));
	}
}


void Particles::gdens_differentiate(){
	double *temp = new double[size];
	
	for (int i=0;i<size;i++){
		temp[i] = gdens[i]/(pow(gamma[i],2.)); 
	}
	
	for (int i=0;i<size-1;i++){
		gdens_diff[i] =  (temp[i+1]-temp[i])/(mass*pow(cee,2.)*(gamma[i+1]-gamma[i]));
	}
	
	gdens_diff[size-1] = gdens_diff[size-2];
	
	delete[] temp;
}

//simple method to check arrays; only meant for debugging
void Particles::array_test(){
	for (int i=0;i<size;i++){
		std::cout << p[i] << " "  << gamma[i] << " " << ndens[i] << " " << ndens[i]*p[i] << std::endl;
	}
}
