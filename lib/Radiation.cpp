#include "Radiation.hpp"

#include <iostream>

//Methods to set viewing angle, beaming and geometry of emission region
void Radiation::set_beaming(double theta,double speed,double doppler){
	angle = theta*pi/180.;
	beta = speed;
	dopfac = doppler;
}

void Radiation::set_inclination(double theta){
	angle = theta*pi/180.;
}

void Radiation::set_geometry(bool sphsw,double l1,double l2){
	if (sphsw == false) {
		r = l1;
		z = l2;
		vol = pi*pow(r,2.)*z;
	} else {
		r = l1;
		z = r;
		vol = (4./3.)*pi*pow(r,3.);
	}
}

//Method to include a counterjet in cyclosycnchrotron/Compton classes
void Radiation::set_counterjet(bool flag){
	counterjet = flag;	
}

//simple method to check arrays; only meant for debugging
void Radiation::array_test(){
	for (int i=0;i<size;i++){
		std::cout << en_phot[i] << " "  << num_phot[i] << " " << num_phot[i]*en_phot[i]*herg << std::endl;
	}
}
