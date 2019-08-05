#ifndef THERMAL_HPP
#define THERMAL_HPP

#include "Particles.hpp"
#include <iostream>

//Class for thermal particles, inherited from the generic Particles class in Particles.hpp
//note: ndens is number density per unit momentum

class Thermal: public Particles {
	private:
		double Temp, thnorm, theta;
	public:
		Thermal(int s,int type, double T);		
			
		void set_p();						
		void set_ndens();
		void set_temp(double T);
		void set_norm(double n);
		
		double K2(double x);
		
		void test();		
};

#endif
