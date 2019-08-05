#ifndef Bknpower_HPP
#define Bknpower_HPP

//Class for non-thermal particles, inherited from the generic Particles class in Particles.hpp
//note: ndens is number density per unit momentum

class Bknpower: public Particles {
	private:
		double pspec1, pspec2, norm;
		double pmin, pbrk, pmax;
	public:
		Bknpower(int s,int type,double s1,double s2,bool flag);		
		
		void set_p(double min,double brk,double ucom,double bfield,double tshift,double bjet,double r,
				   double fsc);
		void set_p(double min,double brk,double gmax);	
		void set_ndens();
		void set_pspec1(double s1);
		void set_pspec2(double s2);
		void set_brk(double brk);
		void set_norm(double n);	
		
		friend double norm_bkn_int(double x,void *p);
		friend double injection_bkn_int(double x,void *p);		
		
		void cooling_steadystate(double ucom,double n0,double bfield,double r,double betaeff);
		double max_p(double ucom,double bfield,double tshift,double bjet,double r,double fsc);
		
		void test();	
};

#endif
