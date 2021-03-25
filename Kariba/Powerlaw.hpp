#ifndef Powerlaw_HPP
#define Powerlaw_HPP

//Class for non-thermal particles, inherited from the generic Particles class in Particles.hpp
//note: ndens is number density per unit momentum

class Powerlaw: public Particles {
    private:
        double pspec, plnorm;
        double pmin, pmax;
    public:
        Powerlaw(int s);		

        void set_p(double min,double ucom,double bfield,double betaeff,double r,double fsc);
        void set_p(double min,double gmax);		
        void set_ndens();
        void set_pspec(double s1);
        void set_norm(double n);			

        void cooling_steadystate(double ucom,double n0,double bfield,double r,double tshift);
        double max_p(double ucom,double bfield,double betaeff,double r,double fsc);

        void test();	
};

#endif
