#ifndef KAPPA_HPP
#define KAPPA_HPP

//Class for kappa distribution of particles, inherited from the generic Particles class in Particles.hh

class Kappa: public Particles {
    private:
        double theta;
        double kappa, knorm;
        double pmin,pmax;
    public:	
        Kappa(int s,int type,double T,double k);		
	
        void set_p(double ucom,double bfield,double betaeff,double r,double fsc);
        void set_p(double gmax);					
        void set_ndens();
        void set_kappa(double k);
        void set_temp(double T);
        void set_norm(double n);

        friend double norm_kappa_int(double x,void *p);
        friend double injection_kappa_int(double x,void *p);

        void cooling_steadystate(double ucom, double n0,double bfield,double r,double gshift);
        double max_p(double ucom,double bfield,double betaeff,double r,double fsc);

        void test();
};

#endif
