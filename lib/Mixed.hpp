#ifndef MIXED_HPP
#define MIXED_HPP

//Class for mixed particles, inherited from the generic Particles class in Particles.hh
//the minimum momentum of the PL is always be assumed to be the averge momentum of the thermal
//note: ndens is number density per unit momentum

//The reasons the integrands are friends is kind of magic. Basically, we need pointers to the private members
//Temp/thnorm/theta/whatever to set the parameters of the integrand functions in the GSL libraries. In C++
//this is not possible for non-static private member functions (because otherwise we could access private
//members from the outside by using a pointer). To quote the compiler error, C++ forbids
//taking the address of an unqualified or parenthesized non-static member function to form a pointer to
//member function. And we need a pointer member function void *p to set up the integrands for GSL libraries.
//By using a friend function we can instead set up a set of parameters whose value is that of a pointer which 
//points at the values stored in the private members; by doing this, we can not change the private members' 
//values (correctly so) but we can still access their numerical value and use it elsewhere.

class Mixed: public Particles {
    private:
        double thnorm, theta;
        double pspec, plnorm;
        double pmin_th, pmax_th, pmin_pl, pmax_pl;
        double plfrac;
    public:
        Mixed(int s,int type,double T,double s1,bool flag);			

        void set_p(double ucom,double bfield,double betaeff,double r,double fsc);
        void set_p(double gmax);	
        void set_ndens();
        void set_temp(double T);
        void set_norm(double n);
        void set_plfrac(double f);				
        void set_pspec(double s1);	

        void cooling_steadystate(double ucom,double n0,double bfield,double r,double betaeff);
        double max_p(double ucom,double bfield,double betaeff,double r,double fsc);

        friend double th_num_dens_int(double x,void *p);
        friend double av_th_p_int(double x,void *p);	
        double count_th_particles();
        double av_th_p();
        double av_th_gamma();	
	
        friend double pl_num_dens_int(double x,void *p);
        friend double av_pl_p_int(double x,void *p);
        double count_pl_particles();
        double av_pl_p();	
        double av_pl_gamma();

        double K2(double x);		

        void test();
};

#endif
