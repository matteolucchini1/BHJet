#ifndef BBODY_HPP
#define BBODY_HPP

#include <iostream>
#include "Radiation.hpp"

//Class black body photons, inherited from Radiation.hh

class BBody: public Radiation {
    private:
        double Tbb;
        double Lbb;
        double normbb;
    public:		
        ~BBody();
        BBody(double T,double L);

        void set_temp_kev(double T);
        void set_temp_k(double T);
        void set_temp_hz(double nu);
        void set_lum(double L);
        void bb_spectrum();

        const double temp_kev();
        const double temp_k();
        const double temp_hz();
        const double lum();
        const double norm();
        const double Urad(double d);

        void test();
};

#endif
