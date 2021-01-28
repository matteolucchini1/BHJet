#ifndef SHSDISK_HPP
#define SHSDISK_HPP

#include "Radiation.hpp"
#include <iostream>

//Class Shakura-Sunyeav disk, inherited from Radiation.hpp

class ShSDisk: public Radiation {
    private:
        double Tin;							//Disk temperature at R=Rin in Kev
        double Mbh, Rg;						//Black hole mass in solar masses
        double Eddrat;						//Disk luminosity in Eddington units
        double Hratio;						//Disk aspect ratio; TODO test with constant hbb or with hbb varying
									        //with disk size. Also test with just one Comptonization zone, and
									        //with one up to  hbb(Rin) and another up to the end of the nozzle

    public:		
        ~ShSDisk();	
        ShSDisk(bool Tflag,double M,double T,double R1,double R2);

        const double hdisk()	{return Hratio;};
        const double tin() 		{return Tin;};
        const double rin()		{return r;};	
        const double Eddratio() {return Eddrat;};

        void disk_spectrum();
        friend double disk_int(double nu,void *p);

        double total_luminosity();
        double scale_freq();

        void test();
};

#endif
