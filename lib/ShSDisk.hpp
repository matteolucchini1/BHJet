#ifndef SHSDISK_HPP
#define SHSDISK_HPP

#include "Radiation.hpp"
#include <iostream>

//Class Shakura-Sunyeav disk, inherited from Radiation.hpp

class ShSDisk: public Radiation {
    private:
        double Tin;							//Disk temperature at R=Rin in Kev
        double Mbh, Rg;						//Black hole mass in solar masses
        double Ldisk; 						//Disk luminosity in Eddington units
        double Hratio;						//Disk aspect ratio; TODO test with constant hbb or with hbb varying
									        //with disk size. Also test with just one Comptonization zone, and
									        //with one up to  hbb(Rin) and another up to the end of the nozzle

    public:		
        ~ShSDisk();	
        ShSDisk();

        const double hdisk()	{return Hratio;};
        const double tin() 		{return Tin;};
        const double rin()		{return r;};	
        const double lum()      {return Ldisk;};

        void disk_spectrum();
        void cover_disk(double f);
        friend double disk_int(double nu,void *p);

        double total_luminosity();
        
        void set_mbh(double M);
        void set_rin(double R);
        void set_rout(double R);
        void set_luminosity(double L);
        void set_tin(double T);

        void test();
};

#endif
