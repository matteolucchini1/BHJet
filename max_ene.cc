#include "agnjet.hh"

using namespace std;

void maxene(double z,double r,double ucom,double fsc,double bfield,double beta,double delz,double brk,
double &emax,double &gemax,double &ebreak,double &gebreak,double &bemax){
    
    double accon,syncon,comcon,escom,tdyn,qutrmb,qutrmc,ub;
    
    //Finding the maximum particle energy by solving: t^{-1}_acc = t^{-1}_syn + t^{-1}_compton + t^{-1}_esc
	//This is equivalent to solving a second degree equation in E: 
    //E**2 + escon/(syncon+comcon)*E - accon/(syncon+comcon)
    ub      = bfield*bfield/(8.*pi);    
    accon	= 3./4.*fsc*cee*charg*bfield;
    syncon	= 4./3.*sigtom*ub/(emgm*emgm*cee*cee*cee);
    comcon	= syncon*ucom/ub;    
    tdyn    = (brk*r)/(beta*cee);
    escom	= 1./tdyn;
    qutrmb	= escom/(syncon+comcon);
    qutrmc	= accon/(syncon+comcon);
    emax	= (-qutrmb+sqrt(qutrmb*qutrmb+4.*qutrmc))/2.;
    //Calculation of synchrotron cooling break energy - do this by setting tsyn = tdyn, where the electron's 
    //dynamical time is proportional to the time it takes them to escape the current jet segment  
    ebreak  = 1. / (syncon * tdyn);
    gemax	= emax/emerg;
    gebreak = ebreak/emerg;

	if(gemax >= 100){
        bemax	= 1. - 1./(2.*gemax*gemax) - 1./(8.*gemax*gemax*gemax*gemax);
    }
    else{
        bemax	= sqrt(gemax*gemax - 1)/gemax;
    }
}
