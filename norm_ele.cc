#include "agnjet.hh"

using namespace std;

void elnorm(int flagNorm,double mxsw,double ntot,double pspec,double plfrac,double eltemp,double gshift,
double emax,double ebreak,double &mjteff,double &emin,double &cnorm,double &enorm,double &betat){
    
    double mjpeak;
    
    //From emax we can generate powerlaw between min and emax. Use emin normalized to scaled thermal peak at
    //shock; after shock then scale the two down together (see p.158, ntbk 6). This is a bit below the peak
    //for lower temperatures
    mjteff	= eltemp*gshift;    
    emin	= 2.23*kboltz*mjteff;
	
    //Creating powerlaw of normalization plfrac*ntot - 3 options based on where the cooling break energy lies
    
    if (flagNorm == 1){    
        if(ebreak <= emin) {//if break energy less than minimum, set normalisation for cooled distribution
            cnorm = plfrac*ntot*(pspec)/(pow(emin, -pspec) - pow(emax, -pspec));  
        }
        else if(ebreak >= emax){//if break greater than max, set normalisation for uncooled distribution
            cnorm	= plfrac*ntot*(1.-pspec)/(pow(emax, 1.-pspec) - pow(emin, 1.-pspec)); 
        }
        else{//If break lies within energy limits, set normalisation for broken lower law
            cnorm = plfrac*ntot/((pow(ebreak,1.-pspec)/(1.-pspec))-(pow(emin,1.-pspec)/(1.-pspec))+
            (pow(ebreak,1.-pspec)/pspec)-(ebreak*pow(emax, -pspec)/pspec));       
    	}    
    } else {        
        if(ebreak <= emin) {//if break energy less than minimum, set normalisation for cooled distribution
            cnorm = (1.-mxsw)*ntot*(pspec)/(pow(emin, -pspec) - pow(emax, -pspec));  
        }
        else if(ebreak >= emax){//if break greater than max, set normalisation for uncooled distribution
            cnorm	= (1.-mxsw)*ntot*(1.-pspec)/(pow(emax, 1.-pspec) - pow(emin, 1.-pspec)); 
        }
        else{//If break lies within energy limits, set normalisation for broken lower law
            cnorm = (1.-mxsw)*ntot/((pow(ebreak,1.-pspec)/(1.-pspec))-(pow(emin,1.-pspec)/(1.-pspec))+
            (pow(ebreak,1.-pspec)/pspec)-(ebreak*pow(emax,-pspec)/pspec));
        }

    }

    //Below is full form, but extra factor of emerg cancels out
    betat	= emerg/(kboltz*mjteff);
    mjpeak	= pow(10.,9.514 - 1.028*log10(mjteff));

    if(ebreak <= emin) {//if break energy < min energy, normalisation at boundary between thermal/power law 
    					//distributions is product of thermal and cooled power law
        enorm	= cnorm*pow(emin,-(pspec+1.))*betat/(k2_fnc(betat)*mjpeak); 
    }
    else{// if break energy > min energy, normalisation at boundary between thermal/power law distributions is 
    	 //product of thermal and uncooled power law
        enorm = cnorm*pow(emin, -pspec)*betat/(k2_fnc(betat)*mjpeak); 
    }    
}


void norm_at_shock(double thmfrac,double heat,double pspec,double cnorm,double enorm,double betat,
double ebreak,double emax,double emin,double einc,int nelec,double lelec[],double ledens[],double etemp[],
double dtemp[],double thmshock[],double plshock[],double elen[],double ntot,double eled[],double shden[],
double shelen[],double &bete){

    int i;
    
    double edtrm,pltrm,game,sum1,tmp1,renorm;
    
    //Make energy array between lelec[0] and log10(emax), add power-law and particle distribution together in 
    //appropriate ratios, then write into array and after shock, read it in and shift appropriately according
    //to shock.    
    //In theory this should be correctly normalized, but it's not. renorm to ntot at end is to make more exact
    gsl_interp_accel *acc_shock     = gsl_interp_accel_alloc();
    gsl_spline *spline_shock        = gsl_spline_alloc(gsl_interp_cspline, nelec);
    gsl_spline_init(spline_shock, lelec, ledens, nelec);
    
    for(i=0; i<nelec; i++){
        etemp[i]= pow(10, lelec[0]+(i+0.5)*einc);
        if(log10(etemp[i]) < lelec[nelec-1]){
            edtrm   = gsl_spline_eval(spline_shock, log10(etemp[i]), acc_shock);
            edtrm	= thmfrac*pow(10, edtrm);
        }
        else{
            edtrm	= 0;
        }
      
        thmshock[i] = edtrm; //thermal component assigned to an array at shock, so we can keep track of it
        if(etemp[i] >= emin){
            if(ebreak <= emin){//if break energy < min energy, power law distribution completely cooled
                pltrm	= cnorm*pow(etemp[i],-(pspec+1.)); 
              }
            else if(ebreak >= emax){ //if break energy > max energy, power law distribution completely uncooled
                pltrm = cnorm*pow(etemp[i],-pspec);
            }
            else{//if break energy within limits, set broken power law, power law index +1 after break energy
                if(etemp[i] < ebreak){
                    pltrm = cnorm*pow(etemp[i],-pspec);
                }
                else{
                    pltrm = cnorm*ebreak*pow(etemp[i],-(pspec+1.));
                }
            }
        }
        else{
            game = etemp[i]/emerg;
            if(game <= 1e3){
                bete	= sqrt(game*game-1.)/game;
            }
            else{
                bete	= 1. - 1./(2.*game*game) - 1.*(8*game*game*game*game);
            }
            pltrm	= enorm * game*game*bete * exp(-game*betat);
           	//this caused the squiggles, betat was a local variable when it should have been global!
        }
        /*if(etemp[i] < 2.*emin/heat && pltrm > edtrm){//avoid discontinuity between PL and thermal components
            pltrm	= edtrm;					 	  //may crash the code if heat is large
        }*/
        
        plshock[i] = pltrm;//power law component assigned to an array at shock, so can keep track of it
        dtemp[i]= log10(edtrm+pltrm);//total distribution sum of edtrm and pltrm (thermal + power law)        
        lelec[i]= log10(etemp[i]);
        elen[i]	= etemp[i];
    } // end of loop for nelec

    gsl_spline_free(spline_shock), gsl_interp_accel_free(acc_shock); // free the gsl_spline

    sum1    = 0.;
    
    for(i=0; i < nelec-1; i++){// integrate distribution, and then renormalise ntot
        tmp1	= (elen[i+1] - elen[i])*0.5*(pow(10, dtemp[i])+pow(10, dtemp[i+1])); 
           sum1	+= tmp1;
    }
    renorm	= ntot/sum1;        
    for(i=0; i<nelec; i++){
        ledens[i]= log10(renorm*pow(10,dtemp[i])); // probably don't need them here
        eled[i]	= elen[i]*pow(10,ledens[i]); // probably don't need them here
        shelen[i]= lelec[i];
        shden[i]= dtemp[i];
        plshock[i] = renorm*plshock[i]; //renormalise power law component at shock, to be read in after shock
        thmshock[i] = renorm*thmshock[i]; //renormalise thermal component at shock, to be read in after shock      
    }
}
