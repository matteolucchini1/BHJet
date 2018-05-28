
/**
 * equipartition function
 * ***********************
 * @param mxsw          maxwellian switch [jet parameter]
 * @param equip         equipartition factor [jet parameter]
 * @param pspec         power-law index [jet parameter]
 *
 * @return cnorm        normalization factor of the power-law distribution <!-- FIXME -->
 * @return ntot0        initial number of electron
 * @return b_en         magnetic energy in nozzle
 * @return b0         	magnetic field in nozzle
 * @return gamax0       normalized Lorentz factor for shifting electrons
 *
 */

#include "agnjet.hh"

using namespace std;

void equipartition(double velsw,double mxsw,int nenpsw,double eta,double equip,double pspec,double nprot0,
double emin,double emax,double endnsmj,double &cnorm,double &ntot0,double &b_en,double &b0){
	
	cnorm   = 0;

    if(nenpsw==0){
       	if(mxsw==1){
			if(equip==1){
				ntot0 = nprot0;
				b_en = ntot0*endnsmj;
			}
			else if(equip==0){
				ntot0 = nprot0;
				b_en = ntot0*(pmgm*cee*cee-endnsmj);
			} else{
				ntot0 = nprot0;
				b_en = equip*ntot0*endnsmj;
			}
		} else{
			if(equip==1){
				ntot0   = nprot0;
				cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
				if(pspec >= 1.99 && pspec <= 2.01){
					b_en = cnorm*log(emax/emin);
				} else{
					b_en = cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
				}
			}else if(equip==0){
            	ntot0 = nprot0;
                cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
               	if(pspec >= 1.99 && pspec <= 2.01){
					b_en = ntot0*pmgm*cee*cee - cnorm*log(emax/emin);
                } else{
					b_en = ntot0*pmgm*cee*cee-cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
                }
			} else{
				ntot0 = nprot0;
            	cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
            	if(pspec >= 1.99 && pspec <= 2.01){
					b_en = equip*cnorm*log(emax/emin);
				} else{
					b_en = equip*cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
            	}
			}
        }
	} else{
		if (velsw <= 1){
			if(mxsw==1){
				ntot0 = nprot0*pmgm*cee*cee/((1.+equip)*endnsmj);
				b_en = equip*ntot0*endnsmj;
			}else if(mxsw == 0){
				if(pspec >= 1.99 && pspec <= 2.01){
					cnorm = nprot0*pmgm*cee*cee/((1.+equip)*log(emax/emin));
					b_en = equip*cnorm*log(emax/emin);
				} else{
					cnorm = nprot0*pmgm*cee*cee*(2.-pspec)/((1.+equip)*(pow(emax,(2.-pspec))-
							pow(emin,(2.-pspec))));
					b_en = equip*cnorm*(pow(emax,(2.-pspec))-pow(emin,(2.-pspec)))/(2.-pspec);
				}
				ntot0 = cnorm*(pow(emax,(1.-pspec))-pow(emin,(1.-pspec)))/(1.-pspec);
			} else {
				ntot0 = nprot0*pmgm*cee*cee/((1.+equip)*endnsmj);
           		b_en = equip * ntot0 * endnsmj;
           		cnorm = (ntot0*endnsmj*(1.-mxsw)*(2.-pspec))/(pow(emax,(2.-pspec))-pow(emin,(2.-pspec)));
			}
		} else { 
			ntot0 = nprot0*eta;
			if(mxsw==1){
				b_en = equip*ntot0*endnsmj;
			} else if (mxsw==0){
				cnorm = ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
				if(pspec >= 1.99 && pspec <= 2.01){				
					b_en = equip*cnorm*log(emax/emin);
				} else {
					b_en = equip*cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec);
				}
			} else {
				cnorm = (1.-mxsw)*ntot0*(1.-pspec)/(pow(emax,1.-pspec)-pow(emin,1.-pspec));
				if(pspec >= 1.99 && pspec <= 2.01){				
					b_en = equip*(mxsw*ntot0*endnsmj+cnorm*log(emax/emin));
				} else {
					b_en = equip*(mxsw*ntot0*endnsmj+cnorm*(pow(emax,2.-pspec)-pow(emin,2.-pspec))/(2.-pspec));
				}
			}
		}
	}
		
	b0 = sqrt(8.*pi*b_en);
}
