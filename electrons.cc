/*
 * Here the electron distribution is calculated based on the input parameters.
 * ----- INPUT ------
 * mxsw:    percentage of thermal (mxsw = 1) vs non thermal electrons (mxsw = 0)
 * pspec:   power law index
 * gamfac:  if mxsw !=1 is used to calculate PL
 * endnsmj: some average number from MJ
 * ntot0:   initial number density of electrons
 * cnorm:   normalization from equipartition
 * eltemp:  electron temperature from init params
 * ----- OUTPUT -----
 * rdlgen:  energy array in logspace
 * rdedn:   energy density in logspace
 * thmbase: thermal distribution at the base
 * plbase:  nonthermal distribution at the base
 */

#include "agnjet.hh"

using namespace std;

void electrons(int mxsw,int nelec,double pspec,double gamfac,double endnsmj,double ntot0,double cnorm,
double eltemp,double &gamax0,double &endens,double rdlgen[],double rdedn[],double thmbase[],double plbase[],
double &bete){
    
    double betat,enorm,elmin,elmax,einc,game,reled,emin,emax,gamin,edtrm,pltrm;
    double *rden = new double[nelec]();    
    int i;
    
    //Gamax0 is set to unity so that gamax/gamax0 gives the correct Mach number factor for energy shift
    gamax0  = 1.;
    emin    = 2.23*kboltz*eltemp;
    gamin   = emin/emerg;
    emax    = gamfac*gamin*emerg;

	//generating maxwellian/power law distirbutions here based on mxsw = 1 or mxsw = 0

    if(mxsw == 1) {
        endens  = endnsmj*ntot0;
        betat   = emerg/(kboltz*eltemp);
        enorm   = ntot0*betat/(emerg*k2_fnc(betat));
        elmin	= emerg;
        elmax	= 50.*kboltz*eltemp;
        einc	= log10(elmax/elmin)/nelec;
        
        for(i=0; i<nelec; i++){
            rdlgen[i] = log10(elmin)+(i+0.5)*einc;
            rden[i]	= pow(10.,rdlgen[i]);
            game	= rden[i]/emerg;
            if(game < 1.e-3){
                bete	= sqrt(game*game - 1.)/game;
            }
            else{
                bete	= 1.-1./(2.*game*game)-1./(8.*game*game*game*game);
            }
            reled	= enorm*game*game*bete*exp(-game*betat);
            rdedn[i]= log10(reled);
        }
    }
    else if(mxsw == 0) {    
        if(pspec >= 1.99 && pspec <= 2.01){
            endens	= cnorm*log(emax/emin);        }
        else{
            endens	= cnorm*(pow(emax,2.-pspec) - pow(emin,2.-pspec))/(2.-pspec);
        }
        einc	= log10(emax/emin)/nelec;
        for(i=0; i<nelec; i++){
            rdlgen[i] = log10(emin) + (i+0.5)*einc;
            rden[i]	= pow(10,rdlgen[i]);
            reled	= cnorm*pow(rden[i],-pspec);
            rdedn[i]= log10(reled);
        }
    }

    //This case represents a mixed thermal/non-thermal injected particle distribution (i.e. mxsw != 0 && 
    //mxsw != 1). We start by setting up a thermal distribution as we would when mxsw = 1, with the energy 
    //array extending from elmin to elmax, where elmin is game = 1, elmax is 50kT. We assign these electrons 
    //to variable edtrm. The energy array is then interpolated to fit the new grid, elmin - emax (emax is set 
    //by gamfac, it is the maximum energy of the power-law distirbution), using the gsl spline. Then we 
    //calculate the power-law distribution, and append both edtrm and pltrm to the full distribution, as well 
    //as setting up the arrays thmcomp and plcomp for future use.

    else {
        
        endens  = endnsmj*ntot0;
        betat   = emerg/(kboltz*eltemp);
        enorm   = mxsw * ntot0 * betat/(emerg*k2_fnc(betat));
        elmin	= emerg;
        elmax	= 50.*kboltz*eltemp;
        einc	= log10(emax/elmin)/nelec;
        
        for(i=0; i<nelec; i++){
            rdlgen[i] = log10(elmin)+(i+0.5)*einc;
            rden[i]	= pow(10.,rdlgen[i]);
            if(rden[i] <= elmax) {
                game	= rden[i]/emerg;
                if(game < 1.e-3){
                    bete = sqrt(game*game - 1.)/game;
                }
                else{
                    bete = 1.-1./(2.*game*game)-1./(8.*game*game*game*game);
                }
                edtrm = enorm*game*game*bete*exp(-game*betat);                
            }
            else {
                edtrm = 0.;
            }
            
            if(rden[i] >= emin) {
                pltrm = cnorm*pow(rden[i],-pspec);
            }
            else{
                pltrm = 0.;
            }
            
            thmbase[i] = edtrm;
            plbase[i] = pltrm;            
            reled = edtrm + pltrm;
            rdedn[i] = log10(reled);            
        }        
    }
    delete[] rden;
}

void pl_and_th_comp(bool isShock,double thmfrac,int &nw,int nelec,double mxsw,double rdlgen[],
double thmbase[],double ebreak,double emin,double emax,double gshift,double gshock,double mjteff,
double ratio_ne,double ntot, double cnorm,double enorm,double pspec,double &pltrm,double etemp[],
double dtemp[],double plcomp[],double thmcomp[],double &bete){
    
    int i;
    
    double game,sum1,tmp1,renorm,betat;
    double *lelec	= new double[nelec]();
    double *elen	= new double[nelec]();
    
    betat   = emerg/(kboltz*mjteff);
    
    for(i=0; i<nelec; i++){
        
        if (!isShock) etemp[i] = pow(10,rdlgen[i]);
        
        if(etemp[i] >= emin){
            if(ebreak <= emin){
                pltrm	= cnorm*pow(etemp[i],-(pspec+1.));//if break energy<min energy, PL completely cooled
                
            }
            else if(ebreak >= emax){
                pltrm = cnorm*pow(etemp[i],-pspec); //if break energy>emax, PL completely uncooled
            }
            else{
                if(etemp[i] < ebreak){
                    pltrm = cnorm*pow(etemp[i],-pspec);	//if break energy inside energy limits, broken power 
                    								   	//law distribution with break at ebreak
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
        }
        plcomp[i] = log10(pltrm);
    }
    
    nw = 0;
    
    if (!isShock){
        for(i=0; i<nelec; i++) {
            lelec[i] = rdlgen[i];
            elen[i] = pow(10.,lelec[i]);
            etemp[i] = log10(gshift * elen[i]);
            if(gshift * elen[i]/emerg <= 1){
                nw = nw + 1;
            }
            if(thmbase[i] > 0){
                thmcomp[i] = log10((ratio_ne)*thmbase[i]);
            }
            else{
                thmcomp[i] = 0.;
            }
            
        }
        
        sum1 = 0.;
        for(i=0; i<nelec-1; i++) {
            if(thmcomp[i] == 0. || thmcomp[i+1] == 0.){
                tmp1 = 0;
            }
            else{
                tmp1 = (pow(10.,etemp[i+1])-pow(10.,etemp[i]))*0.5*(pow(10.,thmcomp[i])+pow(10.,thmcomp[i+1]));
            }
            sum1 += tmp1;
        }
        
        renorm = log10((ntot*mxsw)/sum1);        
        
        for(i=0; i<nelec; i++) {
            if(thmcomp[i] == 0.){
                thmcomp[i] = 0.;
            }
            else{
                thmcomp[i] += renorm;
            }
            
            if(etemp[i] > emin && plcomp[i] > (enorm * game*game*bete * exp(-game*betat))){
                plcomp[i] = thmcomp[i];
            }
            if(thmcomp[i] == 0 || plcomp[i] == 0){
                if(thmcomp[i] == 0.){
                    dtemp[i] = log10(pow(10.,plcomp[i]));
                }
                else{
                    dtemp[i] = log10(pow(10.,thmcomp[i]));
                }
            }
            
            else{
                dtemp[i] = log10(pow(10.,thmcomp[i]) + pow(10.,plcomp[i]));
            }
        }
    } else{
		//Now shift down the thermal component (thmcomp), determine the power law component, sum them, and 
		//renormalise everything accordingly  - this differs from previous version, in that before we were  
		//shifting down the full distribution, wheras now we shift down the thermal component only, and 
		//recalculate the PL component according to the thermal peak of the adiabatically cooled thermal 
		//distribution. In this way  the break can be re-calculated at each step in a simpler way. 
        
        nw	= 0;
        for(i=0; i<nelec; i++){
            etemp[i] = log10(etemp[i]);
            etemp[i]= log10(gshift/gshock*pow(10,etemp[i]));
            if(pow(10,etemp[i])/emerg <= 1.){
                nw++;
            }
        }
        
        //calculating newly shifted down thermal component
        sum1	= 0;
        for(i=0; i<nelec-1; i++){
            tmp1	= (pow(10,etemp[i+1]) - pow(10,etemp[i]))*0.5*(pow(10,thmcomp[i]) + pow(10,thmcomp[i+1]));
            sum1	+= tmp1;
        }
        renorm	= log10(ntot*thmfrac/sum1);
        for(i=0; i<nelec; i++){
            if(thmcomp[i]>0){ //only read in the thermal component (and take log value) if it is non-zero.
                thmcomp[i]+= renorm;
            }
            else{
                thmcomp[i]=0;
            }
            if(pow(10.,etemp[i])*(gshock/gshift) < 2.*emin && plcomp[i]>thmcomp[i]){//smooth out distirbution
                plcomp[i]=thmcomp[i];												//if PLcomponent > thermal 
            }																		
            if(thmcomp[i]>0){ //again, thermal component should only be included in total sum if it's non-zero
                dtemp[i] = log10(pow(10.,thmcomp[i])+pow(10.,plcomp[i]));
            }
            else{
                dtemp[i] = log10(pow(10.,plcomp[i]));
            }
        }
    }// end isShock
    
    delete[] lelec, delete[] elen;
}

void update_ebreak(int nz,int k,double z,double r, double h0,double brk,double bfield,double beta,
double &ebreak,double &gebreak){
    
    
    double delebreak,syncon,ub,tdyn;
    
    //New electron energy/dist arrays for printing to output file
    double *tdyntot = new double[nz]();
    double *tsyntot = new double[nz]();
    double *ebrarr = new double[nz]();

    ub = bfield * bfield/(8. * pi);
    syncon	= 4./3. * sigtom * ub/(emgm*emgm*cee*cee*cee);
    if(z < h0) {
        tdyn = brk*h0/(beta*cee);
        tdyntot[k] = tdyn;
        tsyntot[k] = (1./syncon);
        ebrarr[k] = tsyntot[k]/tdyntot[k];
        ebreak = ebrarr[k];
    }
    else {
        tdyn    = brk*r/(beta*cee);
        tdyntot[k] = tdyn + tdyntot[k-1];
        tsyntot[k] = (1./syncon) + tsyntot[k-1] ;
        delebreak = tsyntot[k]/tdyntot[k];        
        ebrarr[k] += delebreak;
        ebreak = ebrarr[k];
    }
    gebreak = ebreak/emerg;
    delete[] tdyntot, delete[] tsyntot, delete[] ebrarr;
}

//Reads electron distribution and renormalizes it
void read_and_shift_eldists(int nelec,int &nw,double ntot,double ntot0,double gshift,double rdlgen[],
double rdedn[],double elen[],double ledens[],double etemp[],double dtemp[]){
    
    int i;
    double sum1,tmp1,renorm;
    
    for(i=0; i<nelec; i++){
        elen[i]	= pow(10.,rdlgen[i]);
        ledens[i]= log10(ntot/ntot0*pow(10.,rdedn[i]));
    }
    for(i=0; i< nelec; i++){
        etemp[i]= log10(gshift*elen[i]);
        if(gshift*elen[i]/emerg <= 1.){
            nw = nw+1;
        }
        dtemp[i]= ledens[i];
    }
    sum1    = 0.;
    for(i=0; i < nelec-1; i++){
        tmp1	= (pow(10.,etemp[i+1])-pow(10.,etemp[i])) * 0.5 * (pow(10.,dtemp[i])+pow(10.,dtemp[i+1]));
        sum1	+= tmp1;
    }
    renorm	= log10(ntot/sum1);
    for(i=0; i<nelec; i++){
        dtemp[i] += renorm;
    }
}

void lost_ele(int &nw,int nelec,double etemp[],double dtemp[],double ntot,double einc,double elen[],
double eled[],double lelec[],double ledens[]){
    
    int i;
    double totlos =0;
    double tmlos,oldnum,renorm,sum1,tmp1;
    
    if(nw>0){
        if(nw==1){
            totlos  = (pow(10.,etemp[1])-pow(10.,etemp[0]))*pow(10.,dtemp[0]);
        }
        else if(nw > 1){        
            totlos  = 0;
            for(i=0; i < nw-1; i++){
                tmlos	= (pow(10.,etemp[i+1])-pow(10.,etemp[i]))*0.5*(pow(10.,dtemp[i])+pow(10.,dtemp[i+1]));
                totlos	+= tmlos;
            }
        }
        oldnum  = ntot;
        ntot    -= totlos;
        renorm  = log10(ntot/oldnum);
        for(i=0; i<nelec; i++){
            dtemp[i] += renorm;
        }
        
        gsl_interp_accel *acc_shift  = gsl_interp_accel_alloc();
        gsl_spline *spline_shift     = gsl_spline_alloc(gsl_interp_cspline, nelec);
        gsl_spline_init(spline_shift, etemp, dtemp, nelec);
        
        einc    = (etemp[nelec-1]-etemp[nw])/nelec;
        for(i=0; i<nelec; i++){
            lelec[i]= etemp[nw]+(i+0.5)*einc;
            ledens[i] = gsl_spline_eval(spline_shift, lelec[i], acc_shift);
            elen[i] = pow(10., lelec[i]);
            eled[i] = elen[i]*pow(10., dtemp[i]);
        }
        gsl_spline_free(spline_shift), gsl_interp_accel_free(acc_shift);
    }//END nw > 0    
    else{
        for(i=0; i<nelec; i++){
            lelec[i]= etemp[i];
            elen[i] = pow(10.,lelec[i]);
            eled[i] = elen[i]*pow(10.,dtemp[i]);
            ledens[i]= dtemp[i];            
        }
    }//END else nw >0
    
    sum1 = 0.;
    for(i=0; i<nelec-1; i++) {
        tmp1 = (pow(10.,lelec[i+1]) - pow(10.,lelec[i])) * 0.5 * (pow(10.,ledens[i]) + pow(10.,ledens[i+1]));
        sum1 += tmp1;
    }
    renorm	= log10(ntot/sum1);
    for(i=0; i<nelec; i++){
        ledens[i] += renorm;
    }
}
