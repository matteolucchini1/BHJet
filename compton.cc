#include "agnjet.hh"

using namespace std;

/**
 * Main Compton function 
 */

void MICompton(bool isSSC,int Niter,int ncom,int njet,int nz,int nzdum,int k,double r,double ntot,double vol,
double dist,double dopfac[],gsl_spline *spline_eldis1, gsl_interp_accel *acc_eldis1, gsl_spline *spline_com2, 
gsl_interp_accel *acc_com2, double eltemp,double elenmn,double elenmx,double &ephmin,double &ephmax,
double ephxr[],double comspc[],double nucom[],double totcom[]){
    
    //Calling integration Compton function over electrons, returns #/cm**3/s/erg from #/cm**3/erg
    
    double yParR,yParNR,yPar,tau,com,blim,eph;
    int intYpar,m,i, N;
    
    double *comspc_src = new double[ncom]();
    double *nucom_src = new double[ncom]();
    double *comi = new double[ncom]();
    double *comspc_it = new double[njet*nzdum*ncom]();
    
    tau = r*ntot*sigtom;
    if (kboltz*eltemp/emerg < 1.0) {
        yParNR = pow(4.0*kboltz*eltemp/emerg,1)*max(tau*tau,tau);
        yPar = yParNR;
    } else {
        yParR = pow(4.0*kboltz*eltemp/emerg,2)*max(tau*tau,tau);
        yPar = yParR;
    }
    intYpar = int(yPar+0.5);

    //Cleaning up the array for new source term in between slices
    for(i=0; i<ncom; i++){
        comi[i] = 1e-100;
        totcom[i] = 1e-100;
    }

    //If this boolean variable is true, the function does single scattering
    if (isSSC){
        N = 1;
    } else {
        N = Niter;
    }

    for (int it=0; it<N;it++){
        //Cleaning stuff up before filling them
        for(i=0; i<ncom; i++){
            nucom_src[i] = -100.;
            comspc_src[i]= -100.;
        }
        
        //New source terms for higher orders
        if (it > 0) {
            for(i=0; i<ncom; i++){
                if (comi[i]==0) comi[i] = 1e-100;
                    nucom_src[i] = log10(ephxr[i]*herg);
                    comspc_src[i]= log10(comi[i])-log10(herg*herg*cee*ephxr[i]*pi*r*r);
                    
                	//The energy limits must be changed in order to avoid problems with interpolation 
                    ephmin = ephxr[0]*herg;
                    ephmax = ephxr[ncom-1]*herg;
            }
        } else {
            //This is a dummy filling for the it=0 because it's not needed
            for(i=0; i<ncom; i++){
                nucom_src[i] = log10(ephxr[i]*herg);
                comspc_src[i]= -100.;
            }
        }//end of if it
        gsl_interp_accel *acc_comiter      = gsl_interp_accel_alloc();
        gsl_spline *spline_comiter         = gsl_spline_alloc(gsl_interp_akima, ncom);
        gsl_spline_init(spline_comiter, nucom_src, comspc_src, ncom);
        
        for(i=0; i<ncom; i++){
            //Note: ephxr is in Hz
            eph     = ephxr[i]*herg;
            for(m=0; m<njet; m++){
                nucom[(m*nz+k)*ncom+i]= log10(eph*dopfac[m*nz+k]/herg);
            }
            //avoid out of bounds error
            if(i >= 2 && comspc[(0*nz+k)*ncom+(i-1)] <= -20. && 
            (comspc[(0*nz+k)*ncom+(i-1)]-comspc[(0*nz+k)*ncom+(i-2)]) <= 0. ){
                for(m=0; m<njet; m++){
                    comspc[(m*nz+k)*ncom+i]= comspc[(m*nz+k)*ncom+(i-1)]-4.;
                }
                continue;//Out of i(ncom) for loop for a new iteration
            }
            
            blim	= max(log(elenmn/emerg),log(eph/emerg))+1.e-15;
            
            if(blim >= log(elenmx/emerg)){
                com	= 1e-100;
            }
            else{
                if (it==0) {
                    com     = comintegral(blim,log(elenmx/emerg),eph,ephmin,ephmax,spline_eldis1,acc_eldis1,
                    		  spline_com2,acc_com2);
                } else {
                    com     = comintegral(blim,log(elenmx/emerg),eph,ephmin,ephmax,spline_eldis1,acc_eldis1,
                    		  spline_comiter,acc_comiter);
                }
            }
            comi[i] = com*ephxr[i]*herg*herg*vol;
            if (com == 0) comi[i] = 1e-100;
            //This is used for pair production rate estimation
            totcom[i] += com;
            
            // Converting from ergs/cm**2/s/Hz to mJy, writing in Hz and mJy           
            for(m=0; m<njet; m++){
                comspc[(m*nz+k)*ncom+i]=com*ephxr[i]*herg*herg*vol*dopfac[m*nz+k]*dopfac[m*nz+k]/
                						(mjy*4.*pi*dist*dist);                
                comspc_it[(m*nz+k)*ncom+i] = comspc_it[(m*nz+k)*ncom+i]+comspc[(m*nz+k)*ncom+i];      
                comspc[(m*nz+k)*ncom+i]=comspc_it[(m*nz+k)*ncom+i];                
                if(comspc[(m*nz+k)*ncom+i] == 0){
                    comspc[(m*nz+k)*ncom+i]=comspc[(m*nz+k)*ncom+(i-1)]-2.;
                }
                else{
                    comspc[(m*nz+k)*ncom+i]=log10(comspc[(m*nz+k)*ncom+i]);
                }                
            }            
        }//END for(ncom)
        
		if ((it == intYpar) || (it > 8)) {
			gsl_spline_free(spline_comiter), gsl_interp_accel_free(acc_comiter);
			break;
		}
		gsl_spline_free(spline_comiter), gsl_interp_accel_free(acc_comiter);
	} // END Niter
    delete[] comspc_src, delete[] nucom_src,  delete[] comi, delete[] comspc_it;
}

/**
 * This function is the kernel of eq 2.48 in Blumenthal & Gould(1970)
 */
double comfnc(double ein,void *p){
    struct comfnc_params *params = (struct comfnc_params *)p;
    double game     = (params -> game);
    double e1       = (params -> e1);
    gsl_spline *spline_com1 = (params -> spline_com1);
    gsl_interp_accel *acc_com1 = (params -> acc_com1);
    
    double com      = 0;
    double einit, utst, biggam, q, elg;
    double phonum	= 0;
    double tm1, tm2, tm3;
    double btst, eg4;
    
    einit	= exp(ein);
    btst	= einit/(game*emerg);
    eg4	= 4.*einit*game;
    utst	= eg4/(emerg+eg4);
    
    if((e1 < btst && (btst-e1)/btst >= 4.e-8) || (e1 > utst && (e1-utst)/utst >= 4.e-8)){
        com	= 0.;
        return com;
    }
    
    biggam	= eg4/emerg;
    q	= e1/(biggam*(1.-e1));
    elg	= ein*0.434294481903252; // log10(exp(1)) = 0.434294481903252
    phonum  = gsl_spline_eval(spline_com1, elg, acc_com1);

    tm1	= 2.*q*log(q);
    tm2	= (1.+2.*q)*(1.-q);
    tm3	= 0.5*((biggam*q)*(biggam*q))*(1.-q)/(1.+biggam*q);
    
    com	= (tm1+tm2+tm3)*pow(10.,phonum);
    
    return com;
}


/**
 * This function is for SSC flux -- numerically integrates over photon distribution and cross-section
 */
double comint(double gam,void *p){
    struct comint_params *params = (struct comint_params *)p;
    double eph      = (params -> eph);
    double ephmin   = (params -> ephmin);
    double ephmax   = (params -> ephmax);
    gsl_spline *spline_eldis1 = (params -> spline_eldis1);
    gsl_interp_accel *acc_eldis1 = (params -> acc_eldis1);
    gsl_spline *spline_com1 = (params -> spline_com1);
    gsl_interp_accel *acc_com1 = (params -> acc_com1);
    
    double game, econst, gnorm, blim, ulim, value, ellog, e1;
    double den	= 0;
    
    game	= exp(gam);
    e1	= eph/(game*emerg);
    
    econst	= 2.*pi*re0*re0*cee;
    gnorm	= emerg;

    blim	= max(log(eph/(4.*game*(game-eph/emerg))), log(ephmin))+1.e-14;
    ulim	= log(min(eph, ephmax));

    if(ulim <= blim){
        return 0;
    }
    
    // DON'T CHANGE THIS EPS!!!
    gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(1000);
    double result2, error2;
    gsl_function F2;
    struct comfnc_params F2params = {game, e1, spline_com1, acc_com1};
    F2.function     = &comfnc;
    F2.params       = &F2params;
    gsl_integration_qag(&F2, blim, ulim, 1e-1, 1e-1, 1000, 1, w2, &result2, &error2);
    value   = result2;
    gsl_integration_workspace_free(w2);
    
    ellog	= log10(game*emerg);
    den     = gsl_spline_eval(spline_eldis1, ellog, acc_eldis1);
    den     = pow(10., den);
    
    return gnorm*econst*den*value/game;
}

double comintegral(double blim,double ulim,double eph,double ephmin,double ephmax,gsl_spline *spline_eldis1,
gsl_interp_accel *acc_eldis1,gsl_spline *spline_com1,gsl_interp_accel *acc_com1){
    gsl_integration_workspace *w    = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F1;
    struct comint_params F1params   = {eph, ephmin, ephmax, spline_eldis1, acc_eldis1, spline_com1, acc_com1};
    F1.function     = &comint;
    F1.params       = &F1params;
    gsl_integration_qag(&F1, blim, ulim, 1e-1, 1e-1, 1000, 1, w, &result, &error);
    gsl_integration_workspace_free(w);
    
    return result;
}

void ucompton(int nz,int njet,int nsyn,int k,double dist,double r,double nusyn[],double dopfac[],
double synabs[],double &ucom){

    double *phofrq	= new double[nsyn]();
    double *phoint	= new double[nsyn]();
    double sum1,tmp1;

    int i,m;

    sum1    = 0.;
    for(m=0; m < njet; m++){
        for(i=0; i < nsyn; i++){
            phofrq[i]= pow(10.,nusyn[(m*nz+(k-1))*nsyn+i])/dopfac[m*nz +(k-1)];
            //3 powers of dopfac because one was inside synint for angle aberration            
            phoint[i]= mjy * pow(10., synabs[(m*nz+(k-1))*nsyn+i])/pow(dopfac[m*nz +(k-1)],3);
            //Handle differently: not flux but assuming that the synchrotron was produced in this segment. 
            //That's overestimating but hopefully not too much because r, delz also increased.
            phoint[i]*= 4. * dist*dist/(r*r*cee);
		}
        //Rough integration to get approx energy density
        for(i=0; i< nsyn-1; i++){
            tmp1	= (phofrq[i+1]-phofrq[i])*0.5*(phoint[i]+phoint[i+1]);
            sum1	+= tmp1;
        }
    }
    ucom    = sum1;

    delete[] phoint, delete[] phofrq;
}

void seed_sync_and_disk_phtns(bool isVerbose,int disksw,int bbsw,double tin,double rin,
double rout,double gamv,int nsyn,double snumax, double z,double r,double reff2,double hbb,double tbbeff,
double nphot[],double nurad[],double nubb[],double energ[], double phodis[],double ephot[],
double &ephmax,double &ephmin){

    double *jetu	= new double[nsyn]();
    double *disku	= new double[nsyn]();
    
    double bbnumax,bbdil,photmx,bbcon,bbirad;
    double peaksw;
    int i;
    
    //spline_dbl photon density (erg/s/Hz) is passed to Compton routines, we need to divide by
    //herg*cee*energy(erg)*pi*r**2 (erg**2*s**2*cm**3*Hz*s**-1) to get #/cm**3/erg
    //Also adding in the emission from the disk, making the max energy based on tin
    for(i=0; i<nsyn; i++){
        if(nphot[i] == 0){
            nphot[i]= nphot[i-1]-4.;
        }
        else{
            nphot[i]= log10(nphot[i]);
        }
        jetu[i]	= pow(10,nphot[i])/(cee*herg*energ[i]*pi*r*r);
    }

    if(isVerbose){
        //Checking SSC energy density to compare to ucom1 at shock
        double sum1	= 0, tmp1 = 0,toten;
        for(i=0; i<nsyn-1; i++){
            tmp1	= (pow(10, nubb[i+1])-pow(10,nubb[i]))*0.5*(pow(10,nphot[i+1])+pow(10,nphot[i]));
            sum1	+=tmp1;
        }
        toten	= sum1/(cee*pi*r*r);
        cout << "ENDEN FOR SSC: " << toten << endl;
    }
    
    if(disksw==1){
        bbnumax	= log10(1e2*tin*kboltz/herg);
        if(bbnumax > snumax){
            cerr << "Need to change BB+syn for IC" << endl;
            cerr << "BB going to higher freq than syn" << endl;
            exit(1);
        }
        bbdil	= 2.*reff2/(2.*reff2+(z-hbb/2.)*(z-hbb/2.))-reff2/(reff2+(z-hbb/2.)*(z-hbb/2.));
        
        //Stopping Compton loop when too far past photon peak
        photmx	= 1.e-20;
        peaksw	= 0;
        bbcon	= 0;
        bbirad = 0;
        
        for(i=0; i<nsyn; i++){
            if(bbsw==1){
                if (peaksw!=1){
                    if(nubb[i] <= bbnumax){
                        if(photmx > 1.e-20 && bbcon/photmx <= 1.e-13){
                            if(isVerbose){
                                cout << left << setw(15) << "gone too low" << setw(15) << i << setw(15) 
                                	 << photmx << setw(15) << bbcon << endl;
                            }
                            bbcon	= 0.;
                            peaksw	= 1.;
                        }
                        else{
							bbcon	= bbdisk(log(nurad[i]), gamv, tin, rin, rout, z, hbb)/
									  (cee*herg*energ[i]*nurad[i]);                               
                            photmx	= max(photmx,bbcon);
                        }//end else(photmx)
                    }//END if(nubb[i] <= bbnumax)
                    else{
                        bbcon   = 0.;
                    }//END nubb
                    bbirad	= bbdil*2.*herg*pow(nurad[i],3)/(cee*cee*(exp(energ[i]/(kboltz*tbbeff))-1.)
                    		  *cee*herg*energ[i]);
                }//END if(peak!=1)
                
				if(z > hbb/2){
					phodis[i]= log10(jetu[i]+bbcon+bbirad);
					if (jetu[i]+bbcon == 0) phodis[i]= -500;
						disku[i]= bbcon + bbirad;
					} else{
						phodis[i]= log10(jetu[i]+bbcon);
					if (jetu[i]+bbcon == 0) phodis[i]= -500;
					disku[i]= bbcon;
				}
			} else{ //END if(bbsw==1)
				phodis[i]= log10(jetu[i]);
				if (jetu[i] == 0) phodis[i]= -500;
			}//end bbsw
		}//END for-loop nsyn
	}//END if(disksw==1)
    else{
        for(i=0; i<nsyn; i++){
            phodis[i]= log10(jetu[i]);
            if (jetu[i] == 0) phodis[i]= -500;
        }
    }// disksw!=1
    ephmax	= pow(10,ephot[nsyn-1]);
    ephmin	= pow(10,ephot[0]);
    
    delete[] disku,delete[] jetu;
}
