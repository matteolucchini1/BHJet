#include "bhjet.hh"

//This function determines very, very roughly whether the Compton emission from a zone is worth computing or 
//not. The criteria are a) are we in the first zone (which we almost always care about because it's the 
//corona) or b) do we expect the non-thermal SSC luminosity to be sufficiently bright. Note that in XRBs for
// standard parameters this function returns false in most zones; also, for LLAGN (velsw not higher than 1, 
//large BH mass) it assumes you are not trying to compute the gamma ray spectrum. This is because neither
//class of objects actually has any gamma ray detections
bool Compton_check(bool IsShock,int i,double Mbh,double Urad,double velsw,zone_pars &zone){
    double Lumnorm,Ub,Usyn,Lsyn,Lcom;
    Lumnorm = pi*pow(zone.r,2.)*zone.delz*pow(zone.delta,4.)*zone.lepdens*sigtom*cee*zone.avgammasq;
    Ub = pow(zone.bfield,2.)/(8.*pi); 
    Lsyn = Lumnorm*Ub;
    Usyn = Lsyn/(pi*pow(zone.r,2.)*cee*pow(zone.delta,4.));
    Lcom = Lumnorm*(Usyn+Urad);

    bool test1 = (Lcom/Lsyn > 2e-2);
    bool test2 = (Lcom >= 1.26e38*Mbh*1e-7);

    //the logic for these tests is as follows: 
    //1) always include the jet nozzle/corona, and the first region after it just in case
    //2) if large BH mass and not using blhet, assume we're dealing with a LLAGN, hence no gamma ray data
    //available, hence no need to calculate the non-thermal IC spectrum
    //3) in any other zone, make sure that a) the IC emission from the zone is at least 0.02 times that of them
    //synchrotron component and b) make sure this luminosity is not lower than 1e-7 the Eddington luminosity
    //of the black hole. 1e-7 is taken to be an arbitrarily small number, anything lower will just not be 
    //detectable anyway. 

    if (i <= 1){
        return true;
    }
    else if (Mbh > 1.e4 && velsw <= 1) {
        return false;
    } else if(test1 == true && test2 == true && IsShock == true){
        return true;
    }
}

void param_write(const double *par,std::string path){
    std::ofstream file;
    file.open(path.c_str(),std::ios::trunc);	
    
    for (int k=0;k<27;k++) {
        file << par[k] << std::endl;
    }
    
    file.close();
}

//Plots a given array in units of ergs (x axis) and erg/s/Hz (y axis) to the file on the provided path. The 
//overloard with or without const is to be able to pass arrays directly from the radiation libraries
//note: the factor 1+z in the specific luminosity calculation is to ensure that the output spectrum only moves
//to lower frequency, not up/down. 
void plot_write(int size,double *en,double *lum,std::string path,double dist,double redsh){
    std::ofstream file;
    file.open(path.c_str(),std::ios::app);	

    for(int k=0;k<size;k++){
        file << en[k]/(herg*(1.+redsh)) << " " << lum[k]*(1.+redsh)/(4.*pi*pow(dist,2.)*mjy) << std::endl;
    }

    file.close();	
}

void plot_write(int size,const double *en,const double *lum,std::string path,double dist,double redsh){
    std::ofstream file;
    file.open(path.c_str(),std::ios::app);	

    for(int k=0;k<size;k++){
        file << en[k]/(herg*(1.+redsh)) << " " << lum[k]*(1.+redsh)/(4.*pi*pow(dist,2.)*mjy) << std::endl;
    }

    file.close();	
}

//Same as above but for particle distributions
void plot_write(int size,const double *p,const double *g,const double *pdens,const double *gdens,
				std::string path){
			
    std::ofstream file;
    file.open(path.c_str(),std::ios::app);
    for(int k=0;k<size;k++){
     file << p[k] << " " << g[k] << " " << pdens[k] << " " << gdens[k] << std::endl;
    }	

    file.close();
}

//This function takes the observed arrays of the Cyclosyn and Compton classes for jet/counterjet
//sums up the contributions of both and stores them in one array of observed frequencies and one
//of comoving luminosities
void sum_counterjet(int size,const double* input_en,const double* input_lum,double* en,double* lum){
    double en_cj_min,en_j_min,en_cj_max,en_j_max,einc;
    double *en_j = new double [size];
    double *en_cj = new double [size];
    double *lum_j = new double[size];
    double *lum_cj = new double[size];	

    en_j_min = input_en[0];
    en_cj_min = input_en[size];
    en_j_max = input_en[size-1];
    en_cj_max = input_en[2*size-1];	
    einc = (log10(en_j_max)-log10(en_cj_min))/(size-1);

    for(int i=0;i<size;i++){
        en[i] = pow(10.,log10(en_cj_min)+i*einc);
        en_j[i] = input_en[i];
        en_cj[i] = input_en[i+size];
        lum_j[i] = input_lum[i];
        lum_cj[i] = input_lum[i+size];
    }	

    gsl_interp_accel *acc_j = gsl_interp_accel_alloc();
    gsl_spline *spline_j = gsl_spline_alloc(gsl_interp_akima,size);
    gsl_spline_init(spline_j,en_j,lum_j,size);

    gsl_interp_accel *acc_cj = gsl_interp_accel_alloc();
    gsl_spline *spline_cj = gsl_spline_alloc(gsl_interp_akima,size);
    gsl_spline_init(spline_cj,en_cj,lum_cj,size);

    for(int i=0;i<size;i++){
        if(i==0){
            lum[i] = lum_cj[i];
        } else if (i==size-1){
            lum[i] = lum_j[i];
        }else if(en[i] < en_j_min){
            lum[i] = gsl_spline_eval(spline_cj,en[i],acc_cj);
        } else if (en[i] < en_cj_max){
            lum[i] = gsl_spline_eval(spline_j,en[i],acc_j) + gsl_spline_eval(spline_cj,en[i],acc_cj); 
        } else {
            lum[i] = gsl_spline_eval(spline_j,en[i],acc_j);//note: the factor 0.999 is to avoid 
														         	//occasional gsl interpolation errors 
														         	//due to numerical inaccuracies
        }
    }

    gsl_spline_free(spline_j), gsl_interp_accel_free(acc_j);
    gsl_spline_free(spline_cj), gsl_interp_accel_free(acc_cj);

    delete[]en_j,delete[]en_cj,delete[]lum_j,delete[]lum_cj;

    return;
}

//Calculates the redshifted spectrum as seen by the observer, starting from the emitted spectrum in the frame
//comoving with the source. Only applicable to (distant) AGN, not to galactic XRBs.
void output_spectrum(int size,double* en,double* lum,double* spec,double redsh,double dist){
	
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_akima,size);
    gsl_spline_init(input_spline,en,lum,size);	

    for(int k=0;k<size;k++){
        if (en[k]*(1.+redsh) < en[size-1]){
           spec[k] = log10(gsl_spline_eval(input_spline,en[k]*(1.+redsh),acc)*(1.+redsh)/(4.*pi*pow(dist,2.)*mjy));
        }
        else {
           spec[k] = -50.;
        }	
    }	
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

//Used for summing individual zone contributions for a generic spectral component from code: pre/post particle 
//acceleration synchrotron, pre/post particle acceleration Comptonization
//The second function does the same, but sums the disk/corona/bb to the total jet spectrum. The reason for the
//const arryas in input is that the input arrays are directly accessed from the ShSDisk class, which are const
void sum_zones(int size_in,int size_out,double* input_en,double* input_lum,double* en,double* lum){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_akima,size_in);
    gsl_spline_init(input_spline,input_en,input_lum,size_in);	

    for(int i=0;i<size_out;i++){
        if (en[i] > input_en[0] && en[i] < input_en[size_in-1]){
    	    lum[i] = lum[i] + gsl_spline_eval(input_spline,en[i],acc);
        }
    }
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

void sum_ext(int size_in,int size_out,const double* input_en,const double* input_lum,double* en,double* lum){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_akima,size_in);
    gsl_spline_init(input_spline,input_en,input_lum,size_in);	

    for(int i=0;i<size_out;i++){
       if (en[i] > input_en[0] && en[i] < input_en[size_in-1]){
            lum[i] = lum[i] + gsl_spline_eval(input_spline,en[i],acc);
       }
    }
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

//Simple numerical integral to calculate the luminosity between numin and numax of a given array; input units
//are erg for the frequency/energy array, erg/s/Hz for the luminosity array to be integrated, Hz for the 
//integration bounds; size is the dimension of the input arrays
//Note: this uses a VERY rough method and wide bins, so thread carefully
double integrate_lum(int size,double numin,double numax,const double* input_en,const double* input_lum){
    double temp = 0.;
    for (int i=0;i<size-1;i++){
        if (input_en[i]/herg > numin && input_en[i+1]/herg < numax) {
            temp = temp+(1./2.)*(input_en[i+1]/herg-input_en[i]/herg)*(input_lum[i+1]+input_lum[i]);
        }		
    }
    return temp;
}

//Overly simplified estimate of the photon index between numin and numax of a given array; input is the same
//as integrate_lum. Note that this assumes input_lum is a power-law in shape
double photon_index(int size,double numin,double numax,const double* input_en,const double* input_lum){
    int counter_1, counter_2 = 0;	
    double delta_y, delta_x, gamma;	
    for (int i=0;i<size;i++){
        if (input_en[i]/herg < numin){
            counter_1 = i;
        }
        if (input_en[i]/herg < numax){
    	    counter_2 = i;
        }
    }	
    delta_y = log10(input_lum[counter_2])-log10(input_lum[counter_1]);
    delta_x = log10(input_en[counter_2]/herg)-log10(input_en[counter_1]/herg);
    gamma = delta_y/delta_x - 1.;
    return gamma;
}

//Prepares files for above printing functions at the start of the run. There are two reasons this exists:
//1) specifying the units of the output obviously makes things easier to read and
//2) S-lang does not allow to pass ofstream objects to functions, so we need to pass a path and open the file
//from the path inside the write function. As a result, it's impossible to just truncate and clean the files
//at the start of each iteration. This is only relevant for the cyclosyn_zones, compton_zones, and numdens
//files
void clean_file(std::string path,int check){
    std::ofstream file;
    file.open(path.c_str(),std::ios::trunc);

    if (check == 1) {
        file << std::left << std::setw(20) << "#Parameter:" << std::endl;
    } else if (check==2){
        file << std::left << std::setw(20) << "#nu [Hz] " << std::setw(20) << "Flux [mJy] " << std::endl;
    } else if (check==4){
        file << std::left << std::setw(20) << "#p [g cm s-1] " << std::setw(20) << "g [] " << std::setw(20) << 
             " n(p) [# cm^-3 p^-1] " << std::setw(20) << " n(g) [# cm^-3 g^-1]" << std::setw(20)<<  std::endl;
    } else {
        std::cout << "File to be cleaned not supported" << std::endl;
    }	

    file.close();
}

//Used for interpolation by slang code
void jetinterp(double *ear,double *energ,double *phot,double *photar,int ne,int newne){
    int i, iplus1, j, jplus1;
    double emid, phflux;

    j = 0; 	
    for(i=0; i<newne; i++){
        // Middle of bin
        iplus1 = i+1;
        emid = (ear[i] + ear[iplus1])/2.;
        // Linear search if we don't bracket yet
        if(j == -1){
            j = 1;
        }
        while(j <= ne && energ[j] < emid){
            j++;
        }

        jplus1	= j;
        j	= j - 1;

        if(j < 1 || j > ne){
            photar[i]	= 0.;
        }
        else{
            // ph/cm^2/s/keV
            phflux = phot[j] + (phot[jplus1]-phot[j])*(emid-energ[j])/(energ[jplus1]-energ[j]);
            photar[i] = phflux*(ear[iplus1] - ear[i]);
        }
    }
}
