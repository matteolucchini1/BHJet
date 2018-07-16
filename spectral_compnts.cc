#include "agnjet.hh"

using namespace std;

/**
 * spline_dbl arrays for adding to big total array for z<zmax, and add BB flux, nutot in Hz
 */

void spcomponents(int infosw,int plotsw,int ne,int njet,int nz,int nsyn,int ncom,double zed[],double zcut,
double zsh,double ephxr[],double nusyn[],double synabs[],double nucom[],double comspc[],double nutot[],
double complot[],double presyn[],double postsyn[],double fplot[]){

    int k,m,i;
    double sflx,cflx,z;
    
    double *snu	= new double[nsyn]();
    double *sdump	= new double[nsyn]();
    double *cnu	= new double[ncom]();
    double *cdump	= new double[ncom]();
    
    ofstream zonesFile;    
    
    for(m=0; m<njet; m++){
        for(k=0; k<nz; k++){
            z	= log10(zed[k]);
            for(i=0; i<nsyn; i++){
                snu[i]	= nusyn[(m*nz+k)*nsyn+i];
                sdump[i]= synabs[(m*nz+k)*nsyn+i];
            }
            
            gsl_interp_accel *acc_snu       = gsl_interp_accel_alloc();
            gsl_spline *spline_snu          = gsl_spline_alloc(gsl_interp_cspline, nsyn);
            gsl_spline_init(spline_snu, snu, sdump, nsyn);
            
            gsl_interp_accel *acc_cnu       = gsl_interp_accel_alloc();
            gsl_spline *spline_cnu          = gsl_spline_alloc(gsl_interp_akima, ncom);
            
            if(z<= log10(zcut)){
                for(i=0; i<ncom; i++){
                    cnu[i]	= nucom[(m*nz+k)*ncom+i];
                    cdump[i]= comspc[(m*nz+k)*ncom+i];
                }
                gsl_spline_init(spline_cnu, cnu, cdump, ncom);
            }
            for(i=0; i<(ne-1); i++){
                if(nutot[i] >= snu[0] && nutot[i] <= snu[nsyn-1]){
                    sflx    = gsl_spline_eval(spline_snu, nutot[i], acc_snu);
                }
                else{
                    sflx	= -200.;
                }
                
                if(z> log10(zcut)){
                    cflx	= -200.;
                }
                else{
                    if(nutot[i] >= cnu[0] && nutot[i] <= cnu[ncom-1]){
                        cflx    = gsl_spline_eval(spline_cnu, nutot[i], acc_cnu);
                    }
                    else{
                        cflx	= -200.;
                    }
                }                
                if(plotsw == 1){
                    complot[i]= complot[i] + pow(10,cflx);                    
                    if(z < log10(zsh)){
                        presyn[i]= presyn[i] + pow(10,sflx);
                    }
                    else{
                        postsyn[i]= postsyn[i] + pow(10,sflx);
                    }
                }
                fplot[i]= fplot[i] + pow(10,sflx) + pow(10,cflx);
                
				if (infosw == 1) {
					zonesFile.open("outputs/zones.dat", ios::app);
						if(!zonesFile.is_open()){
							cerr << "*** Error: can't open zones.dat" << endl;
							exit(1);
						} else{
							zonesFile << left << setw(20) << nutot[i] << setw(20) << sflx << setw(20) 
							<< cflx << endl;								
						}
					zonesFile.close();
				}
            }//end ne for loop            
            gsl_spline_free(spline_snu), gsl_interp_accel_free(acc_snu);
            gsl_spline_free(spline_cnu), gsl_interp_accel_free(acc_cnu);
        }//end z for loop
    }//end m[njet] for loop
    
    delete[] snu, delete[] sdump, delete[] cnu, delete[] cdump;    
}
