#include "agnjet.hh"

using namespace std;


void create_plotFiles(int plotsw,int infosw){

    ofstream fluxplotFile, complotFile, presynFile, postsynFile, bbplotFile, extphFile, zonesFile;
    
    if(plotsw == 1){
        fluxplotFile.open("outputs/total.dat",ios::trunc);
        if(!fluxplotFile.is_open()){
            cerr << "*** Error: can't open outputs/total.dat" << endl;
            exit(1);
		} else{
                fluxplotFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Total Emission" << endl;
        }
        fluxplotFile.close();
            
        complotFile.open("outputs/com.dat",ios::trunc);
        if(!complotFile.is_open()){
            cerr << "*** Error: can't open outputs/com.dat" << endl;
            exit(1);
        } else{
            complotFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Compton Emission" << endl;
        }
        complotFile.close();
            
        presynFile.open("outputs/presyn.dat",ios::trunc);
        if(!presynFile.is_open()){
            cerr << "*** Error: can't open outputs/presyn.dat" << endl;
            exit(1);
        } else{
            presynFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Presyn Emission" << endl;
        }
        presynFile.close();
            
        postsynFile.open("outputs/postsyn.dat",ios::trunc);
        if(!postsynFile.is_open()){
            cerr << "*** Error: can't open outputs/postsyn.dat" << endl;
            exit(1);
        } else{
            postsynFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Postsyn Emission" << endl;
        }
        postsynFile.close();
            
        bbplotFile.open("outputs/bb.dat",ios::trunc);
        if(!bbplotFile.is_open()){
            cerr << "*** Error: can't open outputs/bb.dat" << endl;
            exit(1);
        } else{
            bbplotFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Black-body Emission" << endl;
        }
        bbplotFile.close();
        
        extphFile.open("outputs/extph.dat",ios::trunc);
        if(!extphFile.is_open()){
            cerr << "*** Error: can't open outputs/extph.dat" << endl;
            exit(1);
        } else{
            extphFile << left << setw(20) << "#nu [Hz]" << setw(20) << "External Black-body Emission" << endl;
        }
        extphFile.close();
            
        if(infosw == 1){
			zonesFile.open("outputs/zones.dat",ios::trunc);
			if(!zonesFile.is_open()){
				cerr << "*** Error: can't open zones.dat" << endl;
				exit(1);
			}
			else{
				zonesFile << left << setw(20) << "#nu [Hz]" << setw(20) << "Total Emission" << endl;
			}
			zonesFile.close();
		}
	}
}

void write_plotFiles(int plotsw,int disksw,int compsw,double tin,double rin,double rout,double dist,
double inclin,double normbb,double tbb,double nutot,double &fplot,double complot,double presyn,double postsyn, 
double bbplot,double extphplot){

    ofstream fluxplotFile, complotFile, presynFile, postsynFile, bbplotFile, extphFile;
    
    double blim,ulim,bbtrm,bbflx,frq;        					
    
    blim	= log(rin);
    ulim	= log(rout);
    frq     = pow(10,nutot);

    //Calculating disk contribution only
   	if(disksw == 1){        
		bbearthint(blim, ulim, frq, tin, rin, dist, inclin, bbflx);
        bbflx = bbflx/mjy;        
        fplot= fplot + bbflx;
        if(plotsw == 1){
            bbplot = bbflx;
        }
    }
    
    if(compsw == 1){
    	bbtrm = normbb*2.*herg*frq*frq*frq/(cee*cee*mjy*4.*pi*pow(dist,2)*(exp(herg*frq/(kboltz*tbb))-1.));    	
    	fplot = fplot + bbtrm;
		extphplot = bbtrm;         
    }
    
    if(plotsw == 1){
		fluxplotFile.open("outputs/total.dat", ios::app);
		if(!fluxplotFile.is_open()){
        	cerr << "*** Error: can't open outputs/total.dat" << endl;
            exit(1);
		}
        else{
        	fluxplotFile << left << setw(20) << nutot << setw(20) << log10(fplot) << endl;
        }
        fluxplotFile.close();   

		complotFile.open("outputs/com.dat", ios::app);
        if(!complotFile.is_open()){
        	cerr << "*** Error: can't open outputs/com.dat" << endl;
        	exit(1);
		}
        else{
        	complotFile << left << setw(20) << nutot << setw(20) << log10(complot) << endl;
        }
        complotFile.close();        
        
		presynFile.open("outputs/presyn.dat", ios::app);
        if(!presynFile.is_open()){
        	cerr << "*** Error: can't open outputs/presyn.dat" << endl;
        	exit(1);
		}
        else{
        	presynFile << left << setw(20) << nutot << setw(20) << log10(presyn) << endl;
		}
        presynFile.close();
        
        postsynFile.open("outputs/postsyn.dat", ios::app);
        if(!postsynFile.is_open()){
        	cerr << "*** Error: can't open outputs/postsyn.dat" << endl;
            exit(1);
		}
        else{
        	postsynFile << left << setw(20) << nutot << setw(20) << log10(postsyn) << endl;
        }
        postsynFile.close();        

		bbplotFile.open("outputs/bb.dat", ios::app);
        if(!bbplotFile.is_open()){
        	cerr << "*** Error: can't open outputs/bb.dat" << endl;
            exit(1);
		}
        else{
        	bbplotFile << left << setw(20) << nutot << setw(20) << log10(bbplot) << endl;
		}
        bbplotFile.close();
        
        extphFile.open("outputs/extph.dat", ios::app);
        if(!extphFile.is_open()){
        	cerr << "*** Error: can't open outputs/extph.dat" << endl;
            exit(1);
		}
        else{
        	extphFile << left << setw(20) << nutot << setw(20) << log10(extphplot) << endl;
		}
        extphFile.close();         
    }//end if plotsw    
}
