// =====================================================================================
// 
//       Filename:  jetwrap.cc
// 
//    Description:  c++ wrapper for agnjet routines. Needs xrbjet.{cc/hh} to work.
// 
//        Version:  1.0
//        Created:  05/31/2012 16:30:16
//       Revision:  none
//       Compiler:  g++ -Wall -m32 -g -O0 -o jetwrap mj_aven.cc thermal.cc agnjet.cc jetwrap.cc -lgsl -lgslcblas -lm
// 
//         Author:  Samia Drappeau (sd), drappeau.samia@gmail.com
//        Company:  
// 
// =====================================================================================

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

#include "agnjet.hh"
using namespace std;

void read_params(string file, double *pars);

extern void bhjet(double* ear, int ne, double *params, double *phot_spect, double *photeng);

int main(){
	int time_elaps	= 0;    
	int npar	= 27;
		
	int ne		= 300;
	double emin	= -10.;
	double emax	= 10;
	double einc	= (emax-emin)/ne;

	double *ebins	= new double[ne]();
	double *param	= new double[npar]();
	double *spec	= new double[ne]();
	double *dumarr	= new double[ne]();

	for(int i=0; i<ne; i++){
		ebins[i]= pow(10,(emin+i*einc));
	}

	read_params("Input/ip.dat", param);
    
	bhjet(ebins,ne,param,dumarr,spec);
    
    for (int i=0; i<npar; i++) {
        cout << "param["<< i << "]= " << param[i] << endl;
    }

    delete[] ebins, delete[] param, delete[] spec, delete[] dumarr;

    time_elaps = clock();
    cout << "Total running time: " << (double) time_elaps/CLOCKS_PER_SEC << " seconds" << endl;

    system("python AGNplot2.py");
	return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

/**
 * Read input file
 *
 * @param file		Parameters file
 *
 * @return pars         Parameters
 *
 */
void read_params(string file, double *pars){
	ifstream inFile;
	inFile.open(file.c_str());
	string line;
	int line_nb = 0;
	if(!inFile){
		cerr << "Can't open input file" << endl;
		exit(1);
	}
	while(getline(inFile, line)){
		//the following line trims white space from the beginning of the string
		line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace)))); 

		if(line[0] == '#'){
                        continue;
                }
                else{
                        pars[line_nb] = atof(line.c_str());
                        line_nb++;
                }
	}
        inFile.close();
}
