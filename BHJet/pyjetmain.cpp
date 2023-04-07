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
#include <fenv.h>
#include <omp.h>


using namespace std;

extern void jetmain(double *ear,int ne,double *param,double *photeng,double *photspec);

extern "C" void pyjetmain(double *ear,int ne,double *param,double *photeng,double *photspec){
    jetmain(ear,ne,param,photeng,photspec);
}				// ----------  end of function main  ----------

