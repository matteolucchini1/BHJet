import ctypes
import pathlib

libname = "/home/sheridan/GIT_PROJECTS/BHJet/BHJet/BHJet/pyjetmain"
c_lib = ctypes.CDLL(libname)

 
npar = int(27)
ne = int(201)
emin = float(-10)
emax = float(10)
einc = (emax-emin)/ne

ebins = (ctypes.c_double * ne)()
param = (ctypes.c_double * npar)()
spec = (ctypes.c_double * (ne-1))()
dumarr = (ctypes.c_double * (ne-1))()
for i in range(ne):
    ebins[i]= pow(10,(emin+i*einc))
    

param[0]=1e9		#//black hole mass
param[1]=2.5		#		//viewing angle
param[2]=543e3		#		//distance (kpc)
param[3]=0              #//redshift
param[4]=.9e-2		#//jet power (Eddington units)	
param[5]=18		#	//r0 (rg)
param[6]=600		#	//shock distance (rg)
param[7]=600		#	//magnetic acceleration distance (rg)
param[8]=6.6e5		#	//maximum distance (rg) - z_max
param[9]=300		#	//electron temperature (electron peak gamma) - eltemp
param[10]=.1            # // % non thermal at dissipation region
param[11]=0             # // change plfrac over distance
param[12]=1.95		#		//non-thermal slope - pspec   
param[13]=12		#		//shock heating
param[14]=10.		#		//dynamic time scale parameter
param[15]=3e-6		#	//particle acceleratino time scale - fsc
param[16]=1		#//plasma beta
param[17]=0.025		#		//acceleration sigma
param[18]=0.01		#	//disk temperature/accretion rate
param[19]=50		#		//disk inner radius (rg)
param[20]=1e3		#		//disk outer radius (rg)
param[21]=3300		#		//compar1
param[22]=3e-10		#	//compar2
param[23]=3.e10		#	//compar3
param[24]=0		#		//compsw
param[25]=15		#		//velsw
param[26]=1		#		//infosw  

c_lib.pyjetmain(ebins,ne-1,param,spec,dumarr)

import matplotlib.pyplot as plt
import numpy as np

plt.title('Total emission from BHJet')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Energy (mJy)')

plt.loglog(spec, dumarr, label='BHJet output')

plt.show()
