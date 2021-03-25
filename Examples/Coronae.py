import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import math
import matplotlib.ticker as plticker
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib import rc, rcParams
import matplotlib.font_manager

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['DejaVu Serif Display']})
plt.rcParams.update({'font.size': 20})

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

mjy = 1.e-26
kev=2.41e17
kpc = 3.e21
conv = 200.*(kpc)**2
#kev = 1.
Mbh = 10.

Ledd = Mbh*1.25e38
Ldisk = 1e-2

Disk = np.genfromtxt("Output/Disk.dat")
IC_Tau260Te90 = np.genfromtxt("Output/IC_Tau260Te90.dat")
IC_Tau076Te90 = np.genfromtxt("Output/IC_Tau076Te90.dat")
IC_Tau019Te90 = np.genfromtxt("Output/IC_Tau019Te90.dat")
IC_Tau260Te900 = np.genfromtxt("Output/IC_Tau260Te900.dat")
IC_Tau076Te900 = np.genfromtxt("Output/IC_Tau076Te900.dat")
IC_Tau019Te900 = np.genfromtxt("Output/IC_Tau019Te900.dat")

CompPS_frq = np.genfromtxt("CompPS_table/CompPS_freq.dat");
CompPS_flx_1 = np.genfromtxt("CompPS_table/Sphere/Te_90/tau260.dat")
CompPS_flx_2 = np.genfromtxt("CompPS_table/Sphere/Te_90/tau076.dat")
CompPS_flx_3 = np.genfromtxt("CompPS_table/Sphere/Te_90/tau019.dat")
CompPS_flx_4 = np.genfromtxt("CompPS_table/Sphere/Te_900/tau260.dat")
CompPS_flx_5 = np.genfromtxt("CompPS_table/Sphere/Te_900/tau076.dat")
CompPS_flx_6 = np.genfromtxt("CompPS_table/Sphere/Te_900/tau019.dat")


fig2, (ax5,ax6) = plt.subplots(1,2,figsize=(15.,6))

#ax5.loglog(Disk.T[0]/kev,Disk.T[0]*Disk.T[1],linewidth=1.5,color="black")
ax5.loglog(IC_Tau260Te90.T[0]/kev,IC_Tau260Te90.T[0]*IC_Tau260Te90.T[1],linewidth=2.5,color=colors[0],label='$\\tau = 2.6$')
ax5.loglog(CompPS_frq/kev,conv*CompPS_flx_1,linewidth=3.5,color=colors[0],linestyle='dashed')
ax5.loglog(IC_Tau076Te90.T[0]/kev,IC_Tau076Te90.T[0]*IC_Tau076Te90.T[1],linewidth=2.5,color=colors[1],label='$\\tau = 0.76$')
ax5.loglog(CompPS_frq/kev,0.2*conv*CompPS_flx_2,linewidth=3.5,color=colors[1],linestyle='dashed')
ax5.loglog(IC_Tau019Te90.T[0]/kev,IC_Tau019Te90.T[0]*IC_Tau019Te90.T[1],linewidth=2.5,color=colors[2],label='$\\tau = 0.19$')
ax5.loglog(CompPS_frq/kev,0.05*conv*CompPS_flx_3,linewidth=3.5,color=colors[2],linestyle='dashed')
ax5.set_ylim([1e-5*Ledd*Ldisk,2.*Ledd*Ldisk])
ax5.set_xlim([2e17/kev,5.e20/kev])
ax5.set_xlabel("Energy (keV)",fontsize=24)
ax5.set_ylabel("Luminosity (Arb. units)",fontsize=24)
ax5.legend(loc='lower left',fontsize=24)

#ax6.loglog(Disk.T[0]/kev,Disk.T[0]*Disk.T[1],linewidth=1.5,color="black")
ax6.loglog(IC_Tau260Te900.T[0]/kev,IC_Tau260Te900.T[0]*IC_Tau260Te900.T[1],linewidth=2.5,color=colors[3],label='$\\tau = 2.6$')
ax6.loglog(CompPS_frq/kev,0.026*conv*CompPS_flx_4,linewidth=3.5,color=colors[3],linestyle='dashed')
ax6.loglog(IC_Tau076Te900.T[0]/kev,IC_Tau076Te900.T[0]*IC_Tau076Te900.T[1],linewidth=2.5,color=colors[4],label='$\\tau = 0.76$')
ax6.loglog(CompPS_frq/kev,0.17*conv*CompPS_flx_5,linewidth=3.5,color=colors[4],linestyle='dashed')
ax6.loglog(IC_Tau019Te900.T[0]/kev,IC_Tau019Te900.T[0]*IC_Tau019Te900.T[1],linewidth=2.5,color=colors[5],label='$\\tau = 0.19$')
ax6.loglog(CompPS_frq/kev,0.8*conv*CompPS_flx_6,linewidth=3.5,color=colors[5],linestyle='dashed')
ax6.set_ylim([1e-5*Ledd*Ldisk,2.*Ledd*Ldisk])
ax6.set_xlim([2e17/kev,1.e21/kev])
ax6.set_xlabel("Energy (keV)",fontsize=24)
ax6.yaxis.tick_right()
ax6.yaxis.set_ticklabels('')
ax6.yaxis.set_label_position("right")
ax6.legend(loc='lower right',fontsize=24)

fig2.tight_layout()
fig2.subplots_adjust(wspace=0)	

plt.savefig("corona_examples.pdf")
plt.show()
