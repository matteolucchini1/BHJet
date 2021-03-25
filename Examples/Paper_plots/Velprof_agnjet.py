import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import math
import matplotlib.ticker as plticker
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib import rc, rcParams
from scipy import interpolate

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['DejaVu Serif Display']})
plt.rcParams.update({'font.size': 20})

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

velsw1_hratio01 = np.genfromtxt("Dynamics_data/agnjet_velsw1_r012_hratio01.dat")
velsw1_hratio1 = np.genfromtxt("Dynamics_data/agnjet_velsw1_r012_hratio1.dat")
velsw1_hratio10 = np.genfromtxt("Dynamics_data/agnjet_velsw1_r012_hratio10.dat")
velsw0_hratio01 = np.genfromtxt("Dynamics_data/agnjet_velsw0_r012_hratio01.dat")
velsw0_hratio1 = np.genfromtxt("Dynamics_data/agnjet_velsw0_r012_hratio1.dat")
velsw0_hratio10 = np.genfromtxt("Dynamics_data/agnjet_velsw0_r012_hratio10.dat")

collimation_velsw1_hratio01 = np.genfromtxt("Dynamics_data/collimation_velsw1_r012_hratio01.dat")
collimation_velsw1_hratio1 = np.genfromtxt("Dynamics_data/collimation_velsw1_r012_hratio1.dat")
collimation_velsw1_hratio10 = np.genfromtxt("Dynamics_data/collimation_velsw1_r012_hratio10.dat")
collimation_velsw0_hratio01 = np.genfromtxt("Dynamics_data/collimation_velsw0_r012_hratio01.dat")
collimation_velsw0_hratio1 = np.genfromtxt("Dynamics_data/collimation_velsw0_r012_hratio1.dat")
collimation_velsw0_hratio10 = np.genfromtxt("Dynamics_data/collimation_velsw0_r012_hratio10.dat")

R0 = 12.

nozzle_height = [5.2,16.,124.]

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(15,12))

ax1.plot(velsw1_hratio01.T[0],velsw1_hratio01.T[1],linewidth=2.5,color=colors[0],label='h=0.1,iso')
ax1.plot(velsw1_hratio1.T[0],velsw1_hratio1.T[1],linewidth=2.5,color=colors[1],label='h=1,iso')
ax1.plot(velsw1_hratio10.T[0],velsw1_hratio10.T[1],linewidth=2.5,color=colors[2],label='h=10,iso')

ax1.axvline(nozzle_height[0],0,4,color=colors[0],linewidth=3.0,linestyle='dashed')
ax1.axvline(nozzle_height[1],0,4,color=colors[1],linewidth=3.0,linestyle='dashed')
ax1.axvline(nozzle_height[2],0,4,color=colors[2],linewidth=3.0,linestyle='dashed')

ax1.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax1.set_ylabel("Lorenz factor $\\gamma$",fontsize=24)
ax1.legend(loc='lower right',fontsize=24)
ax1.set_xscale('log',base=10)
ax1.set_xlim([2.,1.e8])
ax1.set_ylim([1.08,3.3])

ax2.plot(velsw0_hratio01.T[0],velsw0_hratio01.T[1],linewidth=2.5,color=colors[3],label='h=0.1,ad')
ax2.plot(velsw0_hratio1.T[0],velsw0_hratio1.T[1],linewidth=2.5,color=colors[4],label='h=1,ad')
ax2.plot(velsw0_hratio10.T[0],velsw0_hratio10.T[1],linewidth=2.5,color=colors[5],label='h=10,ad')
ax2.axvline(nozzle_height[0],0,4,color=colors[3],linewidth=3.0,linestyle='dashed')
ax2.axvline(nozzle_height[1],0,4,color=colors[4],linewidth=3.0,linestyle='dashed')
ax2.axvline(nozzle_height[2],0,4,color=colors[5],linewidth=3.0,linestyle='dashed')

ax2.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax2.set_yticklabels([])
ax2.legend(loc='lower right',fontsize=24)
ax2.set_xscale('log',base=10)
ax2.set_xlim([2.,1.e8])
ax2.set_ylim([1.08,3.3])


ax3.plot(collimation_velsw1_hratio01.T[0],collimation_velsw1_hratio01.T[1],linewidth=2.5,color=colors[0],label='h=0.1,iso')
ax3.plot(collimation_velsw1_hratio1.T[0],collimation_velsw1_hratio1.T[1],linewidth=2.5,color=colors[1],label='h=1,iso')
ax3.plot(collimation_velsw1_hratio10.T[0],collimation_velsw1_hratio10.T[1],linewidth=2.5,color=colors[2],label='h=10,iso')
ax3.axvline(nozzle_height[0],0,1e3,color=colors[0],linewidth=3.0,linestyle='dashed')
ax3.axvline(nozzle_height[1],0,1e3,color=colors[1],linewidth=3.0,linestyle='dashed')
ax3.axvline(nozzle_height[2],0,1e3,color=colors[2],linewidth=3.0,linestyle='dashed')

ax3.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax3.set_ylabel("Jet Radius ($R_{\\rm g}$)",fontsize=24)
ax3.set_xlim([2.,5.e3])
ax3.set_ylim([10.,1.e3])
ax3.set_xscale('log',base=10)
ax3.set_yscale('log',base=10)


ax4.plot(collimation_velsw0_hratio01.T[0],collimation_velsw0_hratio01.T[1],linewidth=2.5,color=colors[3],label='h=0.1,ad')
ax4.plot(collimation_velsw0_hratio1.T[0],collimation_velsw0_hratio1.T[1],linewidth=2.5,color=colors[4],label='h=1,ad')
ax4.plot(collimation_velsw0_hratio10.T[0],collimation_velsw0_hratio10.T[1],linewidth=2.5,color=colors[5],label='h=10,ad')
ax4.axvline(nozzle_height[0],0,1e3,color=colors[3],linewidth=3.0,linestyle='dashed')
ax4.axvline(nozzle_height[1],0,1e3,color=colors[4],linewidth=3.0,linestyle='dashed')
ax4.axvline(nozzle_height[2],0,1e3,color=colors[5],linewidth=3.0,linestyle='dashed')

ax4.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax4.set_yticklabels([])
ax4.axes.get_yaxis().set_visible(False)
ax4.set_xlim([2.,5.e3])
ax4.set_ylim([10.,1.e3])
ax4.set_xscale('log',base=10)
ax4.set_yscale('log',base=10)

fig.tight_layout()
fig.subplots_adjust(wspace=0)

fig.savefig('agnjet_dynamics.pdf')

plt.show()


