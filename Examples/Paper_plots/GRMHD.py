import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import math
import matplotlib.ticker as plticker
from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.patches import Polygon
import matplotlib.ticker as plticker
from matplotlib import rc, rcParams
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['DejaVu Serif Display']})
plt.rcParams.update({'font.size': 20})

#in order: blue, orange, green, red, purple, brown, pink, grey
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

agnjet_sheath = np.genfromtxt("Profiles_agnjet_sheath.dat")
bljet_sheath = np.genfromtxt("Profiles_bljet_sheath.dat")
bljet_spine = np.genfromtxt("Profiles_bljet_spine.dat")

MF13_B = np.genfromtxt("MF13/Sheath_B.dat")
MF13_n = np.genfromtxt("MF13/Sheath_rho.dat")
MF13_T = np.genfromtxt("MF13/Sheath_theta.dat")

MFSG14_B = np.genfromtxt("MFSG14/Sheath_B.dat")
MFSG14_n = np.genfromtxt("MFSG14/Sheath_rho.dat")
MFSG14_T = np.genfromtxt("MFSG14/Sheath_theta.dat")

MFN16_sheath_B_ModelA = np.genfromtxt("MFN16/Sheath_B_ModelA.dat")
MFN16_sheath_n_ModelA = np.genfromtxt("MFN16/Sheath_rho_ModelA.dat")
MFN16_sheath_T_ModelA = np.genfromtxt("MFN16/Sheath_theta_ModelA.dat")

MFN16_spine_B_ModelA = np.genfromtxt("MFN16/Spine_B_ModelA.dat")
MFN16_spine_n_ModelA = np.genfromtxt("MFN16/Spine_rho_ModelA.dat")
MFN16_spine_T_ModelA = np.genfromtxt("MFN16/Spine_theta_ModelA.dat")

MFN16_sheath_B_ModelG = np.genfromtxt("MFN16/Sheath_B_ModelG.dat")
MFN16_sheath_n_ModelG = np.genfromtxt("MFN16/Sheath_rho_ModelG.dat")
MFN16_sheath_T_ModelG = np.genfromtxt("MFN16/Sheath_theta_ModelG.dat")

MFN16_spine_B_ModelG = np.genfromtxt("MFN16/Spine_B_ModelG.dat")
MFN16_spine_n_ModelG = np.genfromtxt("MFN16/Spine_rho_ModelG.dat")
MFN16_spine_T_ModelG = np.genfromtxt("MFN16/Spine_theta_ModelG.dat")

#TBD: two sets of plots. One comparing sheaths, and one comparing spines

line_x1 = np.zeros(100)+16.
line_x2 = np.zeros(100)+14.
line_y = np.logspace(-6,0,100)

#First compare sheaths:
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(15.,12.))
ax1.loglog(10.**MF13_B.T[0],10.**MF13_B.T[1]*1e1,linewidth=2.5,label="MF 13")
ax1.loglog(10.**MFSG14_B.T[0],10.**MFSG14_B.T[1]*1e-3,linewidth=2.5,label="MFSG 14")
ax1.loglog(10.**MFN16_sheath_B_ModelA.T[0],10.**MFN16_sheath_B_ModelA.T[1],linewidth=2.5,label="MFN 16, sheath Mod. A")
ax1.loglog(10.**MFN16_sheath_B_ModelG.T[0],10.**MFN16_sheath_B_ModelG.T[1],linewidth=2.5,label="MFN 16, sheath Mod. G")
ax1.loglog(agnjet_sheath.T[0],agnjet_sheath.T[2]*1e-9,linewidth=2.5,label='agnjet')
ax1.loglog(bljet_sheath.T[0],bljet_sheath.T[2]*1e-9,linewidth=2.5,label='bljet')
ax1.loglog(line_x1,line_y,linestyle='dashed',color='black',linewidth=2.5)
ax1.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax1.set_ylabel("B-field (arb. units)",fontsize=24)
ax1.set_ylim([5e-4,3e-1])
ax1.set_xlim([5.,0.8e3])
ax1.legend(loc="best",fontsize=18)

ax2.loglog(10.**MF13_n.T[0],10.**MF13_n.T[1],linewidth=2.5)
ax2.loglog(10.**MFSG14_n.T[0],10.**MFSG14_n.T[1]*1e-8,linewidth=2.5)
ax2.loglog(10.**MFN16_sheath_n_ModelA.T[0],10.**MFN16_sheath_n_ModelA.T[1],linewidth=2.5)
ax2.loglog(10.**MFN16_sheath_n_ModelG.T[0],10.**MFN16_sheath_n_ModelG.T[1],linewidth=2.5)
ax2.loglog(agnjet_sheath.T[0],agnjet_sheath.T[3]*1e-20,linewidth=2.5,)
ax2.loglog(bljet_sheath.T[0],bljet_sheath.T[3]*1e-20,linewidth=2.5)
ax2.loglog(line_x1,line_y,linestyle='dashed',color='black',linewidth=2.5)
ax2.set_ylim([1e-6,3e-1])
ax2.set_xlim([5.,0.8e3])
ax2.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax2.set_ylabel("Number density (arb. units)",fontsize=24)

#then compare spines:
ax3.loglog(10.**MFN16_spine_B_ModelA.T[0],10.**MFN16_spine_B_ModelA.T[1],linewidth=2.5,label="MFN 16, spine Mod. A")
ax3.loglog(10.**MFN16_spine_B_ModelG.T[0],10.**MFN16_spine_B_ModelG.T[1],linewidth=2.5,label="MFN 16, spine Mod. G")
ax3.loglog(bljet_spine.T[0],bljet_spine.T[2]*1e-5,linewidth=2.5,label='bljet')
ax3.loglog(line_x2,line_y,linestyle='dashed',color='black',linewidth=2.5)
ax3.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax3.set_ylabel("B-field (arb. units)",fontsize=24)
ax3.set_ylim([1e-4,3e-1])
ax3.set_xlim([5.,0.8e3])
ax3.legend(loc="best",fontsize=18)

ax4.loglog(10.**MFN16_spine_n_ModelA.T[0],10.**MFN16_spine_n_ModelA.T[1],linewidth=2.5)
ax4.loglog(10.**MFN16_spine_n_ModelG.T[0],10.**MFN16_spine_n_ModelG.T[1],linewidth=2.5)
ax4.loglog(bljet_spine.T[0],bljet_spine.T[3]*0.3e-11,linewidth=2.5)
ax4.loglog(line_x2,line_y,linestyle='dashed',color='black',linewidth=2.5)
ax4.set_ylim([1e-5,6e-4])
ax4.set_xlim([5.,2e2])
ax4.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax4.set_ylabel("Number density (arb. units)",fontsize=24)

plt.tight_layout()

plt.savefig("MHD_comparison1.pdf")

CH19=np.load('Chatterjee19/Chatterjee19.npz')
#jet speeds now
fig2, ax5 = plt.subplots(1,1,figsize=(7.5,6.))
rad=CH19['radius']
LF_gam=CH19['gamma']
ax5.loglog(rad,LF_gam,linewidth=2.5,label="CLTM 19")
ax5.loglog(bljet_spine.T[0],bljet_spine.T[4],linewidth=2.5,label='Bljet spine')
ax5.loglog(bljet_sheath.T[0],bljet_sheath.T[4],linewidth=2.5,label='Bljet sheath')
ax5.loglog(agnjet_sheath.T[0],agnjet_sheath.T[4],linewidth=2.5,label='Agnjet sheath')
ax5.set_xlim(1e1,1e5)
ax5.set_ylim(1.,27.5)
ax5.legend(loc='upper right',fontsize=18)
ax5.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax5.set_ylabel("Lorenz factor $\\gamma(z)$",fontsize=24)

plt.tight_layout()

plt.savefig("MHD_comparison2.pdf")

plt.show()
