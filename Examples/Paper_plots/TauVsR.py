import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import math
import matplotlib.ticker as plticker
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib import rc, rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['DejaVu Serif Display']})
plt.rcParams.update({'font.size': 20})

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

Mbh = 10.
Msun = 1.9e33
Eddlum = 1.26e38*Mbh
gconst = 6.67e-8
c = 3e10
Rg = gconst*Mbh*Msun/c**2;
me = 9.1e-28
mp = 1.67e-24
sigt = 6.65e-25
eta = 10.
g0 = 1.09
beta0 = 0.43
avgamma = 3.
pbeta = 0.01
k = avgamma*eta*(me/mp)
print(k)
R = np.logspace(0,3,100)
R = R*Rg

Mdot = np.array([1.,1.e-1,1.e-2,1.e-3,1.e-4])
Mdot = Mdot*Eddlum

const = 2.*math.pi*g0*beta0*avgamma*me*c**3*(1.+1./pbeta+1./k)

ne_1 = Mdot[0]/(const*R**2) 
ne_2 = Mdot[1]/(const*R**2)
ne_3 = Mdot[2]/(const*R**2)
ne_4 = Mdot[3]/(const*R**2)
ne_5 = Mdot[4]/(const*R**2)

tau_1 = ne_1*R*sigt
tau_2 = ne_2*R*sigt
tau_3 = ne_3*R*sigt
tau_4 = ne_4*R*sigt
tau_5 = ne_5*R*sigt

fig, (ax5) = plt.subplots(1,1,figsize=(7.5,6))
ax5.loglog(R/Rg,tau_1,color=colors[0],linewidth=2.5)
ax5.loglog(R/Rg,tau_2,color=colors[1],linewidth=2.5,linestyle=(0,(5,1)))
ax5.loglog(R/Rg,tau_3,color=colors[2],linewidth=2.5,linestyle=(0,(5,3)))
ax5.loglog(R/Rg,tau_4,color=colors[3],linewidth=2.5,linestyle=(0,(5,5)))
ax5.loglog(R/Rg,tau_5,color=colors[4],linewidth=2.5,linestyle=(0,(5,7)))

ax5.fill_between(R/Rg,0.05,3.,color='#d8dcd6')

#ax5.legend(loc='best',fontsize=20)
ax5.set_xlim([1.,1000.])
ax5.set_ylim([1e-4,1e1])
ax5.set_xlabel("Jet radius ($R_{\\rm g}$)",fontsize=24)
ax5.set_ylabel("Thomson optical depth $\\tau$",fontsize=24)
	
fig.tight_layout()





plt.show()
