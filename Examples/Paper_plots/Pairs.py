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

mp = 1.67e-24
me = 9.10e-28

sigma_bljet = np.array([2.,7.,20.])
av_gamma = np.array([1.5,3.,10.])

beta = np.logspace(-4,1,200)

eta_agnjet_1 = (mp/me)*1./(av_gamma[0]*(1.+1./beta))
eta_agnjet_2 = (mp/me)*1./(av_gamma[2]*(1.+1./beta))

eta_bljet_1 = (mp/me)*(sigma_bljet[0]*beta)/(av_gamma[0]*(2.-sigma_bljet[0]*beta*4./3.))
eta_bljet_2 = (mp/me)*(sigma_bljet[1]*beta)/(av_gamma[0]*(2.-sigma_bljet[1]*beta*4./3.))
eta_bljet_3 = (mp/me)*(sigma_bljet[2]*beta)/(av_gamma[0]*(2.-sigma_bljet[2]*beta*4./3.))

eta_bljet_4 = (mp/me)*(sigma_bljet[0]*beta)/(av_gamma[2]*(2.-sigma_bljet[0]*beta*4./3.))

shade = np.logspace(-5,-0.15,200)
shade2 = np.logspace(-0.15,1.1,100)

shade1 = np.logspace(-5,-1.15,200)
shade3 = np.logspace(-0.15,1,100)
shade4 = np.logspace(-0.70,-0.15,100)
shade5 = np.logspace(-1.15,-0.7,100)

fig, ((ax1,ax2)) = plt.subplots(1,2,figsize=(15,6))

ax1.plot(beta,eta_agnjet_1,linewidth=2.5,color=colors[0],linestyle='dashed')
ax1.plot(beta,eta_bljet_1,linewidth=2.5,color=colors[0],label='$\\langle\\gamma\\rangle=1.5$')

ax1.plot(beta,eta_agnjet_2,linewidth=2.5,color=colors[1],linestyle='dashed')
ax1.plot(beta,eta_bljet_4,linewidth=2.5,color=colors[1],label='$\\langle\\gamma\\rangle=10$')

ax1.fill_between(shade,1e2,2e3,color='#d8dcd6')
ax1.fill_between(shade,0,1.,color='#c7fdb5')
ax1.fill_between(shade2,0,2e3,color=colors[1],alpha=0.3)

ax1.set_xlabel("$\\beta_{\\rm p}$",fontsize=24)
ax1.set_ylabel("Pair content $\\eta$",fontsize=24)
ax1.legend(loc='upper left',fontsize=24)
ax1.set_xlim([1e-4,10.])
ax1.set_ylim([0.3,1.5e3])
ax1.set_xscale("log",base=10)
ax1.set_yscale("log",base=10)

ax2.plot(beta,eta_bljet_1,linewidth=2.5,color=colors[3],label='$\\sigma_0=2$')
ax2.plot(beta,eta_bljet_2,linewidth=2.5,color=colors[4],label='$\\sigma_0=7$')
ax2.plot(beta,eta_bljet_3,linewidth=2.5,color=colors[5],label='$\\sigma_0=20$')

ax2.fill_between(shade1,1e2,2e3,color='#d8dcd6')
ax2.fill_between(shade1,0,1.,color='#c7fdb5')
ax2.fill_between(shade3,0,2e3,color=colors[3],alpha=0.4)
ax2.fill_between(shade4,0,2e3,color=colors[4],alpha=0.4)
ax2.fill_between(shade5,0,2e3,color=colors[5],alpha=0.4)

ax2.legend(loc='upper left',fontsize=24)
ax2.set_xlabel("$\\beta_{\\rm p}$",fontsize=24)
#ax2.set_yticklabels([])
#ax2.axes.get_yaxis().set_visible(False)
ax2.set_xlim([1.1e-4,1.])
ax2.set_ylim([0.3,1.5e3])
ax2.set_xscale("log",base=10)
ax2.set_yscale("log",base=10)

fig.tight_layout()
#fig.subplots_adjust(wspace=0)
fig.savefig('bhjet_content.pdf')

plt.show()
