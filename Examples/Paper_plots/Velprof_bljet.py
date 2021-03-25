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

r0 = 12.
z0 = 1.*r0+2.
rho = 0.15
gamma_0 = 1.09
zmax = 1.e8

gamma_acc = np.array([4,8,18])
z_acc = np.array([1e3,1e4,1e5])
sigma_f = np.array([1.,0.5,0.1])

sigma_0_1 = (1.+sigma_f[0])*gamma_acc[0]/gamma_0-1.
sigma_0_2 = (1.+sigma_f[1])*gamma_acc[1]/gamma_0-1.
sigma_0_3 = (1.+sigma_f[2])*gamma_acc[2]/gamma_0-1.

z_1 = np.logspace(0,np.log10(z_acc[0]),30)
z_2 = np.logspace(0,np.log10(z_acc[1]),40)
z_3 = np.logspace(0,np.log10(z_acc[2]),50)
z_1_post = np.logspace(np.log10(z_acc[0]),np.log10(zmax),70)
z_2_post = np.logspace(np.log10(z_acc[1]),np.log10(zmax),60)
z_3_post = np.logspace(np.log10(z_acc[2]),np.log10(zmax),50)

gamma_1 = gamma_0 + (gamma_acc[0]-gamma_0)*(z_1**(1./2.)-z0**(1./2.))/(z_acc[0]**(1./2.)-z0**(1./2.))
gamma_2 = gamma_0 + (gamma_acc[1]-gamma_0)*(z_2**(1./2.)-z0**(1./2.))/(z_acc[1]**(1./2.)-z0**(1./2.))
gamma_3 = gamma_0 + (gamma_acc[2]-gamma_0)*(z_3**(1./2.)-z0**(1./2.))/(z_acc[2]**(1./2.)-z0**(1./2.))
gamma_1_post = np.zeros(70)+gamma_acc[0]
gamma_2_post = np.zeros(60)+gamma_acc[1]
gamma_3_post = np.zeros(50)+gamma_acc[2]

sigma_1 = (gamma_0/gamma_1)*(1.+sigma_0_1)-1.
sigma_2 = (gamma_0/gamma_2)*(1.+sigma_0_2)-1.
sigma_3 = (gamma_0/gamma_3)*(1.+sigma_0_3)-1.
sigma_1_post = np.zeros(70)+sigma_f[0]
sigma_2_post = np.zeros(60)+sigma_f[1]
sigma_3_post = np.zeros(50)+sigma_f[2]

theta_1 = rho/gamma_1
theta_2 = rho/gamma_2
theta_3 = rho/gamma_3
r_1 = r0 + (z_1-z0)*np.tan(theta_1)
r_2 = r0 + (z_2-z0)*np.tan(theta_2)
r_3 = r0 + (z_3-z0)*np.tan(theta_3)
r_1_post = r_1[29] + (z_1_post-z_acc[0])*np.tan(rho/gamma_acc[0])
r_2_post = r_2[39] + (z_2_post-z_acc[1])*np.tan(rho/gamma_acc[1])
r_3_post = r_3[49] + (z_3_post-z_acc[2])*np.tan(rho/gamma_acc[2])

fig, ((ax1,ax2)) = plt.subplots(1,2,figsize=(15,6))

ax1.plot(z_1,gamma_1,linewidth=2.5,color=colors[0])
ax1.plot(z_2,gamma_2,linewidth=2.5,color=colors[1])
ax1.plot(z_3,gamma_3,linewidth=2.5,color=colors[2])
ax1.plot(z_1_post,gamma_1_post,linewidth=2.5,color=colors[0])
ax1.plot(z_2_post,gamma_2_post,linewidth=2.5,color=colors[1])
ax1.plot(z_3_post,gamma_3_post,linewidth=2.5,color=colors[2])

ax1.plot(z_1,sigma_1,linewidth=2.5,color=colors[0],linestyle='dashed')
ax1.plot(z_2,sigma_2,linewidth=2.5,color=colors[1],linestyle='dashed')
ax1.plot(z_3,sigma_3,linewidth=2.5,color=colors[2],linestyle='dashed')
ax1.plot(z_1_post,sigma_1_post,linewidth=2.5,color=colors[0],linestyle='dashed')
ax1.plot(z_2_post,sigma_2_post,linewidth=2.5,color=colors[1],linestyle='dashed')
ax1.plot(z_3_post,sigma_3_post,linewidth=2.5,color=colors[2],linestyle='dashed')

ax1.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax1.set_ylabel("Lorenz factor $\\gamma$/Magnetisation $\\sigma$",fontsize=24)
ax1.set_xscale("log",base=10)
ax1.set_xlim([z0,zmax])
ax1.set_ylim([0.01,gamma_acc[2]+1.])

ax2.plot(z_1,r_1,linewidth=2.5,color=colors[0],label='$z_{\\rm acc}=10^{3},\\gamma_{\\rm acc}=4$')
ax2.plot(z_2,r_2,linewidth=2.5,color=colors[1],label='$z_{\\rm acc}=10^{4},\\gamma_{\\rm acc}=8$')
ax2.plot(z_3,r_3,linewidth=2.5,color=colors[2],label='$z_{\\rm acc}=10^{5},\\gamma_{\\rm acc}=18$')
ax2.plot(z_1_post,r_1_post,linewidth=2.5,color=colors[0])
ax2.plot(z_2_post,r_2_post,linewidth=2.5,color=colors[1])
ax2.plot(z_3_post,r_3_post,linewidth=2.5,color=colors[2])

ax2.set_xlabel("Distance from BH ($R_{\\rm g}$)",fontsize=24)
ax2.set_ylabel("Jet radius ($R_{\\rm g}$)",fontsize=24)
ax2.legend(loc="upper left",fontsize=24)
ax2.set_xscale("log",base=10)
ax2.set_yscale("log",base=10)
ax2.set_xlim([z0,1e5])
ax2.set_ylim([0.8*r0,3e3])

plt.tight_layout()
fig.savefig('bljet_dynamics.pdf')
plt.show()
