import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import math
import matplotlib.ticker as plticker
from matplotlib import cm
from matplotlib.colors import LogNorm

total = np.genfromtxt("outputs/total.dat")
bb = np.genfromtxt("outputs/bb.dat")
extph = np.genfromtxt("outputs/extph.dat")
pre = np.genfromtxt("outputs/presyn.dat")
post = np.genfromtxt("outputs/postsyn.dat")
compton = np.genfromtxt("outputs/com.dat")

zones = np.genfromtxt("outputs/zones.dat")

fr = np.zeros(298)
sf = np.zeros(298)
cf = np.zeros(298)

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(23,6))

i = 0
n = 70
colors = pl.cm.jet(np.linspace(0.25,1,n))

zmin = 6.0
zmax = 3.e5
#type in zmin/zmax from agnjet

zheight=np.zeros(5) # the array that we use to print the labels in colorbar
for i in range(0,len(zheight)):
    zheight[i] = pow(10.,math.log10(zmax/zmin)*(i*.25) + math.log10(zmin))
    print("z[",i,"] = %.2e" %zheight[i], "r_g")

ax1.plot(10.**bb.T[0], 10**(bb.T[0]+bb.T[1]-3.-23.), linewidth=1.5, color='#FB02ED', linestyle='dashed')
ax1.plot(10.**extph.T[0], 10**(extph.T[0]+extph.T[1]-3.-23.), linewidth=1.5, color='#4b006e', linestyle='dashed')
ax1.plot(10.**pre.T[0], 10**(pre.T[0]+pre.T[1]-3.-23.), linewidth=2.5, color='#3AEFF4',zorder=151)
ax1.plot(10.**post.T[0], 10**(post.T[0]+post.T[1]-3.-23.), linewidth=2.5, color='#52F856')
ax1.plot(10.**compton.T[0], 10**(compton.T[0]+compton.T[1]-3.-23.), linewidth=2.5, color='#0102EE')
ax1.plot(10.**total.T[0], 10**(total.T[0]+total.T[1]-3.-23.), linewidth=1.5, color='#0B0B0B')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1
	if ((j%10) == 0): #this plots one in ten zones for clarity
		ax1.plot(10.**fr,10.**(fr+sf-3.-23.),linewidth=3.5,color=colors[j],linestyle='dashed',zorder=150-j)
		ax1.plot(10.**fr,10.**(fr+cf-3.-23.),linewidth=3.5,color=colors[j],linestyle='dashed',zorder=150-j)
ax1.set_xlim([1e9,.99e26])
ax1.set_ylim([1.e-14,5.0e-11])
ax1.set_ylabel('$\\nu$F$\\nu$ (erg s$\\rm ^{-1}$ cm$\\rm ^{-2})$', fontsize=18)
ax1.set_xlabel('Frequency (Hz)', fontsize=18)
ax1.set_yscale('log', basey=10)
ax1.set_xscale('log', basex=10)
ax1.tick_params([1.e12,1.e15,1.e18,1.e21,1.e24],fontsize=16)
ax1.tick_params([1.e-14,1.e-13,1.e-12,1.e-11,1.e-10],fontsize=16)

i = 0

ax2.plot(10.**bb.T[0], 10**bb.T[1], linewidth=1.5, color='#FB02ED', linestyle='dashed')
ax2.plot(10.**extph.T[0], 10**extph.T[1], linewidth=1.5, color='#4b006e', linestyle='dashed')
ax2.plot(10.**pre.T[0], 10**pre.T[1], linewidth=2.5, color='#3AEFF4')
ax2.plot(10.**post.T[0], 10**post.T[1], linewidth=2.5, color='#52F856')
ax2.plot(10.**compton.T[0], 10**compton.T[1], linewidth=2.5, color='#0102EE')
ax2.plot(10.**total.T[0], 10**total.T[1], linewidth=1.5, color='#0B0B0B')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1	
	if ((j%10) == 0):					
		ax2.plot(10.**fr,10**(sf),linewidth=2.5,color=colors[j],linestyle='dashed')
		ax2.plot(10.**fr,10**(cf),linewidth=2.5,color=colors[j],linestyle='dashed')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_xlim([1e9,.99e26])
ax2.set_ylim([.99e-11,3.e3])
ax2.set_ylabel('F$\\nu$ (mJy)', fontsize=16)
ax2.set_xlabel('Freq (Hz)', fontsize=16)
ax2.set_yscale('log', basey=10)
ax2.set_xscale('log', basey=10)
ax2.tick_params([1.e12,1.e15,1.e18,1.e21,1.e24],fontsize=16)

fig.subplots_adjust(wspace=0)

#We create the colorbar:
sm  = plt.cm.ScalarMappable(cmap="jet")                                       # we use the color_map jet
sm.set_array(colors)                                                          # based on the array 'colors'
bar = fig.colorbar(sm,pad=0.10, ax=[ax1,ax2], ticks=[0, 0.25, 0.5, .75, 1])  # at distance pad from y axis
bar.set_label('z ($r_g$)', fontsize=18)
bar.ax.set_yticklabels(['$6.0\\times 10^0$', '$0.9\\times10^2$', '$1.3\\times10^3$', '$2.0\\times10^4$',
                        ' $3.0\\times10^5$'], fontsize = 15)
plt.show()
