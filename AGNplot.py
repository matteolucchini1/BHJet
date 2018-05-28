import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import math
import matplotlib.ticker as plticker

total = np.genfromtxt("outputs/total.dat")
bb = np.genfromtxt("outputs/bb.dat")
pre = np.genfromtxt("outputs/presyn.dat")
post = np.genfromtxt("outputs/postsyn.dat")
compton = np.genfromtxt("outputs/com.dat")

zones = np.genfromtxt("outputs/zones.dat")

fr = np.zeros(298)
sf = np.zeros(298)
cf = np.zeros(298)

i = 0
n = 70
colors = pl.cm.jet(np.linspace(0.25,1,n))

plt.figure(1)
#plt.title('$\\gamma_{\\rm max} = 5$, $z_{\\rm acc} = z_{\\rm sh} = 630$, $f_{\\rm heat} = 10$, $\\sigma(z_{\\rm diss}) = 0.024$, $N_{\\rm j} = 5*10^{-3}$',fontsize=24)
plt.plot(10.**bb.T[0], 10**(bb.T[0]+bb.T[1]-3.-23.), linewidth=2.5, color='#FB02ED', linestyle='dashed')
plt.plot(10.**pre.T[0], 10**(pre.T[0]+pre.T[1]-3.-23.), linewidth=2.5, color='#3AEFF4',zorder=151)
plt.plot(10.**post.T[0], 10**(post.T[0]+post.T[1]-3.-23.), linewidth=2.5, color='#52F856')
plt.plot(10.**compton.T[0], 10**(compton.T[0]+compton.T[1]-3.-23.), linewidth=2.5, color='#0102EE')
plt.plot(10.**total.T[0], 10**(total.T[0]+total.T[1]-3.-23.), linewidth=1.5, color='#0B0B0B')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1
	if ((j%5) == 0):
		#print("xfig sucks")
		plt.plot(10.**fr,10**(fr+sf-3.-23.),linewidth=2.5,color=colors[j],linestyle='dashed',zorder=150-j)
		plt.plot(10.**fr,10**(fr+cf-3.-23.),linewidth=2.5,color=colors[j],linestyle='dashed',zorder=150-j)
plt.xlim([1e9,1e26])
plt.ylim([1.e-14,4.0e-10])
plt.ylabel('Flux ($\\nu$F$\\nu$) ergs s$\\rm ^{-1}$ cm$\\rm ^{-2}$', fontsize=18)
plt.xlabel('Frequency (Hz)', fontsize=18)
plt.yscale('log', basey=10)
plt.xscale('log', basex=10)
plt.xticks([1.e12,1.e15,1.e18,1.e21,1.e24],fontsize=16)
plt.yticks([1.e-14,1.e-13,1.e-12,1.e-11,1.e-10],fontsize=16)

i = 0

plt.figure(2)
#plt.title('$\\gamma_{\\rm max} = 15$, $z_{\\rm acc} = 2.5*10^{5}$, $z_{\\rm sh} = 25$, $f_{\\rm heat} = 1$, $\\sigma(z_{\\rm diss}) = 1$, $N_{\\rm j} = 1*10^{-2}$',fontsize=24)
plt.plot(bb.T[0], 10**bb.T[1], linewidth=1.5, color='#FB02ED', linestyle='dashed')
plt.plot(pre.T[0], 10**pre.T[1], linewidth=2.5, color='#3AEFF4')
plt.plot(post.T[0], 10**post.T[1], linewidth=2.5, color='#52F856')
plt.plot(compton.T[0], 10**compton.T[1], linewidth=2.5, color='#0102EE')
plt.plot(total.T[0], 10**total.T[1], linewidth=1.0, color='#0B0B0B')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1					
	plt.plot(fr,10**(sf),linewidth=1.5,color=colors[j],linestyle='dashed')
	plt.plot(fr,10**(cf),linewidth=1.5,color=colors[j],linestyle='dashed')
plt.xlim([8.9,20.])
plt.ylim([1.e-6,3.e2])
plt.ylabel('Flux (mJy)', fontsize=16)
plt.xlabel('Freq (Hz)', fontsize=16)
plt.yscale('log', basey=10)

plt.show()
