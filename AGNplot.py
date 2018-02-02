import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import math

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
plt.plot(bb.T[0], 10**(bb.T[0]+bb.T[1]-3.-23.), linewidth=2.5, color='red',linestyle='dashed')
plt.plot(pre.T[0], 10**(pre.T[0]+pre.T[1]-3.-23.), linewidth=2.5, color='teal',linestyle='dashed')
plt.plot(post.T[0], 10**(post.T[0]+post.T[1]-3.-23.), linewidth=2.5, color='green',linestyle='dashed')
plt.plot(compton.T[0], 10**(compton.T[0]+compton.T[1]-3.-23.), linewidth=2.5, color='blue',linestyle='dashed')
plt.plot(total.T[0], 10**(total.T[0]+total.T[1]-3.-23.), linewidth=2.5, color='orange',linestyle='dashed')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1
	plt.plot(fr,10**(fr+sf-3.-23.),linewidth=1.5,color=colors[j])
	plt.plot(fr,10**(fr+cf-3.-23.),linewidth=1.5,color=colors[j])
plt.xlim([7.,25.])
plt.ylim([1.e-20,1.e-8])
plt.grid(True)
plt.ylabel('Flux ($\\nu$F$\\nu$) ergs s$\\rm ^{-1}$ cm$\\rm ^{-2}$', fontsize=16)
plt.xlabel('Freq (Hz)', fontsize=16)
plt.yscale('log', basey=10)

i = 0

plt.figure(2)
plt.plot(bb.T[0], 10**bb.T[1], linewidth=2.5, color='red',linestyle='dashed')
plt.plot(pre.T[0], 10**pre.T[1], linewidth=2.5, color='teal',linestyle='dashed')
plt.plot(post.T[0], 10**post.T[1], linewidth=2.5, color='green',linestyle='dashed')
plt.plot(compton.T[0], 10**compton.T[1], linewidth=2.5, color='blue',linestyle='dashed')
plt.plot(total.T[0], 10**total.T[1], linewidth=2.5, color='orange',linestyle='dashed')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1					
	plt.plot(fr,10**(sf),linewidth=1.5,color=colors[j])
	plt.plot(fr,10**(cf),linewidth=1.5,color=colors[j])
plt.xlim([7.,25.])
plt.ylim([1.e-7,4000.])
plt.grid(True)
plt.ylabel('Flux (mJy)', fontsize=16)
plt.xlabel('Freq (Hz)', fontsize=16)
plt.yscale('log', basey=10)

plt.show()
