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
#plt.title('$\\gamma_{\\rm max} = 15$, $\\gamma(z) \\propto z^{0.5}$, $\\theta(z) = 0.15/\\gamma(z)$, $\\sigma(z_{\\rm diss}) = 0.01$, $\\beta_{0} = 0.05$, $\\beta_{\\rm pl} = 0.0005$',fontsize=24)
plt.plot(bb.T[0], 10**(bb.T[0]+bb.T[1]-3.-23.), linewidth=2.5, color='red')
plt.plot(pre.T[0], 10**(pre.T[0]+pre.T[1]-3.-23.), linewidth=2.5, color='teal')
plt.plot(post.T[0], 10**(post.T[0]+post.T[1]-3.-23.), linewidth=2.5, color='green')
plt.plot(compton.T[0], 10**(compton.T[0]+compton.T[1]-3.-23.), linewidth=2.5, color='blue')
plt.plot(total.T[0], 10**(total.T[0]+total.T[1]-3.-23.), linewidth=2.5, color='orange')
for j in range(69):
	if (zones.T[0][i+1] - zones.T[0][i] < 0.):
		i = i + 1
	while(zones.T[0][i+1] > zones.T[0][i]):
		fr[i-299*j] = zones.T[0][i]
		sf[i-299*j] = zones.T[1][i]
		cf[i-299*j] = zones.T[2][i]
		i = i + 1
	plt.plot(fr,10**(fr+sf-3.-23.),linewidth=1.5,color=colors[j],linestyle='dashed')
	plt.plot(fr,10**(fr+cf-3.-23.),linewidth=1.5,color=colors[j],linestyle='dashed')
plt.xlim([8.5,19.])
plt.ylim([1.e-14,0.98e-10])
plt.grid(True)
plt.ylabel('Flux ($\\nu$F$\\nu$) ergs s$\\rm ^{-1}$ cm$\\rm ^{-2}$', fontsize=16)
plt.xlabel('Freq (Hz)', fontsize=16)
plt.yscale('log', basey=10)
#extend plot down to line between 49 and 50 on gedit, enough to match the = in the colors definition before plotting code
i = 0

plt.figure(2)
#plt.title('$\\gamma_{\\rm max} = 15$, $\\gamma(z) \\propto z^{0.5}$, $\\theta(z) = 0.15/\\gamma(z)$, $\\sigma(z_{\\rm diss}) = 0.01$, $\\beta_{0} = 0.4$, $\\beta_{\\rm pl} = 0.0005$',fontsize=24)
plt.plot(bb.T[0], 10**bb.T[1], linewidth=2.5, color='red')
plt.plot(pre.T[0], 10**pre.T[1], linewidth=2.5, color='teal')
plt.plot(post.T[0], 10**post.T[1], linewidth=2.5, color='green')
plt.plot(compton.T[0], 10**compton.T[1], linewidth=2.5, color='blue')
plt.plot(total.T[0], 10**total.T[1], linewidth=2.5, color='orange')
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
plt.xlim([8.5,19.])
plt.ylim([1.e-6,3.e4])
plt.grid(True)
plt.ylabel('Flux (mJy)', fontsize=16)
plt.xlabel('Freq (Hz)', fontsize=16)
plt.yscale('log', basey=10)

plt.show()
