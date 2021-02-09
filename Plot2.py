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
rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rcParams.update({'font.size': 20})

pars = np.genfromtxt("Input/ip.dat")
data = np.genfromtxt("1H0323342.dat")


kevconv = 1.
#kevconv = 2.41*10**17

#fluxconv = 1.
fluxconv = 4.*math.pi*(pars[2]*3.*10**21)**2

if (pars[0] >= 1.e3):
	ulim_f = 1.e27
	blim_f = 0.5e9
	ulim_fl = 1.e-9
	blim_fl = 1.e-14
	ulim_fd = 1.e5
	blim_fd = 1.e-13
else:
	ulim_f = 0.9e17
	blim_f = 0.5e9
	ulim_fl = 1e-9
	blim_fl = 1.e-16
	ulim_fd = 30.
	blim_fd = 1.01e-7
	
zmin = 2.0
zmax = pars[9]

zheight=np.zeros(5)

zheight=np.zeros(5) # the array that we use to print the labels in colorbar
for n in range(0,len(zheight)):
    zheight[n] = np.log10(np.int((pow(10.,math.log10(zmax/zmin)*(n*.25) + math.log10(zmin)))))

zheight = np.around(zheight,decimals=2)

z = pars[3]
plotcheck = pars[26]

if (plotcheck >= 1):
	Presyn = np.genfromtxt("Output/Presyn.dat")
	Postsyn = np.genfromtxt("Output/Postsyn.dat")
	Precom = np.genfromtxt("Output/Precom.dat")
	Postcom = np.genfromtxt("Output/Postcom.dat")
	Disk = np.genfromtxt("Output/Disk.dat")
	Corona = np.genfromtxt("Output/Corona.dat")
	BB = np.genfromtxt("Output/BB.dat")
	Total = np.genfromtxt("Output/Total.dat")

mjy = 1.e-26

i = 0
j = 0
nzones = 80

colors = pl.cm.magma(np.linspace(0.,0.9,nzones))

size_cyclo_arr = np.zeros(nzones)
size_com_arr = np.zeros(nzones)

if (plotcheck >=2):
	Cyclosyn_zones = np.genfromtxt("Output/Cyclosyn_zones.dat")
	Compton_zones = np.genfromtxt("Output/Compton_zones.dat")
	Numdens = np.genfromtxt("Output/Numdens.dat")
	size3 = len(Numdens.T[0])/nzones
	i = 0
	j = 0
	k = 0
	for j in range(len(Cyclosyn_zones.T[0])-1):
		if (Cyclosyn_zones.T[0][j+1]>Cyclosyn_zones.T[0][j]):
			k = k + 1
		else:
			size_cyclo_arr[i] = k + 1
			i = i + 1
			k = 0		
	i = 0
	j = 0
	k = 0
	for j in range(len(Compton_zones.T[0])-1):
		if (Compton_zones.T[0][j+1]>Compton_zones.T[0][j]):
			k = k + 1
		else:
			size_com_arr[i] = k + 1
			i = i + 1
			k = 0
	g = np.zeros(size3)
	ng = np.zeros(size3)
	p = np.zeros(size3)
	nparr = np.zeros(size3)

i = 0
j = 0
l = 0
totindex1 = 0
totindex2 = 0
fig3, (ax1) = plt.subplots(1,1,figsize=(9.,6))
for i in range(nzones-1):
	nu_cyclosyn = np.zeros(int(size_cyclo_arr[i]))
	lnu_cyclosyn = np.zeros(int(size_cyclo_arr[i]))	
	nu_compton = np.zeros(int(size_com_arr[i]))
	lnu_compton = np.zeros(int(size_com_arr[i]))	
	for j in range(int(size_cyclo_arr[i])):
		nu_cyclosyn[j] = Cyclosyn_zones.T[0][totindex1+j]/kevconv
		lnu_cyclosyn[j] = Cyclosyn_zones.T[1][totindex1+j]*nu_cyclosyn[j]*mjy*kevconv
	for l in range(int(size_com_arr[i])):	
		nu_compton[l] = Compton_zones.T[0][totindex2+l]/kevconv
		lnu_compton[l] = Compton_zones.T[1][totindex2+l]*nu_compton[l]*mjy*kevconv
	totindex1 = totindex1 + int(size_cyclo_arr[i])
	totindex2 = totindex2 + int(size_com_arr[i])	
	if ((i%1) == 0):
		ax1.plot(nu_cyclosyn,lnu_cyclosyn*fluxconv,linewidth=2.5,color=colors[i],zorder=nzones-i,linestyle='dashed')
		ax1.plot(nu_compton,lnu_compton*fluxconv,linewidth=2.5,color=colors[i],zorder=nzones-i,linestyle='dashed')
ax1.plot(Presyn.T[0]/kevconv,Presyn.T[1]*Presyn.T[0]*mjy*fluxconv,linewidth=2.5,color='cyan',zorder=nzones+1)
ax1.plot(Postsyn.T[0]/kevconv,Postsyn.T[1]*Postsyn.T[0]*mjy*fluxconv,linewidth=2.5,color='green',zorder=nzones+1)
ax1.plot(Precom.T[0]/kevconv,Precom.T[1]*Precom.T[0]*mjy*fluxconv,linewidth=2.5,color='blue',zorder=nzones+1)
ax1.plot(Postcom.T[0]/kevconv,Postcom.T[1]*Postcom.T[0]*mjy*fluxconv,linewidth=2.5,color='purple',zorder=nzones+1)
#ax1.plot(Corona.T[0]/kevconv,Corona.T[1]*Corona.T[0]*mjy*fluxconv,linewidth=2.5,color='#0165fc',zorder=nzones+1)
ax1.plot(Disk.T[0]/kevconv,Disk.T[1]*Disk.T[0]*mjy*fluxconv,linewidth=2.5,color='red',zorder=nzones+1)
ax1.plot(BB.T[0]/kevconv,BB.T[1]*BB.T[0]*mjy*fluxconv,linewidth=2.5,color='orange',zorder=nzones+1)
ax1.plot(Total.T[0]/kevconv,Total.T[1]*Total.T[0]*mjy*fluxconv,linewidth=1.5,color='grey',zorder=nzones+1)
ax1.errorbar(data.T[0]/kevconv,data.T[2],yerr=data.T[3],color='black',fmt='o')
ax1.set_ylim([blim_fl*fluxconv,ulim_fl*fluxconv])
ax1.set_xlim([blim_f/kevconv,ulim_f/kevconv])
ax1.set_xscale('log', basex=10)
ax1.set_yscale('log', basex=10)
if (kevconv !=1):
	ax1.set_xlabel('Energy (kev)',fontsize=18)
else:
	ax1.set_xlabel('Frequency (Hz)',fontsize=18)
if (fluxconv !=1):
	ax1.set_ylabel('Luminosity (erg/s)',fontsize=18)
else:
	ax1.set_ylabel('Flux (erg/s/cm^2)',fontsize=18)

i = 0
j = 0
l = 0
totindex1 = 0
totindex2 = 0
fig4, (ax2) = plt.subplots(1,1,figsize=(9.,6))
for i in range(nzones-1):
	nu_cyclosyn = np.zeros(int(size_cyclo_arr[i]))
	lnu_cyclosyn = np.zeros(int(size_cyclo_arr[i]))
	nu_compton = np.zeros(int(size_com_arr[i]))
	lnu_compton = np.zeros(int(size_com_arr[i]))
	for j in range(int(size_cyclo_arr[i])):
		nu_cyclosyn[j] = Cyclosyn_zones.T[0][totindex1+j]/kevconv
		lnu_cyclosyn[j] = Cyclosyn_zones.T[1][totindex1+j]
	for l in range(int(size_com_arr[i])):	
		nu_compton[l] = Compton_zones.T[0][totindex2+l]/kevconv
		lnu_compton[l] = Compton_zones.T[1][totindex2+l]
	totindex1 = totindex1 + int(size_cyclo_arr[i])
	totindex2 = totindex2 + int(size_com_arr[i])	
	if ((i%5)== 0):
		ax2.plot(nu_cyclosyn,lnu_cyclosyn,linewidth=2.5,color=colors[i],zorder=nzones-i)#,linestyle='dashed')
		#ax2.plot(nu_compton,lnu_compton,linewidth=2.5,color=colors[i],zorder=nzones-i,linestyle='dashed')
ax2.plot(Presyn.T[0]/kevconv,Presyn.T[1],linewidth=3.5,color='black',zorder=nzones+1,linestyle='dashed')
#ax2.plot(Postsyn.T[0]/kevconv,Postsyn.T[1],linewidth=2.5,color='green',zorder=nzones+1)
#ax2.plot(Precom.T[0]/kevconv,Precom.T[1],linewidth=2.5,color='blue',zorder=nzones+1)
#ax2.plot(Postcom.T[0]/kevconv,Postcom.T[1],linewidth=2.5,color='purple',zorder=nzones+1,linestyle='dotted')
#ax2.plot(Corona.T[0]/kevconv,Corona.T[1],linewidth=2.5,color='#0165fc',zorder=nzones+1)
#ax2.plot(Disk.T[0]/kevconv,Disk.T[1],linewidth=1.5,color='red',zorder=nzones+1)
#ax2.plot(Total.T[0]/kevconv,Total.T[1],linewidth=1.5,color='black',zorder=nzones+1)
ax2.set_ylim([blim_fd,ulim_fd])
ax2.set_xlim([blim_f/kevconv,ulim_f/kevconv])
ax2.set_xscale('log', basex=10)
ax2.set_yscale('log', basex=10)
#ax2.yaxis.tick_right()
#ax2.yaxis.set_label_position("right")
if (kevconv !=1):
	ax2.set_xlabel('Energy (kev)',fontsize=24)
else:
	ax2.set_xlabel('Frequency (Hz)',fontsize=24)
ax2.set_ylabel('Flux density (arb. units)',fontsize=24)
ax2.tick_params(axis='both', which='major', labelsize=20)

if(plotcheck>=2): #We create the colorbar:
	divider = make_axes_locatable(plt.gca())
	cax = divider.append_axes("right", "5%", pad="3%")
	sm = plt.cm.ScalarMappable(cmap="magma")                                
	sm.set_array(colors)                                                
	bar = fig.colorbar(sm,cax=cax,ticks=[0.01,0.25,0.5,.75,0.99])  
	bar.set_label('$log_{10}(z)\,(r_g$)', fontsize=18)
	bar.ax.set_yticklabels(zheight,fontsize = 15)

plt.tight_layout()

i = 0
j = 0
if(plotcheck>=2):
	plt.figure(2)
	fig5, (ax3) = plt.subplots(1,1,figsize=(9.,6))	
	for j in range(nzones-1):
		for i in range(size3):
			p[i] = Numdens.T[0][i+size3*j]
			nparr[i] = Numdens.T[2][i+size3*j]
		if ((j%10) == 0): #this plots one in five zones for clarity
			ax3.plot(p,nparr*p,linewidth=3.5,color=colors[j],zorder=150-j)
	ax3.set_xlabel('Momentum (erg*s/cm)',fontsize=22)
	ax3.set_ylabel('Number density ($\\#/cm^3$)',fontsize=22)
	ax3.tick_params(axis='both', which='major', labelsize=20)
	ax3.set_xscale('log', basex=10)
	ax3.set_yscale('log', basex=10)

	i = 0
	j = 0
	plt.figure(3)
	fig6, (ax4) = plt.subplots(1,1,figsize=(9.,6))	
	for j in range(nzones-1):
		for i in range(size3):
			g[i] = Numdens.T[1][i+size3*j]
			ng[i] = Numdens.T[3][i+size3*j]
		if ((j%3) == 0): #this plots one in five zones for clarity
			ax4.plot(g,ng*g,linewidth=3.5,color=colors[j],zorder=150-j)
	ax4.set_xlabel('Electron $\\gamma$',fontsize=18)
	ax4.set_ylabel('Number density ($\\#/cm^3$)',fontsize=18)	
	ax4.set_xscale('log', basex=10)
	ax4.set_yscale('log', basex=10)
	ax4.yaxis.tick_right()
	ax4.yaxis.set_label_position("right")
	
	fig2.subplots_adjust(wspace=0)
	#We create the colorbar:
	sm2 = plt.cm.ScalarMappable(cmap="magma")                                
	sm2.set_array(colors)                                                
	bar2 = fig2.colorbar(sm2,pad=0.10,ax=[ax3],ticks=[0.01,0.25,0.5,.75,0.99])  
	bar2.set_label('$log_{10}(z)\,(r_g$)', fontsize=18)
	bar2.ax.set_yticklabels(zheight,fontsize = 15)

plt.show()
