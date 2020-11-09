import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from pylab import *

infile='3Dplot.dat'			  #First input file
outname='var'	  		          #Name output files will take
xlbl='Amino Acid Number'
ylbl='$k_R/k$'
ttl=''

mi=[]
mj=[]
ol=[]
i=-1

#############################################################################
# Read arguments from terminal, and assign input files and a name that all output files will contain. 
#############################################################################
for x in range(1,len(sys.argv)):
	if sys.argv[x] == '-i':
		infile = sys.argv[x+1]
		
	if sys.argv[x]=='-xlabel':
		xlbl = sys.argv[x+1]
		
	if sys.argv[x]=='-ylabel':
		ylbl = sys.argv[x+1]
	
	if sys.argv[x]=='-title':
		ttl = sys.argv[x+1]
		
	if sys.argv[x]=='-h':
		print '\n\nProgram to plot overlap data...\n\nOPTIONS:\n'\
		'-i = Name of input file (Default=overlap.dat)\n'\
		'-xlabel = Label for x axis (Default=mode i)\n'\
		'-ylabel = Label for y axis (Default=mode j)\n'\
		'-title = Title for plot\n'
		exit()
		
inlines=open(infile,'r').readlines()

if inlines[-1]=='\n':
	inlines[-1:]=[]

i=i+1
mi.append([])
mj.append([])
ol.append([])

for line in inlines:
	if line=='\n':
		i=i+1
		mi.append([])
		mj.append([])
		ol.append([])
		
	else:
		mi[i].append(int(line.split()[0]))
		mj[i].append(float(line.split()[1]))
		ol[i].append(float(line.split()[2]))
		
		if float(line.split()[1])==1.0:
			mid=float(line.split()[2])
		
mi=np.array(mi)
mj=np.array(mj)
ol=np.array(ol)

valmax=ol.max()
valmin=ol.min()

space=min(valmax-mid,mid-valmin)

space=float(int(space*1000))/1000
mid=float(int(mid*1000))/1000

fig=plt.figure(1, figsize=(15,8))
ax=fig.add_subplot(111)
cmain=ax.pcolor(mi-0.5,log(mj),ol,vmin=mid-space, vmax=mid+space,cmap=plt.cm.RdBu_r)
ax.set_title(ttl)
ax.set_xlabel(xlbl)
#ax.set_xticks([20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
#ax.set_xticklabels(['\n 20', '\n 40', '\n 60', '\n 80', '\n 100', '\n 120', '\n 140', '\n 160', '\n 180', '\n 200'])

ax.set_xlim(mi.min(), mi.max())
ax.set_ylabel(ylbl)

ax.set_ylim(log(mj.min()), log(mj.max()))
ax.set_yticks([log(0.25), log(0.5), log(0.75), log(1), log(2), log(3), log(4)])
ax.set_yticklabels(['0.25', '0.50', '0.75', '1.00', '2.00', '3.00', '4.00'])

cbar=fig.colorbar(cmain,aspect=10)#,ticks=[0.95,1,1.05,1.10,1.15,1.2,1.25,1.3])
#cbar.ax.set_yticklabels(['$0.95$','$1.00$','$1.05$','$1.10$','$1.15$','$1.20$','$1.25$','$1.30$'])

fig.text(.82, .95, '$K_2$/$K_1$', horizontalalignment='center')

plt.rcParams.update({'font.size': 22,'font.family': 'serif'})


plt.savefig(outname+'.pdf',format='pdf')
plt.show()

print 'DONE'
