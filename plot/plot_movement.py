import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from pylab import *

infile='movement.dat'			  #First input file
outname='movement'	  		          #Name output files will take
xlbl='Amino Acid Number'
ylbl='Amino Acid Number'
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
		
	if sys.argv[x]=='-help':
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
		mj[i].append(int(line.split()[1]))
		ol[i].append(float(line.split()[2]))
		
mi=np.array(mi)
mj=np.array(mj)
ol=np.array(ol)

fig=plt.figure(1, figsize=(11,8))
ax=fig.add_subplot(111)
cmain=ax.pcolor(mi,mj,ol,vmin=0, vmax=int(1.8*ol.mean()),cmap=plt.cm.gist_heat_r)#1.8
ax.set_title(ttl)
ax.set_xlabel(xlbl)

ax.set_xlim(mi.min(), mi.max())

ax.set_ylabel(ylbl)
ax.set_ylim(mj.max(), mj.min())

cbar=fig.colorbar(cmain,aspect=10)

fig.text(.83, .95, '$<\delta r_{ij}^2> $', horizontalalignment='center')

plt.rcParams.update({'font.size': 22})


plt.savefig(outname+'.png',format='png')
plt.show()

print 'DONE'
