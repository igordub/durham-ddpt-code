import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
from pylab import *

infile='summary.dat'				#First input file
outname='overlap'			#Name output files will take
xlbl='Eigenvector 1'
ylbl='Eigenvector 2'
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

for line in inlines:
	if line=='\n' or line.startswith('#'):
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

fig=plt.figure(1, figsize=(11,9))
ax=fig.add_subplot(111,autoscale_on=False)
cmain=ax.pcolor((mi-0.5),(mj-0.5),ol,vmin=0, vmax=0.8,cmap=plt.cm.gist_yarg)#cmap=plt.cm.YlOrRd)
ax.set_aspect(1)

ax.set_title(ttl)

ax.set_xlabel(xlbl)
ax.set_xlim(mi.min()-0.5, mi.max()-0.5)

ax.set_ylabel(ylbl)
ax.set_ylim(mj.min()-0.5, mj.max()-0.5)

cbar=fig.colorbar(cmain,aspect=10,shrink=0.9,ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

fig.text(.83, .92, 'Overlap', horizontalalignment='center')

plt.rcParams.update({'font.size': 22})

plt.savefig(outname+'.png',format='png')
#plt.show()

print 'DONE'
