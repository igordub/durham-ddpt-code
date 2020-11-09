import numpy as np
import os,sys
import matplotlib.pyplot as plt

# Default values
infile='fiteigen.dat'
outname='eigen_proj'
cbartitl='Time / ns'
val1=1
val2=2
ttl=''
scale=1
aset=0
bset=0

# User input values
for x in range(1,len(sys.argv)):
    if sys.argv[x] == '-i':
        infile = sys.argv[x+1]
        
    if sys.argv[x]=='-1':
        val1 = int(sys.argv[x+1])

    if sys.argv[x]=='-2':
        val2 = int(sys.argv[x+1])    

    if sys.argv[x]=='-scale':
	scale = int(sys.argv[x+1])

    if sys.argv[x]=='-cbartitle':
	cbartitl = sys.argv[x+1]
		
    if sys.argv[x]=='-min':
        amin = float(sys.argv[x+1])
        aset = 1

    if sys.argv[x]=='-max':
        amax = float(sys.argv[x+1])
        bset = 1
	
    if sys.argv[x]=='-title':
        ttl = sys.argv[x+1]

    if sys.argv[x]=='-out':
        outname = sys.argv[x+1]       
		
    if sys.argv[x]=='-help':
        print '\n\nProgram to plot eigen_proj data...\n\nOPTIONS:\n'\
            '-i = Name of input file (Default=fiteigen.dat)\n'\
            '-1 = First eigenvector number x-axis (Default=1)\n'\
            '-2 = Second eigenvector number y-axis (Default=2)\n'\
            '-scale = Devision scale for colour bar (Default=1)\n'\
            '-cbartitle = Title for colour bar (Default="Time / ns")\n'\
            '-min = scatter axis min value\n'\
            '-max = scatter axis max value\n'\
            '-out = output filename (Default=eigen_proj)\n'\
            '-title = Title for plot\n'
        exit()


xlbl='Eigenvector '+ str(val1)
ylbl='Eigenvector '+ str(val2)


inlines=open(infile,'r').readlines()

x=[]
y=[]
c=[]
d=[]

for line in inlines:
    x.append(float(line.split()[val1]))
    y.append(float(line.split()[val2]))
    d=line.split()[0]
    if d=='#':
        d=[]
    else:
        c.append(int(d))

del x[0]
del y[0]
c=np.array(c)

fig = plt.figure(1, figsize=(10,10))

from mpl_toolkits.axes_grid1 import make_axes_locatable

a=int(c.max()/10)
b=a*10

# the scatter plot:
axScatter = fig.add_subplot(111)
ax=axScatter.scatter(x, y, c=c, vmin=0,vmax=b)
axScatter.set_aspect(1.)
cbar=fig.colorbar(ax,ticks=[0, b/4, b/2, 3*b/4, b],shrink=0.2,aspect=5)
cbar.ax.set_yticklabels([0,b/(4*scale),b/(2*scale), 3*b/(4*scale), b/scale])
axScatter.set_xlabel(xlbl)
axScatter.set_ylabel(ylbl)
fig.text(.5, .85, ttl, horizontalalignment='center') 
axScatter.xaxis.grid(True,'major')
axScatter.yaxis.grid(True,'major')

# create new axes on the right and on the top of the current axes
# The first argument of the new_vertical(new_horizontal) method is
# the height (width) of the axes to be created in inches.
divider = make_axes_locatable(axScatter)
axHistx = divider.append_axes("top", 1.2, pad=0.28, sharex=axScatter)
axHisty = divider.append_axes("right", 1.2, pad=0.28, sharey=axScatter)

# make some labels invisible
plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
         visible=False)

# now determine nice limits by hand:
binwidth = 0.25
xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )


lim = ( int(xymax/binwidth) + 1) * binwidth
bins = np.arange(-lim, lim + binwidth, binwidth)
axHistx.hist(x, bins=bins)
axHisty.hist(y, bins=bins, orientation='horizontal')

# the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
# thus there is no need to manually adjust the xlim and ylim of these
# axis.

#axHistx.axis["bottom"].major_ticklabels.set_visible(False)
for tl in axHistx.get_xticklabels():
    tl.set_visible(False)
axHistx.set_yticks([0, len(x)/amax, len(x)/(amax/2)]) 
axHistx.set_ylabel('Number')
axHistx.xaxis.grid(True,'major')
axHistx.yaxis.grid(True,'major')

#axHisty.axis["left"].major_ticklabels.set_visible(False)
for tl in axHisty.get_yticklabels():
    tl.set_visible(False)
axHisty.set_xticks([0, len(y)/amax, len(y)/(amax/2)])
axHisty.set_xlabel('Number')
axHisty.xaxis.grid(True,'major')
axHisty.yaxis.grid(True,'major')


fig.text(.83, .62, cbartitl, horizontalalignment='center')
plt.rcParams.update({'font.size': 16 ,'font.family': 'serif'})

if aset==1 & bset==1:
    axScatter.set_xlim(amin, amax)
    axScatter.set_ylim(amin, amax)

plt.draw()
plt.savefig(outname+'.png',format='png')
plt.show()
