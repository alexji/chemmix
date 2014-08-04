import numpy as np
import pylab as plt
from optparse import OptionParser
from plot_util import density_contour,plot1dhist
import util

def plot_cfrac_panel(ax,histdict,CFecrit=0.75,scale=1.0):
    iC = 0; iFe = 5
    H,edgesC,edgesFe = histdict[(iFe,iC)]
    pdf = util.hist2d2pdf(H,edgesC,edgesFe)
    H = H.T #put Fe on the X axis
    midFe = (edgesFe[1:]+edgesFe[:-1])/2.
    midC  = (edgesC[1:]+edgesC[:-1])/2.
    cfrac = np.zeros(len(midFe))
    for i,Fe in enumerate(midFe):
        iiCrich = midC-Fe > CFecrit
        cfrac[i] = np.sum(H[i,iiCrich])/np.sum(H[i,:])
    ax.plot(edgesFe[:-1],cfrac*scale,'b-',drawstyle='steps')

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--save',action='store_true',dest='save',default=False)
    options,args = parser.parse_args()
    envname,sfrname,postfix = args
    
    k2max,k3max,cgrid = util.load_ckk(envname,sfrname)

    fig,axarr = plt.subplots(k2max+1,k3max+1,sharex=True,sharey=True,figsize=(20,20))
    fig.subplots_adjust(hspace=0,wspace=0)
    histlist = util.load_histlist(envname,sfrname,postfix)
    for tasknum,histdict in enumerate(histlist):
        kII,kIII = divmod(tasknum,k3max+1)
        ax = axarr[kII,kIII]
        if kII==0: ax.set_title(r'$k_{III}='+str(kIII)+r'$')
        if kIII==0: ax.set_ylabel(r'$k_{II}='+str(kII)+r'$')
        if tasknum==0:
            continue
        plot_cfrac_panel(ax,histdict)
        h,x = histdict[(5,5)]
        h = h.astype(float)/np.sum(h)
        plot1dhist(h,x,ax=ax,color='red')
        ax.set_xlim((x[0],x[-1]))
        ax.set_ylim((0,1))
    if options.save:
        plt.savefig('PLOTS/'+envname+'_'+sfrname+'_'+postfix+'_cfrac.png')
    else:
        plt.show()
