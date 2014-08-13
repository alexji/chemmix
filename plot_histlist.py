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

def plotone(envname,sfrname,postfix,options):
    k2max,k3max,cgrid = util.load_ckk(envname,sfrname)

    scale = 0.5

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
        plot_cfrac_panel(ax,histdict,scale=scale)
        h,x = histdict[(5,5)]
        h = h.astype(float)/np.sum(h)
        plot1dhist(h,x,ax=ax,color='red')
        ax.set_xlim((x[0],x[-1]))
        ax.set_ylim((0,scale*1.1))
        ax.set_yticklabels([])
    if options.save:
        plt.savefig('PLOTS/'+envname+'_'+sfrname+'_'+postfix+'_cfrac.png')
    else:
        plt.show()

def plottwo(envname,sfrname,options):
    k2max,k3max,cgrid = util.load_ckk(envname,sfrname)
    postfixlist = ['p%3.2f' % p for p in [.1,.3,.5,.7,.9]]
    scale=0.5
    
    fig,axarr = plt.subplots(len(postfixlist),k3max,sharex=True,sharey=True,figsize=(20,8))
    fig.subplots_adjust(hspace=0,wspace=0)
    for irow,postfix in enumerate(postfixlist):
        histlist = util.load_histlist(envname,sfrname,postfix)
        for icol in range(k3max):
            ilist = icol+1
            histdict = histlist[ilist]
            ax = axarr[irow,icol]
            if irow==0: ax.set_title(r'$k_{III}='+str(ilist)+r'$')
            if icol==0: ax.set_ylabel(postfixlist[irow])

            plot_cfrac_panel(ax,histdict,scale=scale)
            h,x = histdict[(5,5)]
            h = h.astype(float)/np.sum(h)
            plot1dhist(h,x,ax=ax,color='red')
            ax.set_xlim((x[0],x[-1]))
            ax.set_ylim((0,scale*1.1))
            ax.set_yticklabels([])
    if options.save:
        plt.savefig('PLOTS/'+envname+'_'+sfrname+'_pgrid.png')
    else:
        plt.show()
            
def plotthree(envname,sfrname,options):
    k2max,k3max,cgrid = util.load_ckk(envname,sfrname)
    postfixlist = ['p%3.2f' % p for p in [.1,.3,.5,.7,.9]]
    scale=0.5
    
    fig,axarr = plt.subplots(len(postfixlist),k3max,sharex=True,sharey=True,figsize=(20,8))
    fig.subplots_adjust(hspace=0,wspace=0)
    for irow,postfix in enumerate(postfixlist):
        histlist = util.load_histlist(envname,sfrname,postfix)
        for icol in range(k3max):
            ilist = icol+1
            histdict = histlist[ilist]
            ax = axarr[irow,icol]
            if irow==0: ax.set_title(r'$k_{III}='+str(ilist)+r'$')
            if icol==0: ax.set_ylabel(postfixlist[irow])

            ax.plot([-4,-3],[scale*0.75,scale*0.3],'ko-')
            plot_cfrac_panel(ax,histdict,scale=scale)
            h,x = histdict[(5,5)]
            h = h.astype(float)/np.sum(h)
            plot1dhist(h,x,ax=ax,color='red')
            ax.set_xlim((-6,-2))
            ax.set_ylim((0,scale*1.1))
            ax.set_yticklabels([])
    if options.save:
        plt.savefig('PLOTS/'+envname+'_'+sfrname+'_pgridzoom.png')
    else:
        plt.show()
            
def plotfour(envname,sfrname,options):
    k2max,k3max,cgrid = util.load_ckk(envname,sfrname)
    postfixlist = ['p%3.2f' % p for p in [.1,.3,.5,.7,.9]]
    scale = 0.5

    carr = cgrid[0,1:]
    print "Fraction only Pop III:", np.sum(carr)

    carr = carr/np.sum(carr)

    fig,axarr = plt.subplots(len(postfixlist),1,sharex=True,sharey=True,figsize=(4,8))
    fig.subplots_adjust(hspace=0,wspace=0)
    for irow,postfix in enumerate(postfixlist):
        ax = axarr[irow]
        histlist = util.load_histlist(envname,sfrname,postfix)
        Htot = 0.
        htot = 0.
        for icol in range(k3max):
            histdict = histlist[icol+1]
            ax.set_ylabel(postfixlist[irow])
            H,edgesC,edgesFe = histdict[(5,0)] #Fe,C hist
            Htot = Htot + H*carr[icol]
            h,x = histdict[(5,5)] #Fe hist
            htot = htot + h*carr[icol]
        histdict = {(5,0): (H,edgesC,edgesFe)}
        plot_cfrac_panel(ax,histdict,scale=scale)
        plot1dhist(htot.astype(float)/np.sum(htot),x,ax=ax,color='red')
        ax.plot([-4,-3],[scale*0.75,scale*0.3],'ko-')
        ax.set_xlim((-6,-2))
        ax.set_ylim((0,scale*1.1))
        ax.set_yticklabels([])
    if options.save:
        plt.savefig('PLOTS/'+envname+'_'+sfrname+'_pweightzoom.png')
    else:
        plt.show()
        
#
#            ax.plot([-4,-3],[scale*0.75,scale*0.3],'ko-')
#            plot_cfrac_panel(ax,histdict,scale=scale)
#            h,x = histdict[(5,5)]
#            h = h.astype(float)/np.sum(h)
#            plot1dhist(h,x,ax=ax,color='red')
#            ax.set_xlim((-6,-2))
#            ax.set_ylim((0,scale*1.1))
#            ax.set_yticklabels([])
#    if options.save:
#        plt.savefig('PLOTS/'+envname+'_'+sfrname+'_pgridzoom.png')
#    else:
#        plt.show()


if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--save',action='store_true',dest='save',default=False)
    parser.add_option('--plotone',action='store_true',dest='plotone',default=False)
    parser.add_option('--plottwo',action='store_true',dest='plottwo',default=False)
    parser.add_option('--plotthree',action='store_true',dest='plotthree',default=False)
    parser.add_option('--plotfour',action='store_true',dest='plotfour',default=False)
    options,args = parser.parse_args()

    if options.plotone:
        envname,sfrname,postfix = args
        plotone(envname,sfrname,postfix,options)
    elif options.plottwo:
        envname,sfrname = args
        plottwo(envname,sfrname,options)
    elif options.plotthree:
        envname,sfrname = args
        plotthree(envname,sfrname,options)
    elif options.plotfour:
        envname,sfrname = args
        plotfour(envname,sfrname,options)
