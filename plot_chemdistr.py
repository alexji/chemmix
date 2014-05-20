import numpy as np
import pylab as plt
from matplotlib.ticker import NullFormatter
import h5py
from optparse import OptionParser

import karlsson
from tophat import TopHat

def plot_2d_hist(xdat,ydat,
                 nxbins=50,nybins=50,
                 xlabel=None,ylabel=None,
                 return_data=False):
    """
    Plots 2D histogram of a data array, using columns given by colx and coly
    Based on http://jessresearch.blogspot.com/2012/04/pretty-plots-2d-histogram-with-1d.html
    """

    xmin = np.min(xdat); xmax = np.max(xdat)
    ymin = np.min(ydat); ymax = np.max(ydat)

    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    # Set up the geometry of the three plots
    rect_histxy = [left, bottom, width, height] # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
    # Set up the size of the figure
    fig = plt.figure(figsize=(9.5,9))
    # Make the three axes
    axHistxy = plt.axes(rect_histxy) # temperature plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram
    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHistx.yaxis.set_major_formatter(nullfmt)
    axHisty.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Plot 2D histogram
    xbins = np.linspace(xmin,xmax,nxbins)
    ybins = np.linspace(ymin,ymax,nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    #aspectratio = 1.0*(xmax)/(1.0*ymax)
    aspectratio='auto'
    H, xedges,yedges = np.histogram2d(ydat,xdat,bins=(ybins,xbins))
    X = xcenter; Y = ycenter
    axHistxy.imshow(H, extent = [xmin,xmax,ymin,ymax],
                    interpolation='nearest',
                    origin='lower',aspect=aspectratio)
    axHistxy.set_xlabel(xlabel)
    axHistxy.set_ylabel(ylabel)

    # Plot 1D histograms
    axHistx.hist(xdat,bins=xbins,color='blue')
    axHisty.hist(ydat,bins=ybins,orientation='horizontal',color='red')
    axHistx.set_xlim(axHistxy.get_xlim())
    axHisty.set_ylim(axHistxy.get_ylim())

    if return_data:
        return fig,H
    return fig


if __name__=="__main__":
### PICK ONE!
#    filename = 'CHEMGRIDS/minihalo_chemgrid_N06i.hdf5'
#    plotprefix = 'minihalo_N06i'
#    Mhalo,vturb,lturb,nSN,trecovery = karlsson.params_minihalo()
    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_N06i.hdf5'
    plotprefix = 'atomiccoolinghalo_N06i'
    Mhalo,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalo()

    vturb *= 3.16/3.08 * .001 #km/s to kpc/yr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    print filename
    print "Mhalo",Mhalo,"vturb",vturb,"lturb",lturb
    print "Dt",Dt,"uSN",uSN
    th = TopHat(Mhalo=Mhalo,nvir=0.1,fb=0.1551,mu=1.4)
    RHO = th.get_rho_of_t_fn()
    VMIX = karlsson.get_Vmixfn_K08(RHO,Dt=Dt)
    WISM = karlsson.wISM_K05
    PSI = lambda t: uSN
    MUFN = karlsson.get_mufn(VMIX,PSI)

    parser = OptionParser()
    parser.add_option('-k','--kmax',
                      action='store',type='int',dest='kmax',
                      default=10,
                      help='pick kmax')
    parser.add_option('--show',
                      action='store_true',dest='showplots',
                      default=False,
                      help='Call plt.show() in addition to saving all the figures')
    options,args = parser.parse_args()


    f = h5py.File(filename,'r')
    chemarr = np.array(f['chemgrid'])
    f.close()

    kmax = options.kmax
    print "kmax to plot:",kmax
    Nstar,numyields,kmax_tmp = chemarr.shape
    assert kmax<=kmax_tmp
    chemarr = chemarr[:,:,0:kmax]
    
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore') #ck gives a convergence warning which isn't important
        ck = karlsson.calc_ck(kmax,WISM,MUFN,PSI)

    for k in range(kmax):
        chemarr[:,:,k] *= ck[k]
    chemarr = np.log10(np.sum(chemarr,2))
    
    elemnames = ['C', 'O', 'Mg', 'Si', 'Ca', 'Fe']
    asplund09solar = np.array([8.43,8.69,7.60,7.51,6.34,7.50]) - 12 #log10(Z/H) solar
    for z in range(numyields):
        chemarr[:,z] -= asplund09solar[z]
        print elemnames[z],np.min(chemarr[:,z]),np.max(chemarr[:,z])
    
    ## BEGIN PLOTS
    plt.figure()
    plt.plot(np.arange(kmax)+1,ck.cumsum(),'ks-')
    plt.xlabel('k'); plt.ylabel(r'cumulative $c_k$')
    plt.xlim(0,kmax+1); plt.ylim((0,1))
    plt.xticks(np.arange(kmax)+1)
    plt.savefig(plotprefix+'_k'+str(kmax)+'_ck.png')

    plt.figure(figsize=(7.5,8.5))
    plt.subplots_adjust(hspace=0.25)
    #bins = np.arange(-4.6,-1.05,.05)
    bins = np.arange(-6.5,0.05,.05)
    for i in range(6):
        plt.subplot(3,2,i+1)
        plt.hist(chemarr[:,i],bins=bins,normed=True)
        plt.xticks([-6,-5,-4,-3,-2,-1,0])
        plt.xlabel('['+elemnames[i]+'/H]')
        plt.yticks([0,1,2,3,4]); plt.ylim((0,4))
    plt.savefig(plotprefix+'_k'+str(kmax)+'_sixelem.png',bbox_inches='tight')

    plot_2d_hist(chemarr[:,5],chemarr[:,0]-chemarr[:,5],
                 nxbins=100,nybins=100,
                 xlabel='[Fe/H]',ylabel='[C/Fe]')
    plt.savefig(plotprefix+'_k'+str(kmax)+'_CFe-FeH.png',bbox_inches='tight')

#    plot_2d_hist(chemarr[:,5],chemarr[:,2]-chemarr[:,5],
#                 nxbins=100,nybins=100,
#                 xlabel='[Fe/H]',ylabel='[Mg/Fe]')
#    plt.savefig(plotprefix+'_k'+str(kmax)+'_nomoto06_K05_MgFe-FeH.png',bbox_inches='tight')
#
#    plot_2d_hist(chemarr[:,5],chemarr[:,3],
#                 nxbins=100,nybins=100,
#                 xlabel='[Fe/H]',ylabel='[Si/H]')
#    plt.savefig(plotprefix+'_k'+str(kmax)+'_nomoto06_K05_SiH-FeH.png',bbox_inches='tight')
#
    if options.showplots:
        plt.show()
