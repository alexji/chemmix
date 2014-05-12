import numpy as np
import pylab as plt
from matplotlib.ticker import NullFormatter

from karlsson import *

from optparse import OptionParser

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
    #tmin = 0.0; tmax = 1000.0; dt = 0.03
    tmin = 0.0; tmax = 1000.0; dt = 0.03
    plotprefix = 'dt0.03'
    plotprefix = 'dt0.03_Cenhanced_p0.75'
    #plotprefix = 'dt0.03_Cenhanced_p0.5'
    #plotprefix = 'dt0.03_Cenhanced_p0.25'
    Nstars = 10**5
    print "tmax,dt:",tmax,dt
    print "Nstars:",Nstars
    fileprefix="Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    #chemarr = np.load(fileprefix+'_chemgrid_N06_N'+str(Nstars)+'.npy')
    chemarr = np.load(fileprefix+'_chemgrid_N06_Cenhanced_p0.75_N'+str(Nstars)+'.npy')
    #chemarr = np.load(fileprefix+'_chemgrid_N06_Cenhanced_p0.5_N'+str(Nstars)+'.npy')
    #chemarr = np.load(fileprefix+'_chemgrid_N06_Cenhanced_p0.25_N'+str(Nstars)+'.npy')

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

    kmax = options.kmax
    print "KMAX:",kmax
    Nstar,numyields,kmax_tmp = chemarr.shape
    assert kmax<=kmax_tmp
    chemarr = chemarr[:,:,0:kmax]
    
    aLMS = 0.835
    #mufn = get_mufn_K05A()
    #psi_K05 = uSN_K05A
    #ck = calc_ck(kmax,wISM_K05,mufn,psi_K05)
    mufn = get_mufn_K05A()
    psi_K05 = uSN_K05A
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore') #ck gives a convergence warning which isn't important
        ck = calc_ck(kmax,wISM_K05,mufn,psi_K05)

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
    plt.savefig(plotprefix+'_k'+str(kmax)+'_nomoto06_K05_sixelem.png',bbox_inches='tight')

    plot_2d_hist(chemarr[:,5],chemarr[:,0]-chemarr[:,5],
                 nxbins=100,nybins=100,
                 xlabel='[Fe/H]',ylabel='[C/Fe]')
    plt.savefig(plotprefix+'_k'+str(kmax)+'_nomoto06_K05_CFe-FeH.png',bbox_inches='tight')

    plot_2d_hist(chemarr[:,5],chemarr[:,2]-chemarr[:,5],
                 nxbins=100,nybins=100,
                 xlabel='[Fe/H]',ylabel='[Mg/Fe]')
    plt.savefig(plotprefix+'_k'+str(kmax)+'_nomoto06_K05_MgFe-FeH.png',bbox_inches='tight')

    plot_2d_hist(chemarr[:,5],chemarr[:,3],
                 nxbins=100,nybins=100,
                 xlabel='[Fe/H]',ylabel='[Si/H]')
    plt.savefig(plotprefix+'_k'+str(kmax)+'_nomoto06_K05_SiH-FeH.png',bbox_inches='tight')

    if options.showplots:
        plt.show()
