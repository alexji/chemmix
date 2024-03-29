import numpy as np
import pylab as plt
from matplotlib.ticker import NullFormatter
import h5py
from optparse import OptionParser

import karlsson

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
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_N06i.hdf5'
#    plotprefix = 'atomiccoolinghalo_N06i'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo_chemgrid_N06ip0.100f100.0.hdf5'
#    plotprefix = 'minihalo_N06ip0.100f100.0'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_N06ip0.100f100.0.hdf5'
#    plotprefix = 'atomiccoolinghalo_N06ip0.100f100.0'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo_chemgrid_HW10E1.2S4m0.hdf5'
#    plotprefix = 'minihalo_HW10E1.2S4m0'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_HW10E1.2S4m0.hdf5'
#    plotprefix = 'atomiccoolinghalo_HW10E1.2S4m0'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo_chemgrid_mixN06HW10p0.1.hdf5'
#    plotprefix = 'minihalo_NHWE1.2S4m0'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_mixN06HW10p0.1.hdf5'
#    plotprefix = 'atomiccoolinghalo_NHWE1.2S4m0'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo_chemgrid_mixN06HW10p0.5.hdf5'
#    plotprefix = 'minihalo_NHWE1.2S4m0'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_mixN06HW10p0.5.hdf5'
#    plotprefix = 'atomiccoolinghalo_NHWE1.2S4m0'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo_chemgrid_mixN06HW10p0.5.hdf5'
#    plotprefix = 'minihalo_NHWE1.2S4m0'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_mixN06HW10p0.5.hdf5'
#    plotprefix = 'atomiccoolinghalo_NHWE1.2S4m0'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo_chemgrid_I05T07_p0.5.hdf5'
#    plotprefix = 'minihalo_I05T07_p0.5'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_I05T07_p0.5.hdf5'
#    plotprefix = 'atomiccoolinghalo_I05T07_p0.5'
#    paramfn = karlsson.params_atomiccoolinghalo
#    filename = 'CHEMGRIDS/minihalo4_chemgrid_I05T07_p0.5.hdf5'
#    plotprefix = 'minihalo4_I05T07_p0.5'
#    paramfn = karlsson.params_minihalo
#    filename = 'CHEMGRIDS/atomiccoolinghalo4_chemgrid_I05T07_p0.5.hdf5'
#    plotprefix = 'atomiccoolinghalo4_I05T07_p0.5'
#    paramfn = karlsson.params_atomiccoolinghalo

    filename = 'CHEMGRIDS/atomiccoolinghalo_lowmass_chemgrid_I05N06_p0.5Mmax.hdf5'
    plotprefix = 'atomiccoolinghalo_lowmass_I05N06_p0.5'
    paramfn = karlsson.params_atomiccoolinghalo_lowmass

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

    print filename
    f = h5py.File(filename,'r')
    chemarr = np.array(f['chemgrid'])
    f.close()

    kmax = options.kmax
    print "kmax to plot:",kmax
    elemnames = ['C', 'O', 'Mg', 'Si', 'Ca', 'Fe']
    chemarr,ck = karlsson.weight_chemgrid(kmax,chemarr,paramfn,
                                          elemnames=elemnames,verbose=True)
    #Nstar,numyields = chemarr.shape
    #asplund09solar = np.array([8.43,8.69,7.60,7.51,6.34,7.50]) - 12 #log10(Z/H) solar
    #for z in range(numyields):
    #    chemarr[:,z] -= asplund09solar[z]
    #    print elemnames[z],np.min(chemarr[:,z]),np.max(chemarr[:,z])
    
    ## BEGIN PLOTS
    plt.figure()
    plt.plot(np.arange(kmax)+1,ck.cumsum(),'ks-')
    plt.xlabel('k'); plt.ylabel(r'cumulative $c_k$')
    plt.xlim(0,kmax+1); plt.ylim((0,1))
    plt.xticks(np.arange(kmax)+1)
    plt.savefig('PLOTS/'+plotprefix+'_k'+str(kmax)+'_ck.png')

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
    plt.savefig('PLOTS/'+plotprefix+'_k'+str(kmax)+'_sixelem.png',bbox_inches='tight')

    plot_2d_hist(chemarr[:,5],chemarr[:,0]-chemarr[:,5],
                 nxbins=100,nybins=100,
                 xlabel='[Fe/H]',ylabel='[C/Fe]')
    plt.savefig('PLOTS/'+plotprefix+'_k'+str(kmax)+'_CFe-FeH.png',bbox_inches='tight')

#    plot_2d_hist(chemarr[:,5],chemarr[:,2]-chemarr[:,5],
#                 nxbins=100,nybins=100,
#                 xlabel='[Fe/H]',ylabel='[Mg/Fe]')
#    plt.savefig('PLOTS/'+plotprefix+'_k'+str(kmax)+'_nomoto06_K05_MgFe-FeH.png',bbox_inches='tight')
#
#    plot_2d_hist(chemarr[:,5],chemarr[:,3],
#                 nxbins=100,nybins=100,
#                 xlabel='[Fe/H]',ylabel='[Si/H]')
#    plt.savefig('PLOTS/'+plotprefix+'_k'+str(kmax)+'_nomoto06_K05_SiH-FeH.png',bbox_inches='tight')

    if options.showplots:
        plt.show()
