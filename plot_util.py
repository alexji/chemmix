import numpy as np
import pylab as plt
from matplotlib.ticker import NullFormatter
import scipy.optimize as so
from karlsson import warn_bins

def plot_chemarr_hists(chemgrid,elemarr,binslist):
    Nstar,Nelem = chemgrid.shape
    assert Nelem == len(elemarr) and Nelem == len(binslist)
    f,axarr = plt.subplots(Nelem,Nelem,sharex=True,sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    minarr = []; maxarr = []
    for bins in binslist:
        minarr.append(bins[0]); maxarr.append(bins[-1])

    for irow in range(Nelem):
        for icol in range(irow+1):
            ax = axarr[irow,icol]
            xdata = chemgrid[:,icol]
            xbins = binslist[icol]
            ydata = chemgrid[:,irow]
            ybins = binslist[irow]
            
            density_contour(xdata, ydata, xbins, ybins, ax=ax)
            ax.set_xlim(minarr[icol],maxarr[icol])
            ax.set_ylim(minarr[irow],maxarr[irow])
            
            if irow==Nelem-1: ax.set_xlabel(elemarr[icol])
            if icol==0:       ax.set_ylabel(elemarr[irow])
        for icol in range(irow+1,Nelem):
            axarr[irow,icol].axis('off')
    
    plt.show()
    return f

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


def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level
 
def density_contour(xdata, ydata, xbins, ybins, ax=None, my2dhist=None, **contour_kwargs):
    """ Create a density contour plot. 
    Modified from https://gist.github.com/adrn/3993992
    
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    xbins : int/array
        Number of bins along x dimension/bin edges
    ybins : int/array
        Number of bins along y dimension/bin edges
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """
    
    if my2dhist != None:
        H,xedges,yedges = my2dhist
        nbins_x = len(xedges)-1; nbins_y = len(yedges)-1
    else:
        H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(xbins,ybins), normed=True)
        try: #xbins, ybins are arrays
            nbins_x = len(xbins)-1; nbins_y = len(ybins)-1
        except: #xbins, ybins are integers
            nbins_x = xbins; nbins_y = ybins
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    
    pdf = (H*(x_bin_sizes*y_bin_sizes).T)
    # if H was unnormalized, need to redo pdf
    if np.abs(np.sum(pdf)-1) > 10**-6: pdf = H/float(np.sum(H))
    
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [one_sigma, two_sigma, three_sigma]
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
    
    return contour
 
def plot1dhist(h,x,ax=None,**kwargs):
    if ax==None: ax = plt.gca()
    xplot = (x[:-1]+x[1:])/2.
    ax.plot(xplot,h,drawstyle='steps-mid',**kwargs)

def test_density_contour():
    norm = np.random.normal(10., 15., size=(12540035, 2))
    f,(ax1,ax2) = plt.subplots(2)
    density_contour(norm[:,0], norm[:,1], 100, 75, ax=ax1)
    density_contour(norm[:,0], norm[:,1], np.arange(-80,81), np.arange(-80,81), ax=ax2)
    plt.show()

