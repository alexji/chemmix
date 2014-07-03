import numpy as np
import pylab as plt
import pickle
from optparse import OptionParser

from plot_util import density_contour,plot1dhist
import karlsson
import util

if __name__=="__main__":
    import yields

    parser = OptionParser()
    options,args = parser.parse_args()
    #envname = args[0]; sfrname = args[1]
    envname = 'atomiccoolinghalo'
    sfrname = 'fidTS300'
    postfix = 'testp0.50'

    histdict = util.load_chemgridhist(envname,sfrname,postfix)
    yII = histdict['yII']
    yIII= histdict['yIII']
    #binlist = histdict['binlist']
    elemnames = yII.elemnames
    numyields = len(elemnames)

    fig,axarr = plt.subplots(numyields,numyields,sharex=True)
    for irow,erow in enumerate(elemnames):
        for icol,ecol in enumerate(elemnames):
            ax = axarr[irow,icol]
            key = (irow,icol)
            if icol > irow:
                ax.axis('off')
            elif irow == icol:
                h,x = histdict[key]
                plot1dhist(h,x,ax=ax)
                ax.set_xlim((np.min(x),np.max(x)))
                ax.set_ylim((0,np.max(h)*1.1))
            else:
                my2dhist = histdict[key]
                h,x,y = my2dhist
                density_contour(0,0,0,0,ax=ax,my2dhist=my2dhist)
                ax.set_xlim((np.min(x),np.max(x)))
                ax.set_ylim((np.min(y),np.max(y)))
            if irow==numyields-1: ax.set_xlabel(elemnames[icol])
            if icol==0:           ax.set_ylabel(elemnames[irow])
    plt.show()
