import numpy as np
import pylab as plt
from optparse import OptionParser

from plot_util import density_contour,plot1dhist
import karlsson
import util

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--save',action='store_true',dest='save',default=False)
    options,args = parser.parse_args()
    envname,sfrname,postfix = args
    #envname = args[0]; sfrname = args[1]
    #envname = 'atomiccoolinghalo'
    #sfrname = 'fidTS300'
    #postfix = 'testp0.50'

    #histdict = util.load_chemgridhist(envname,sfrname,postfix)
    histlist = util.load_histlist(envname,sfrname,postfix)
    histdict = histlist[1] #kIII==1
    #yII = histdict['yII']
    #yIII= histdict['yIII']
    #binlist = histdict['binlist']
    #elemnames = yII.elemnames
    elemnames = ['C','O','Mg','Si','Ca','Fe']
    numyields = len(elemnames)
    #TODO
    #import pickle
    #f = open('/spacebase/data/alexji/karlsson-model/tmp.hist','r')
    #histdict = pickle.load(f)
    #f.close()
    #elemnames = ['C','O','Mg','Si','Ca','Fe']
    #numyields = 6

    fig,axarr = plt.subplots(numyields,numyields,sharex=True)
    for irow,erow in enumerate(elemnames):
        for icol,ecol in enumerate(elemnames):
            ax = axarr[irow,icol]
            key = (irow,icol)
            if icol > irow: #reverse
                #ax.axis('off')
                key = (icol,irow)
                my2dhist = histdict[key]
                h,y,x = my2dhist
                density_contour(0,0,0,0,ax=ax,my2dhist=(h.T,x,y))
                ax.set_xlim((np.min(x),np.max(x)))
                ax.set_ylim((np.min(y),np.max(y)))
            elif irow == icol:
                h,x = histdict[key]
                if np.sum(h) > 10: h = h/float(np.sum(h))
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
    if options.save: 
        plt.savefig('PLOTS/chemhists_'+util.default_filename(envname,sfrname,postfix)+'.png')
    else: plt.show()
