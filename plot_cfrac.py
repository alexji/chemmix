import numpy as np
import pylab as plt
from optparse import OptionParser

import karlsson
import util

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--save',action='store_true',dest='save',default=False)
    parser.add_option('-c','--cfecrit',action='store',type="float",dest='CFecrit',default=0.75)
    options,args = parser.parse_args()
    envname,sfrname,postfix = args
    
    histdict = util.load_chemgridhist(envname,sfrname,postfix)
    yII = histdict['yII']
    elemnames = np.array(yII.elemnames)
    numyields = len(elemnames)
    CFecrit = options.CFecrit
    iC  = np.where('C'==elemnames)[0][0]  #0
    iFe = np.where('Fe'==elemnames)[0][0] #5
    assert iC < iFe 

    histdict = util.load_chemgridhist(envname,sfrname,postfix)
    H,edgesC,edgesFe = histdict[(iFe,iC)]
    pdf = util.hist2d2pdf(H,edgesC,edgesFe)

    H = H.T #put Fe on the X axis
    midFe = (edgesFe[1:]+edgesFe[:-1])/2.
    midC  = (edgesC[1:]+edgesC[:-1])/2.

    cfrac = np.zeros(len(midFe))
    for i,Fe in enumerate(midFe):
        iiCrich = midC-Fe > CFecrit
        cfrac[i] = np.sum(H[i,iiCrich])/np.sum(H[i,:])
    
    fig = plt.figure()
    ax = plt.gca()
    #ax.plot(midFe,cfrac,'k.-',drawstyle='steps-mid')
    ax.plot(edgesFe[:-1],cfrac,'k-',drawstyle='steps')
    ax.set_ylim((0,1.1))
    ax.set_xlabel('[Fe/H]'); ax.set_ylabel('cfrac')
    if options.save: 
        plt.savefig('PLOTS/cfrac_'+util.default_filename(envname,sfrname,postfix)+'.png')
    else: plt.show()
        
