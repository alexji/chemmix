import numpy as np
import numpy.random as random
import util
import backgroundmodel as bgm
import yields

from multiprocessing import Pool
import functools
import time
import sys
import pickle

def draw_from_distr(N,x,pdf,seed=None,eps=10**-10):
    if N==0: return np.array([])
    if seed!=None:
        random.seed(seed)
    unifarr = random.rand(N)
    if len(pdf.shape)==1:
        cdf = np.cumsum(pdf)
    else:
        cdf = np.cumsum(pdf.reshape((-1,)))
    assert np.abs(cdf[-1] - 1.0) < eps, 'cdf[-1] != 1.0 (to tolerance of %e): %f' % (eps,cdf[-1])
    cdf[-1]=1
    output = np.array(x[np.searchsorted(cdf,unifarr)])
    if len(pdf.shape)==1:
        return output
    else:
        return output.reshape(pdf.shape)

def _hist_chemgridkk(yieldnumratio,binlist,elemnames,verbose=False):
    Nstars,numyields = yieldnumratio.shape
    assert numyields == len(elemnames) and numyields == len(binlist)
    # do NOT normalize histograms because will be weighted later
    outdict = {}
    for irow,erow in enumerate(elemnames):
        for icol,ecol in enumerate(elemnames):
            if icol > irow: continue
            if irow==icol: 
                outdict[(irow,icol)] = np.histogram(yieldnumratio[:,icol],
                                                    bins=binlist[icol])
            else:
                myhist2d = np.histogram2d(yieldnumratio[:,icol],
                                          yieldnumratio[:,irow],
                                          bins=[binlist[icol],binlist[irow]])
                outdict[(irow,icol)] = myhist2d
    return outdict

def _montecarloyields(tasknum,k2max,k3max,Mplot,
                      fMkkII,fMkkIII,yII,yIII,
                      masstonum,Nstars,binlist,
                      timethis=True,tmplabel=''):
    savetmpfiles = (tmplabel != '') #TODO
    if timethis: start = time.time()
    numyields = yII.numyields

    kII,kIII = divmod(tasknum,k3max+1)
    if kII==0 and kIII==0: return {}
    this_fMkkII  = fMkkII[kII,kIII,:]
    this_fMkkIII = fMkkIII[kII,kIII,:]
    yieldratioII = None; yieldratioIII = None
    if kII > 0:
        mixmassII = draw_from_distr(kII*Nstars,Mplot,this_fMkkII)
        yieldarrII = yII.draw_yields(kII*Nstars)
        yieldratioII = np.reshape((yieldarrII.T/mixmassII).T,(Nstars,kII,numyields))
    if kIII > 0:
        mixmassIII = draw_from_distr(kIII*Nstars,Mplot,this_fMkkIII)
        yieldarrIII = yIII.draw_yields(kIII*Nstars)
        yieldratioIII = np.reshape((yieldarrIII.T/mixmassIII).T,(Nstars,kIII,numyields))

    if yieldratioII != None and yieldratioIII != None:
        yieldratio = np.concatenate((yieldratioII,yieldratioIII),axis=1)
    elif yieldratioII == None:
        yieldratio = yieldratioIII
    elif yieldratioIII == None:
        yieldratio = yieldratioII
    yieldratio = np.sum(yieldratio,axis=1) #mass ratios
    yieldratio  = yieldratio*masstonum     #number ratios
    if timethis: 
        print "task %i (k2=%i, k3=%i): %f" % (tasknum,kII,kIII,time.time()-start); sys.stdout.flush()

    # compute/return histograms
    yieldratio = yields.convert_to_solar(yII.elemnames,yieldratio)
    outdict = _hist_chemgridkk(yieldratio,binlist,yII.elemnames)
    outdict['yII'] = yII
    outdict['yIII']= yIII
    outdict['binlist']=binlist
    return outdict

def montecarloyields(envname,sfrname,postfix,
                     yII,yIII,binlist,
                     XH=0.75,
                     Nstars=10**5,numprocs=1):
    Mplot  = np.load(util.fname_Mplotkk(envname,sfrname))
    fMkkII = np.load(util.fname_fMkkII(envname,sfrname))
    fMkkIII= np.load(util.fname_fMkkIII(envname,sfrname))
    
    assert fMkkII.shape == fMkkIII.shape
    k2max,k3max,numbins = fMkkII.shape; k2max -= 1; k3max -= 1
    assert numbins == len(Mplot)

    masstonum = 1.0/(yII.elemmass * XH)

    start = time.time()
    fname = util.fname_chemgridhistlist(envname,sfrname,postfix)
    pool = Pool(numprocs)

    myfunc = functools.partial(_montecarloyields,k2max=k2max,k3max=k3max,
                               Mplot=Mplot,fMkkII=fMkkII,fMkkIII=fMkkIII,
                               yII=yII,yIII=yIII,masstonum=masstonum,Nstars=Nstars,
                               binlist=binlist)
    dictlist = pool.map(myfunc,range((k2max+1)*(k3max+1)))
    pool.close()
    print "Finished! Time %f" % (time.time()-start)
    f = open(fname,'w')
    pickle.dump(dictlist,f)
    f.close()

if __name__=="__main__":
    envname='atomiccoolinghalo'
    sfrname='fixTS500'

    binwidth = .1
    binlist = [np.arange(-10,0+binwidth,binwidth) for i in range(6)]
    salimf  = util.relative_imf(np.array([13,15,18,20,25,30,40.]),2.35,norm=True) #salpeter imf
    popIIy  = yields.ryN06(salimf)
    popIIIy_arr = [yields.ryI05N06(p) for p in [0,.1,.3,.5,.7,.9,1.]]

    for popIIIy in popIIIy_arr:
        postfix = 'p%0.2f' % (popIIIy.p)
        montecarloyields(envname,sfrname,postfix,
                         popIIy,popIIIy,binlist,
                         numprocs=8)
