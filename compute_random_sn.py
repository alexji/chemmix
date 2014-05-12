import numpy as np
import karlsson

import time
from multiprocessing import Pool
import functools
import numpy.random as random

from yields import interp_nomoto06_yields

def relative_imf(marr,alpha):
    Mmax = float(marr[-1])
    return np.array([(M/Mmax)**(-alpha) for M in marr])

def run_one_star(j,k,Mplot,fMk,sntypearr,imfpdf,yieldfn,masstonum,gaussianprofile=False,carbonenhanceprob=0.0,carbonenhancefactor=100.):
    # draw k mixing masses, SN types
    mixing_masses = karlsson.draw_from_distr(k,Mplot,fMk)
    sn_types = karlsson.draw_from_distr(k,sntypearr,imfpdf)
    # get the k * numyields element production array
    yields = np.array([yieldfn(sn_type) for sn_type in sn_types])
    # with a certain probability, add carbon enhancement
    if carbonenhanceprob > 0:
        yields[random.rand(k) < carbonenhanceprob, 0] *= carbonenhancefactor
    # weighted sum to get 1 * numyields, which goes into output
    if gaussianprofile:
        pass #TODO!!!
    else:
        weightedyields = masstonum*np.sum((yields.transpose()/mixing_masses),1).transpose()
    return weightedyields

def compute_random_sn_yields(FILENAME,kmax,Mbins,XH=0.75,Nstars=10**5,numprocs=1):
    pass

if __name__=="__main__":
    #tmin = 0.0; tmax = 1000.0; dt = 0.1
    tmin = 0.0; tmax = 1000.0; dt = 0.03
    FILENAME="Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    Nstars = 10**5
    numyields = 6
    kmin = 1
    kmax = 10
    XH = 0.75
    numprocs=20
    Mbins = 10**np.arange(2,8.01,.01)
    Mplot = 10**(np.arange(2,8.0,.01)+.005)
    elemmass = np.array([12.0,16.0,24.3,28.1,40.1,55.8]) #mZ/mH
    masstonum = 1.0/(elemmass * XH)

    output = np.zeros((Nstars,numyields,kmax))

    pool = Pool(numprocs)

    imfarr = np.arange(13,41)
    imfpdf = relative_imf(imfarr,2.35)
    imfpdf = imfpdf/np.sum(imfpdf)
    sntypearr = np.arange(1,28+1)
    
    #for k in xrange(kmin,kmax+1):
    #    print "Starting type",k
    #    start = time.time()
    #    fMk = np.load(FILENAME+'_fM'+str(k)+'.npy')
    #    starmaker = functools.partial(run_one_star,k=k,Mplot=Mplot,fMk=fMk,sntypearr=sntypearr,imfpdf=imfpdf,yieldfn=interp_nomoto06_yields,masstonum=masstonum)
    #    this_output = pool.map(starmaker,xrange(Nstars))
    #    output[:,:,k-1] = np.array(this_output)
    #    print "  time:",time.time()-start
    #np.save(FILENAME+"_chemgrid_N06_N"+str(Nstars),output)
    
    for k in xrange(kmin,kmax+1):
        print "Starting type",k
        start = time.time()
        fMk = np.load(FILENAME+'_fM'+str(k)+'.npy')
        starmaker = functools.partial(run_one_star,k=k,Mplot=Mplot,fMk=fMk,sntypearr=sntypearr,imfpdf=imfpdf,yieldfn=interp_nomoto06_yields,masstonum=masstonum,carbonenhanceprob=0.001)
        this_output = pool.map(starmaker,xrange(Nstars))
        output[:,:,k-1] = np.array(this_output)
        print "  time:",time.time()-start
    np.save(FILENAME+"_chemgrid_N06_Cenhanced_p0.001_N"+str(Nstars),output)
