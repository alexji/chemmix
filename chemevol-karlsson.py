import numpy as np
from karlsson import *

from yields import interp_nomoto06_yields

import time
from multiprocessing import Pool
import functools

def relative_imf(marr,alpha):
    Mmax = float(marr[-1])
    return np.array([(M/Mmax)**(-alpha) for M in marr])

def run_one_star(j,k,Mplot,fMk,sntypearr,imfpdf,yieldfn):
    # draw k mixing masses
    mixing_masses = draw_from_distr(k,Mplot,fMk)
    # draw k SN types
    sn_types = draw_from_distr(k,sntypearr,imfpdf)
    # get the k * numyields element production array
    yields = np.array([yieldfn(sn_type) for sn_type in sn_types])
    # weighted sum to get 1 * numyields, which goes into output
    weightedyields = np.sum((yields.transpose()/mixing_masses),1).transpose()
    return weightedyields
    #output[j,:,k-1] = weightedyields    

if __name__=="__main__":
    NUMPROCS = 8

    kmax = 10
    Nstars = 10**6
    XH= 0.75
    numyields = 6
    
    output = np.zeros((Nstars,numyields,kmax))

    Mbins = np.linspace(0,4,201)*10**6
    Mplot = (Mbins[:-1]+Mbins[1:])/2.0

    imfarr = np.arange(13,41)
    imfpdf = relative_imf(imfarr,2.35)
    imfpdf = imfpdf/np.sum(imfpdf)
    sntypearr = np.arange(1,28+1)

    pool = Pool(NUMPROCS)

    #totaloutput = []

    for k in xrange(1,kmax+1):
        # read fMk
        fMk = np.load('fM'+str(k)+'_K05_10000_dt0.1.npy')
        starmaker = functools.partial(run_one_star,k=k,Mplot=Mplot,fMk=fMk,sntypearr=sntypearr,imfpdf=imfpdf,yieldfn=interp_nomoto06_yields)
        this_output = pool.map(starmaker,xrange(Nstars))
        output[:,:,k-1] = np.array(this_output)
        #print np.array(output).shape
        #totaloutput.append(output)
        #output[j,:,k-1] = weightedyields
    np.save('k10_nomoto06_K05_N'+str(Nstars)+'_res10000_dt0.1',output)
