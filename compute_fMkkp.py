import numpy as np
import pylab as plt
import time
import h5py
from optparse import OptionParser
import functools
from multiprocessing import Pool
import sys

import karlsson
import util

def run_compute_ckkp(envname,sfrname,maxkmax,
                     tmin=0,tmax=1000,
                     cutoff=10**-7,numprocs=1):
    uII,uIII,muII,muIII,wII,wIII = util.load_sfr(envname,sfrname)
    print "Starting calc_ckkp"
    start = time.time()
    cgrid = karlsson.calc_ckkp(maxkmax,wII,muII,wIII,muIII,uII,tmin=tmin,tmax=tmax,
                               numprocs=numprocs)
    print "calc_ckkp time:",time.time()-start
    maxkmax = cgrid.shape[0]

    kmax = -1
    for irow in range(maxkmax):
        if np.nansum(cgrid[irow,:]) > cutoff: kmax=irow+1
    print "kmax:",kmax

    filename_ckkp = karlsson.get_ckkp_filename(envname,sfrname)
    f = h5py.File(filename_ckkp,'w')
    f.attrs['kmax']=kmax
    f['cgrid']=cgrid
    f.close()

def run_compute_fMkkp(envname,sfrname,
                      logMmin=0,logMmax=8,logdM=.01,
                      numprocs=1):
    Vmix = util.load_Vmix(envname)
    kmax,cgrid = util.load_ckkp(envname,sfrname)
    uII,uIII,muII,muIII,wII,wIII = util.load_sfr(envname,sfrname)

    Mbins,Mplot = karlsson.logMbinsMplot(logMmin,logMmax,logdM)
    np.save(karlsson.get_Mbinskkp_filename(envname,sfrname),Mbins)
    np.save(karlsson.get_Mplotkkp_filename(envname,sfrname),Mplot)
    
    tarr,tauarr,Mmixgrid = util.load_mmixgrid(envname)
    print "Starting get_DMlist with",numprocs,"processors"
    start = time.time()
    DMlist = karlsson.get_DMlist(Mbins,tarr,tauarr,Mmixgrid,numprocs=numprocs)
    print "get_DMlist time:",time.time()-start
    sys.stdout.flush()
    
    filename = 'MMIXDISTR/'+envname+'_'+sfrname+'_fMkkp.npy'
    start = time.time()
    fMkkp = karlsson.calc_fMkkp_onearr(kmax,DMlist,Vmix,
                                       wII,muII,uII,wIII,muIII,uIII,
                                       numprocs=numprocs)
    print "calc_fMkkp_onearr:",time.time()-start
    np.save(filename,fMkkp)

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('-j','--numprocs',action='store',type='int',dest='numprocs',default=1)
    parser.add_option('--maxkmax',action='store',type='int',dest='maxkmax',default=150)
    options,args = parser.parse_args()
    
    envname = args[0]; sfrname = args[1]
    print "envname:",envname
    print "sfrname:",sfrname
    run_compute_ckkp(envname,sfrname,options.maxkmax,numprocs=options.numprocs)
    #run_compute_fMkkp(envname,sfrname,numprocs=options.numprocs)
