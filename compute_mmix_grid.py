import numpy as np
import time

import karlsson
from tophat import TopHat

import h5py
from optparse import OptionParser

def run_compute_mmix_grid(filename,Mhalo,zvir,vturb,lturb,tmax,dt,tmin=0.0,logMdil=5):
    """
    Assume a TopHat density function
    Mhalo: Msun
    vturb: km/s
    lturb: kpc
    tmax, dt: Myr
    """
    vturb *= 3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturb * lturb / 3.0
    print filename
    print "Mhalo",Mhalo,"zvir",zvir
    print "vturb",vturb,"lturb",lturb,"Dt",Dt
    print "tmax",tmax,"dt",dt

    th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
    RHO = th.get_rho_of_t_fn()
    VMIX = karlsson.get_Vmixfn_K08(RHO,Dt=Dt,Mdil=10**logMdil)
    MMIXGRID_FILENAME = karlsson.get_mmixgrid_filename(filename)

    f = h5py.File(MMIXGRID_FILENAME,"w")
    print "computing mmixgrid with midpoint (may take a long time)"
    start = time.time()
    tarr,tauarr,Mmixgrid = karlsson.calc_Mmix_grid_fast(VMIX,RHO,
        tmin=tmin,tmax=tmax,dt=dt)
    print "calc_Mmix_grid_fast time:",time.time()-start
    f.attrs['tmin'] = tmin
    f.attrs['tmax'] = tmax
    f.attrs['dt'] = dt
    f.attrs['rho'] = 'tophat'
    f.attrs['vmix'] = 'K08_Dt'+str(Dt)
    f.attrs['Mhalo'] = Mhalo

    f['tarr']=tarr
    f['tauarr']=tauarr
    f['Mmixgrid']=Mmixgrid
    f.close()
    print "Wrote",MMIXGRID_FILENAME

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("--minihalo",action='store_true',dest='minihalo',default=False)
    parser.add_option("--atomiccoolinghalo",action='store_true',dest='atomiccoolinghalo',default=False)
    parser.add_option("--mdil",action='store',type='int',dest='logMdil',default=5)
    options,args = parser.parse_args()
    logMdil = options.logMdil
    
    ## Compute Minihalo
    if options.minihalo:
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_minihalo()
        fact = 1.25
        tmax = 100. * fact #Myr
        dt = 0.003*fact
        if dt == 0.01*fact:
            filename='lores_minihalo'
        if dt == 0.003*fact:
            filename='minihalo'
        if dt == 0.001*fact:
            filename='hires_minihalo'
        if logMdil != 5: filename += str(logMdil)
        run_compute_mmix_grid(filename,Mhalo,zvir,vturb,lturb,tmax,dt,logMdil=float(logMdil))

    ## Compute Fat Minihalo
    #filename='fatminihalo'

    ## Compute Atomic Cooling Halo
    if options.atomiccoolinghalo:
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalo()
        tmax = 1000. #Myr
        dt = 0.03
        if dt == 0.1:
            filename='lores_atomiccoolinghalo'
        if dt == 0.03:
            filename='atomiccoolinghalo'
        if dt == 0.01:
            filename='hires_atomiccoolinghalo'
        if logMdil != 5: filename += str(logMdil)
        run_compute_mmix_grid(filename,Mhalo,zvir,vturb,lturb,tmax,dt,logMdil=float(logMdil))
