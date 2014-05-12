import numpy as np
import time

import karlsson
from tophat import TopHat

import h5py

def run_compute_mmix_grid(filename,Mhalo,vturb,lturb,tmax,dt,tmin=0.0):
    """
    Assume a TopHat density function
    Mhalo: Msun
    vturb: km/s
    lturb: kpc
    tmax, dt: Myr
    """
    vturb *= 3.16/3.08 * .001 #km/s to kpc/yr
    Dt =  vturb * lturb / 3.0
    print filename
    print "Mhalo",Mhalo,"vturb",vturb,"lturb",lturb,"Dt",Dt
    print "tmax",tmax,"dt",dt

    th = TopHat(Mhalo=Mhalo,nvir=0.1,fb=0.1551,mu=1.4)
    RHO = th.get_rho_of_t_fn()
    VMIX = karlsson.get_Vmixfn_K08(RHO,Dt=Dt)
    FILENAME = filename#+"Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)

    f = h5py.File(FILENAME+'.hdf5',"w")
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


if __name__=="__main__":
    ## Compute Minihalo
    filename='mmixgrid_minihalo'
    Mhalo = 10.**6
    vturb = 2. #km/s
    lturb = .01 #kpc
    tmax = 100. #Myr
    dt = 0.01
    run_compute_mmix_grid(filename,Mhalo,vturb,lturb,tmax,dt)

    ## Compute Fat Minihalo
    filename='mmixgrid_fatminihalo'

    ## Compute Atomic Cooling Halo
    filename='mmixgrid_atomiccoolinghalo'
    Mhalo = 10.**8
    vturb = 10. #km/s
    lturb = .1 #kpc
    tmax = 1000. #Myr
    dt = 0.1
    run_compute_mmix_grid(filename,Mhalo,vturb,lturb,tmax,dt)

    
