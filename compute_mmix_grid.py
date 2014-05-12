import numpy as np
import time

import karlsson
from tophat import TopHat

import h5py

if __name__=="__main__":
    #tmin = 0.0; tmax = 1000.0; dt = 0.1
    #VMIX = karlsson.Vmix_K05
    #RHO = karlsson.rho_K05A

    #FILENAME = "TESTa_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    #print "starting compute_mmix_grid (may take a long time)"
    #start = time.time()
    #tarr,tauarr,Mmixgrid = karlsson.calc_Mmix_grid(VMIXDOT,RHO,
    #    tmin=tmin,tmax=tmax,dt=dt,numprocs=1)
    #np.savez(FILENAME,tarr=tarr,tauarr=tauarr,Mmixgrid=Mmixgrid)
    #print "compute_mmix_grid time:",time.time()-start
    
    #FILENAME = "TESTb_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    #print "starting compute_mmix_grid_fast (may take a long time)"
    #start = time.time()
    #tarr,tauarr,Mmixgrid = karlsson.calc_Mmix_grid_fast(VMIX,RHO,
    #    tmin=tmin,tmax=tmax,dt=dt)
    #np.savez(FILENAME,tarr=tarr,tauarr=tauarr,Mmixgrid=Mmixgrid)
    #print "compute_mmix_grid_fast time:",time.time()-start

    #FILENAME = "TESTc_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    #print "starting compute_mmix_grid_fast_trap (may take a long time)"
    #start = time.time()
    #tarr,tauarr,Mmixgrid = karlsson.calc_Mmix_grid_fast_trap(VMIX,RHO,
    #    tmin=tmin,tmax=tmax,dt=dt)
    #np.savez(FILENAME,tarr=tarr,tauarr=tauarr,Mmixgrid=Mmixgrid)
    #print "compute_mmix_grid_fast_trap time:",time.time()-start

    #tmin = 0.0; tmax = 1000.0; dt = 0.1
    #th = TopHat()
    #RHO = th.get_rho_of_t_fn()
    #VMIX = karlsson.get_Vmixfn_K08(RHO)
    #FILENAME = "Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    #print "computing K08 with midpoint (may take a long time)"
    #start = time.time()
    #tarr,tauarr,Mmixgrid = karlsson.calc_Mmix_grid_fast(VMIX,RHO,
    #    tmin=tmin,tmax=tmax,dt=dt)
    #np.savez(FILENAME,tarr=tarr,tauarr=tauarr,Mmixgrid=Mmixgrid)
    #print "compute_mmix_grid_fast time:",time.time()-start

    tmin = 0.0; tmax = 100.0; dt = 0.03
    th = TopHat()
    RHO = th.get_rho_of_t_fn()
    VMIX = karlsson.get_Vmixfn_K08(RHO)
    FILENAME = "Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    f = h5py.File(FILENAME+'.hdf5',"w")
    print "computing K08 hires with midpoint (may take a long time)"
    start = time.time()
    tarr,tauarr,Mmixgrid = karlsson.calc_Mmix_grid_fast(VMIX,RHO,
        tmin=tmin,tmax=tmax,dt=dt)
    print "compute_mmix_grid_fast time:",time.time()-start
    f.attrs['tmin'] = tmin
    f.attrs['tmax'] = tmax
    f.attrs['dt'] = dt
    f.attrs['rho'] = 'tophat'
    f.attrs['vmix'] = 'K08'
    f['tarr']=tarr
    f['tauarr']=tauarr
    f['Mmixgrid']=Mmixgrid
    f.close()
    #np.savez(FILENAME,tarr=tarr,tauarr=tauarr,Mmixgrid=Mmixgrid)
