import numpy as np
import time
import h5py
from optparse import OptionParser

import util
import karlsson


def run_compute_mmix_grid(envname,tmin,tmax,dt):#filename,Mhalo,zvir,vturb,lturb,tmax,dt,tmin=0.0,logMdil=5):
    Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil = karlsson.envparams(envname)
    th,RHO,VMIX = util.load_Vmix(envname,get_thrho=True)
    MMIXGRID_FILENAME = karlsson.get_mmixgrid_filename(envname)

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
    f.attrs['envname'] = envname
    f.attrs['vmix'] = 'K08_Dt'+str(Dt)+'_logMdil'+str(logMdil)
    f.attrs['Mhalo'] = Mhalo
    f.attrs['zvir'] = zvir

    f['tarr']=tarr
    f['tauarr']=tauarr
    f['Mmixgrid']=Mmixgrid
    f.close()
    print "Wrote",MMIXGRID_FILENAME

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("--hires",action='store_true',dest='hires',default=False)
    parser.add_option("--ultrahires",action='store_true',dest='ultrahires',default=False)
    options,args = parser.parse_args()

    dt = 0.1
    if options.hires: dt = .03
    if options.ultrahires: dt = .001

    envname = args[0]
    tmin = 0
    if 'minihalo' in envname:
        fact = 1.25
        tmax = 100. * fact
        dt = dt/10. * fact
    elif 'atomiccoolinghalo' in envname:
        tmax = 1000.
    else:
        raise ValueError("Invalid envname: "+str(envname))
        
    print "Running",envname
    print "tmin %5.4f tmax %5.1f dt %5.4f" % (tmin,tmax,dt)
    run_compute_mmix_grid(envname,tmin,tmax,dt)
