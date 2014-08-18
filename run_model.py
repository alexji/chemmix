import numpy as np
import util
import backgroundmodel as bgm
import sfrmodel as sfu
from optparse import OptionParser
import time
import h5py

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('-j','--numprocs',action='store',type='int',dest='numprocs',default=1)
    parser.add_option('-m','--mmixgrid',action='store_true',dest='mmixgrid',default=False)
    parser.add_option('-s','--sfr',action='store_true',dest='sfr',default=False)
    parser.add_option('-c','--ckk',action='store_true',dest='ckk',default=False)
    parser.add_option('-f','--fMkk',action='store_true',dest='fMkk',default=False)
    parser.add_option('-y','--chemhist',action='store_true',dest='chemhist',default=False)
    parser.add_option('--k2max',action='store',type='int',dest='k2max',default=15)
    parser.add_option('--k3max',action='store',type='int',dest='k3max',default=15)
    parser.add_option('--hires',action='store_true',dest='hires',default=False)
    options,args = parser.parse_args()
    envname,sfrname=args
    print envname,sfrname
    if options.hires and ('_hires' not in envname): envname=envname+'_hires'

    k2max = options.k2max; k3max = options.k3max

    if options.mmixgrid:
        print "Computing mmixgrid:",envname
        mmixgrid = bgm.compute_mmix_grid(envname,numprocs=options.numprocs,hires=options.hires)
        np.save(util.fname_mmixgrid(envname),mmixgrid)
        print "Saved mmixgrid to "+util.fname_mmixgrid(envname)
    if options.sfr:
        print "Computing sfr:",envname,sfrname
        start = time.time()
        tarr = bgm.gettarr(envname)
        sfu.compute_umufns(envname,sfrname,tarr,numprocs=options.numprocs)
        print "Done! Time: %f" % (time.time()-start)
    if options.ckk:
        print "Computing ckk",envname,sfrname
        print "k2max: %i, k3max: %i" % (k2max,k3max)
        start = time.time()
        cgrid = bgm.compute_ckk(envname,sfrname,k2max,k3max,numprocs=options.numprocs,hires=options.hires)
        f = h5py.File(util.fname_ckk(envname,sfrname),'w')
        f.attrs['k2max']=k2max; f.attrs['k3max']=k3max
        f['cgrid'] = cgrid
        f.close()
        print "Done! Time: %f" % (time.time()-start)
        print "Saved cgrid to "+util.fname_ckk(envname,sfrname)
        
    if options.fMkk:
        print "Computing fMkk"
        print "k2max: %i, k3max: %i" % (k2max,k3max)
        _k2max,_k3max,cgrid = util.load_ckk(envname,sfrname)
        assert _k2max==k2max and _k3max==k3max
        bgm.compute_fMkk(envname,sfrname,k2max,k3max,numprocs=options.numprocs,hires=options.hires)
        
    if options.chemhist:
        print "chemhist not implemented yet"
