import numpy as np
import h5py
import karlsson
import solveforu as sfu
from tophat import TopHat

def load_mmixgrid(envname):
    f = h5py.File(karlsson.get_mmixgrid_filename(envname),'r')
    tarr = np.array(f['tarr']); tauarr = tarr
    Mmixgrid = np.array(f['Mmixgrid'])
    f.close()
    return tarr,tauarr,Mmixgrid

def load_Vmix(envname,get_th=False,get_thrho=False):
    Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil = karlsson.envparams(envname)
    Dt,uSN = karlsson.get_Dt_uSN(vturb,lturb,nSN,trecovery)
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    RHO = th.get_rho_of_t_fn()
    Vmix = karlsson.get_Vmixfn_K08(RHO,Dt=Dt,Mdil=10**logMdil)
    if get_thrho: return th,RHO,Vmix
    if get_th: return th,Vmix
    return Vmix

def get_sfrparams(sfrname,verbose=False):
    try:
        name,ts = sfrname.split('TS')
    except ValueError:
        raise ValueError("Invalid sfrname: "+sfrname+" (requires 'TS' to denote tstop)")
    ts = int(ts)
    u0II = .4;   u0III = .04
    ttII = 110; ttIII = 100

    if 'fid' in name:
        pass
    elif 'utriple' in name:
        u0III *= 3.; u0II *= 3.
    elif 'uIIItriple' in name:
        u0III *= 3.
    elif 'uthird' in name:
        u0III /= 3.; u0II /= 3.
    elif 'uIIIthird' in name:
        u0III /= 3.
    elif 'uten' in name:
        u0III *= 10.; u0II *= 10.
    elif 'uIIIten' in name:
        u0III *= 10.
    elif 'utenth' in name:
        u0III /= 10.; u0II /= 10.
    elif 'uIIItenth' in name:
        u0III /= 10.
    else:
        raise ValueError("Invalid sfrname: "+sfrname)

    if verbose:
        print name
        print "u0II %3.2e u0III %3.2e" % (u0II,u0III)
        print "ttII %i ttIII %i" % (tII,ttIII)
        print "ts %i" % (ts)
    return u0II,u0III,ttII,ttIII,ts

def get_sfrlabel(envname,sfrname):
    return envname+'__'+sfrname

def load_sfr(envname,sfrname):
    label = get_sfrlabel(envname,sfrname)
    uII   = sfu.loaduIIfn(label)
    uIII  = sfu.loaduIIIfn(label)
    muII  = sfu.loadmuIIfn(label)
    muIII = sfu.loadmuIIIfn(label)
    wII   = karlsson.wISM_K05
    wIII  = karlsson.wISM_III
    return uII,uIII,muII,muIII,wII,wIII

def load_ckkp(envname,sfrname):
    filename_ckkp = karlsson.get_ckkp_filename(envname,sfrname)
    f = h5py.File(filename_ckkp,'r')
    kmax = f.attrs['kmax']
    cgrid = f['cgrid']
    f.close()
    return kmax,cgrid

def load_fMkkp(envname,sfrname,k,kp):
    fMkkp_filename = karlsson.get_fMkkp_filename(envname,sfrname,k,kp)
    return np.load(fMkkp_filename)
