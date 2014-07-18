import numpy as np
import h5py
import karlsson
import pickle
import solveforu as sfu
from tophat import TopHat


import backgroundmodel


def load_Vmix(envname,get_th=False,get_thrho=False,getVmixdot=False,getarrays=False):
    Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,Dt,uSN,Vmax = backgroundmodel.envparams(envname)
    if getarrays: myfn = backgroundmodel.get_Vmixarr
    else: myfn = backgroundmodel.get_Vmixfn

    tarr = backgroundmodel.gettarr(envname,hires=True)
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    rhofn = th.get_rho_of_t_fn()
    if getVmixdot:
        Vmix,Vmixdot = myfn(tarr,rhofn,Dt,logMdil,Vmax,getVmixdot=True)
        if get_thrho: return th,rhofn,Vmix,Vmixdot
        if get_th: return th,Vmix,Vmixdot
        return Vmix,Vmixdot
    else:
        Vmix = myfn(tarr,rhofn,Dt,logMdil,Vmax)
        if get_thrho: return th,rhofn,Vmix
        if get_th: return th,Vmix
        return Vmix


def default_filename(envname,sfrname,postfix=''):
    return envname+'_'+sfrname+postfix

def count_kkp(kmax,kpmax):
    return kmax + kpmax*(kpmax+1)/2 + (kmax-kpmax)*kpmax

def load_mmixgrid(envname):
    f = h5py.File(karlsson.get_mmixgrid_filename(envname),'r')
    tarr = np.array(f['tarr']); tauarr = tarr
    Mmixgrid = np.array(f['Mmixgrid'])
    f.close()
    return tarr,tauarr,Mmixgrid

def get_sfrparams(sfrname,envname=None,verbose=False):
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
    elif 'flat' in name:
        assert envname != None
        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil = karlsson.envparams(envname)
        Dt,uSN = karlsson.get_Dt_uSN(vturb,lturb,nSN,trecovery)
        ttII = ttIII + trecovery
        u0III = uSN
        if 'flatfix' in name:
            pass
        elif 'flat10x' in name:
            u0III = uSN; u0II = 10*u0III
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

def load_ckkp(envname,sfrname,full_grid=False):
    filename_ckkp = karlsson.get_ckkp_filename(envname,sfrname)
    f = h5py.File(filename_ckkp,'r')
    kmax = f.attrs['kmax']
    kpmax = f.attrs['kpmax']
    cgrid = np.array(f['cgrid'])
    f.close()
    if full_grid:
        return kmax,kpmax,cgrid
    else:
        return kmax,kpmax,cgrid[0:kmax,0:(kpmax+1)]

def load_chemgridhist(envname,sfrname,postfix):
    f = open(karlsson.get_chemgridhist_filename(envname,sfrname,postfix=postfix),'r')
    mydict = pickle.load(f)
    f.close()
    return mydict

def hist2d2pdf(H,xedges,yedges):
    nbins_x = len(xedges)-1; nbins_y = len(yedges)-1
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    pdf = (H*(x_bin_sizes*y_bin_sizes).T)
    return pdf

