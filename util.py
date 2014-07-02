import numpy as np
import h5py
import karlsson
import solveforu as sfu
from tophat import TopHat

def default_filename(envname,sfrname):
    return envname+'_'+sfrname

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

def load_ckkp(envname,sfrname,full_grid=False):
    filename_ckkp = karlsson.get_ckkp_filename(envname,sfrname)
    f = h5py.File(filename_ckkp,'r')
    kmax = f.attrs['kmax']
    cgrid = np.array(f['cgrid'])
    f.close()
    if full_grid:
        return kmax,cgrid
    else:
        return kmax,cgrid[0:kmax,0:(kmax+1)]

class interp1d(object):
    """
    Mimics scipy.interpolate interp1d but is pickleable for multiprocessing module
    """
    def __init__(self,x,y,copy=True,bounds_error=True,fill_value=np.nan,assume_sorted=False):
        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value
        x = np.array(x,copy=self.copy)
        y = np.array(y,copy=self.copy)
        if not assume_sorted:
            ind = np.argsort(x)
            x = x[ind]
            y = y[ind] #np.take(y,ind,axis=-1)
        if not issubclass(y.dtype.type,np.inexact):
            y = y.astype(np.float_)
        self.x = x
        self.y = y
        self._y = y
    
    def _evaluate(self,x_new):
        x_new = np.asarray(x_new)
        out_of_bounds = self._check_bounds(x_new)
        y_new = self._call(x_new)
        if len(y_new) > 0:
            y_new[out_of_bounds] = self.fill_value
        return y_new

    def _call(self,x_new):
        x_new_indices = np.searchsorted(self.x,x_new)
        x_new_indices = x_new_indices.clip(1, len(self.x)-1).astype(int)
        lo = x_new_indices - 1
        hi = x_new_indices
        x_lo = self.x[lo]
        x_hi = self.x[hi]
        y_lo = self._y[lo]
        y_hi = self._y[hi]
        slope = (y_hi-y_lo) / (x_hi-x_lo)#[:,None]
        return slope*(x_new-x_lo) + y_lo
        #return slope*(x_new-x_lo)[:,None] + y_lo
    
    def _check_bounds(self,x_new):
        below_bounds = x_new < self.x[0]
        above_bounds = x_new > self.x[-1]
        if self.bounds_error and below_bounds.any():
            raise ValueError("A value in x_new is below the interpolation "
                "range.")
        if self.bounds_error and above_bounds.any():
            raise ValueError("A value in x_new is above the interpolation "
                "range.")
        out_of_bounds = np.logical_or(below_bounds, above_bounds)
        return out_of_bounds
    
    def __call__(self,x):
        x = np.asarray(x)
        if not np.issubdtype(x.dtype,np.inexact):
            x = x.astype(float)
        x = x.ravel() #x_shape = x.shape
        return self._evaluate(x)
