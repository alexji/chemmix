import numpy as np
import h5py
import pickle
import sfrmodel as sfu
from tophat import TopHat
import scipy.special

import backgroundmodel

def load_Vmix(envname,get_th=False,get_thrho=False,getVmixdot=False,getarrays=False):
    Mhalo,zvir,vturb,lturb,nSN,trecovery,Dt,uSN,Vmax = backgroundmodel.envparams(envname)
    if getarrays: myfn = backgroundmodel.get_Vmixarr
    else: myfn = backgroundmodel.get_Vmixfn

    tarr = backgroundmodel.gettarr(envname,hires=True)
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    rhofn = th.get_rho_of_t_fn()
    if getVmixdot:
        Vmix,Vmixdot = myfn(tarr,th,Dt,Vmax,getVmixdot=True)
        if get_thrho: return th,rhofn,Vmix,Vmixdot
        if get_th: return th,Vmix,Vmixdot
        return Vmix,Vmixdot
    else:
        Vmix = myfn(tarr,th,Dt,Vmax)
        if get_thrho: return th,rhofn,Vmix
        if get_th: return th,Vmix
        return Vmix

def load_mmixgrid(envname):
    return np.load(fname_mmixgrid(envname))

def load_sfr(envname,sfrname):
    label = sfu.get_sfrlabel(envname,sfrname)
    uII   = sfu.loaduIIfn(label)
    uIII  = sfu.loaduIIIfn(label)
    muII  = sfu.loadmuIIfn(label)
    muIII = sfu.loadmuIIIfn(label)
    wII   = poisson
    wIII  = modifiedpoisson
    return uII,uIII,muII,muIII,wII,wIII

def load_ckk(envname,sfrname):
    filename_ckk = fname_ckk(envname,sfrname)
    f = h5py.File(filename_ckk,'r')
    k2max = f.attrs['k2max']
    k3max = f.attrs['k3max']
    cgrid = np.array(f['cgrid'])
    f.close()
    return k2max,k3max,cgrid

def load_histlist(envname,sfrname,postfix):
    f = open(fname_chemgridhistlist(envname,sfrname,postfix),'r')
    mylist = pickle.load(f)
    f.close()
    return mylist

def _poisson(k,mu):
    return np.exp(k*np.log(mu)-scipy.special.gammaln(k+1)-mu)
def _poisson_mu0(k):
    if k==0: return 1
    return 0
def poisson(k,mu):
    try: #mu is an array; assume k is never an array
        output = np.empty(len(mu))
        imu0 = (mu==0)
        output[imu0]  = _poisson_mu0(k)
        output[~imu0] = _poisson(k,mu[~imu0])
        return output
    except: #mu is a scalar
        if mu==0: return _poisson_mu0(k)
        return _poisson(k,mu)
def modifiedpoisson(k,mu):
    """ NOT normalized!"""
    return np.exp(-(mu-k)**2/4.) * poisson(k,mu) #np.exp(-mu) * mu**k/factorial(k) 

def relative_imf(marr,alpha,norm=False):
    Mmax = float(marr[-1])
    imfpdf = np.array([(M/Mmax)**(-alpha) for M in marr])
    if norm:
        imfpdf = imfpdf/np.sum(imfpdf)
    return imfpdf

def hist2d2pdf(H,xedges,yedges):
    nbins_x = len(xedges)-1; nbins_y = len(yedges)-1
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    pdf = (H*(x_bin_sizes*y_bin_sizes).T)
    return pdf


def default_filename(envname,sfrname,postfix=''):
    return envname+'_'+sfrname+postfix
def logMbinsMplot(logMmin,logMmax,logdM):
    Mbins = 10.**np.arange(logMmin,logMmax+logdM/2.,logdM)
    Mplot = 10.**(np.arange(logMmin,logMmax,logdM)+logdM/2.)
    return Mbins,Mplot


def fname_mmixgrid(envname):
    return "MMIXGRIDS/"+envname+"_mmixgrid.npy"
def fname_ckk(envname,sfrname):
    return 'MMIXDISTR/'+envname+'_'+sfrname+'_ckk.hdf5'
def fname_fMkkII(envname,sfrname):
    return "MMIXDISTR/"+envname+'_'+sfrname+'_fMkkII.npy'
def fname_fMkkIII(envname,sfrname):
    return "MMIXDISTR/"+envname+'_'+sfrname+'_fMkkIII.npy'
def fname_Mbinskk(envname,sfrname):
    return "MMIXDISTR/"+envname+'_'+sfrname+'_Mbins.npy'
def fname_Mplotkk(envname,sfrname):
    return "MMIXDISTR/"+envname+'_'+sfrname+'_Mplot.npy'
def fname_chemgridhistlist(envname,sfrname,postfix):
    return "CHEMGRIDS/"+envname+'_'+sfrname+'_'+postfix+'.list'
