##########################
# tophat.py
# Gives the solution to a spherical tophat collapse model
# Default constants chosen for Planck cosmology (h=.67, Om=.32, Ol=.68)
# with Mhalo = 10^8 Msun, nvir = 0.1/cc (simple atomic cooling halo with zvir ~ 10)

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy.integrate import quad

def cosmo_age(z,Ol=0.7,Om=0.3,h0=0.7):
    H0inv = 9778./h0 #Myr
    return H0inv * quad(lambda zp: 1./((1+zp)*np.sqrt(Om*(1+zp)**3 + Ol)), z, np.inf)[0]

class TopHat:
    def __init__(self,Mhalo=10**8,zvir=10,fb=0.1551,mu=1.4,
                 critfactor=200,
                 Ol = 0.7, Om = 0.3, h0 = 0.7,
                 verbose=False):
        """
        @param Mhalo: mass of DM halo in solar masses (default 10^8)
        @param zvir: redshift of virialization
        @param fb: baryon fraction (default Planck)
        @param mu: gas effective mass (default 1.4)
        """
        self.Mhalo = Mhalo
        self.zvir = zvir
        self.fb = fb
        self.mu = mu
        self.critfactor = critfactor
        self.Ol = Ol; self.Om = Om; self.h0 = h0

        # Assume rho_vir = critfactor * rho_crit
        nc0 = h0**2 * 1.123*10**-5 #protons/cc
        nvir = fb*critfactor*nc0*(Om*(1+zvir)**3 + Ol)/mu #particles/cc
        self.nvir = nvir

        # calculate cycloid A and B; A in kpc, B in Myr
        # A**3 = Rvir**3 = (3M/(4 pi rho_vir)) = (3 M fb / (4 pi mu mp nvir))
        self.A = (2.838 * 10**56 * Mhalo * fb/(mu * nvir))**(1./3) / (3.086 * 10**21) #kpc
        self.Rvir = self.A
        # B**2 = A**3/GM = (3 fb)/(4 pi G mu mp nvir)
        self.B = 46.37 * np.sqrt(fb/(mu*nvir)) #Myr
        self.tvir = cosmo_age(zvir,Ol,Om,h0)
        # vvir**2 = GM/Rvir
        self.vvir = .002074 * np.sqrt(Mhalo/self.Rvir) #km/s

        self.verbose = verbose
        if verbose:
            print "TopHat: A =",self.A,"kpc, B =",self.B,"Myr"

    def cycloidR(self,eta):
        return self.A*(1 - np.cos(eta))
    def cycloidt(self,eta):
        return self.B*(eta - np.sin(eta))

    def n_of_eta(self,eta):
        #try:
        #    assert eta >= 0 and eta <= 2*np.pi, 'eta: '+str(eta)
        #except ValueError:
        #    assert (eta >= 0).all() and (eta <= 2*np.pi).all(), 'eta:' +str(eta)
        try:
            narr = np.zeros(len(eta)); eta = np.array(eta)
            ix = eta < 3*np.pi/2
            narr[ix]  = self.nvir / (1 - np.cos(eta[ix]))**3
            narr[~ix] = self.nvir
            return narr
        except ValueError:
            if eta < 3*np.pi/2:
                return self.nvir / (1-np.cos(eta))**3
            else:
                return self.nvir

    def get_eta_from_t(self,t):
        try:
            return brentq(lambda eta: t-self.cycloidt(eta),0,2*np.pi)
        except TypeError:
            return [brentq(lambda eta: tt-self.cycloidt(eta),0,2*np.pi) for tt in t]
        except ValueError:
            print "Error: t out of bounds for this B. (t,B) = ("+str(t)+","+str(self.B)+")"
            raise
    def n_of_t(self,t):
        """ return density in cm^-3 """
        A = self.A; B = self.B
        tmax  = self.cycloidt(2*np.pi)
        tflat = self.cycloidt(3*np.pi/2)
        try:
            assert t >= 0# and t <= tmax
        except ValueError:
            assert (t >= 0).all()# and (t <= tmax).all()
        eta = self.get_eta_from_t(t)
        return self.n_of_eta(eta)
    def get_n_of_t_fn(self,Nres=300):
        """ return interpolated density fn in cm^-3 """
        tmax = self.B*2*np.pi
        #tarr = np.concatenate((np.logspace(-4,np.log10(tmax/Nres),Nres/2),np.linspace(tmax/Nres,tmax,Nres/2)))
        tarr = np.logspace(-4,np.log10(tmax),Nres,endpoint=False)
        narr = self.n_of_t(tarr)
        tarr = np.concatenate((tarr,[100000.])); narr = np.concatenate((narr,[self.nvir]))
        if self.verbose:
            print "minimum t:",tarr[0]
        return interp1d(tarr,narr)
    def rho_of_t(self,t):
        """ return mass density in g/cm^3 """
        # proton mass in grams = 1.673 * 10^-24
        return self.mu * 1.673*10**-24 * self.n_of_t(t)
    def get_rho_of_t_fn(self,Nres=300):
        tmax = self.B*2*np.pi
        #tarr = np.concatenate((np.logspace(-4,np.log10(tmax/Nres),Nres/2),np.linspace(tmax/Nres,tmax,Nres/2)))
        tarr = np.logspace(-4,np.log10(tmax),Nres,endpoint=False)
        rhoarr = self.rho_of_t(tarr)
        tarr = np.concatenate((tarr,[100000.])); rhoarr = np.concatenate((rhoarr,[self.mu*1.678*(10**-24)*self.nvir]))
        if self.verbose:
            print "minimum t:",tarr[0]
        return interp1d(tarr,rhoarr)

    ## TODO cosmology dependent n_of_z
    #def n_of_z(self,z):
    #    pass

if __name__=="__main__":
    th = TopHat(verbose=True)
    n = th.get_n_of_t_fn()
    tarr = np.linspace(10,700,300)
    import pylab as plt
    plt.plot(tarr,n(tarr))
    plt.gca().set_yscale('log')
    plt.show()
