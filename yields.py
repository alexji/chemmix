### IMPORTANT: all yields functions:
### If they receive the argument YIELD_LENGTH,
###  should return the number of yields they output
### If they receive the argument NUM_SN_TYPES,
###  should return the number of SN types

import numpy as np
import os

class yieldsbase(object):
    """ Base functions for yield classes """
    def clean_sn_type_input(self,sn_type):
        try:
            arrlen = len(sn_type)
            sn_type = np.array(sn_type)
        except TypeError:
            arrlen = 1
            sn_type = np.array([[sn_type]])
        assert np.all(sn_type >= 1) & np.all(sn_type <= self.numtypes)
        return sn_type
    def get_yields(self,sn_type):
        sn_type = self.clean_sn_type_input(sn_type)
        yields = [self.yieldfn(this_sn_type) for this_sn_type in sn_type]
        return np.array(yields)
    def __call__(self,sn_type):
        return self.get_yields(sn_type)
    def __repr__(self):
        return self.name
    def map_elemnames_to_elemmass(self,elemnames):
        elemmap = {'C': 12.0, 'O': 16.0,
                   'Mg': 24.3, 'Si': 28.1,
                   'Ca': 40.1, 'Fe': 55.8}
        return np.array([elemmap[elem] for elem in elemnames])
    def map_elemnames_to_asplund09(self,elemnames):
        """ Solar log10(nZ/nH) """
        elemmap = {'C': 8.43, 'O': 8.69,
                   'Mg': 7.60, 'Si': 7.51,
                   'Ca': 6.34, 'Fe': 7.50}
        return np.array([elemmap[elem]-12.0 for elem in elemnames])

class nomoto06yields(yieldsbase):
    def __init__(self):
        self.name = "Nomoto 06 Yields"
        self.shortname = "N06"
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = self.map_elemnames_to_elemmass(self.elemnames)
        self.numtypes  = 7
        self.massarr = np.array([13.,15,18,20,25,30,40])
    def yieldfn(self,sn_type):
        if sn_type==1: #13 Msun
            return [ 0.07410008,  0.45000175,  0.0864267,   0.08257,     0.00293784,  0.071726  ]
        if sn_type==2: #15 Msun
            return [ 0.17200006,  0.77300646,  0.068896,    0.073588,    0.00443338,  0.07238   ]
        if sn_type==3: #18 Msun
            return [ 0.218,       1.38000491,  0.158456,    0.116787,    0.00441815,  0.072278  ]
        if sn_type==4: #20 Msun
            return [ 0.21100001,  2.11000009,  0.150354,    0.099692,    0.00624737,  0.072287  ]
        if sn_type==5: #25 Msun
            return [ 0.29400001,  2.79000216,  0.1200898,   0.3513464,   0.02481727,  0.073777  ]
        if sn_type==6: #30 Msun
            return [ 0.33700001,  4.81000002,  0.226373,    0.248843,    0.0174063,   0.074573  ]
        if sn_type==7: #40 Msun
            return [ 0.429,       8.38000021,  0.478553,    1.02666,     0.03733042,  0.08000101]


class nomoto06interpyields(yieldsbase):
    def __init__(self,datafile="YIELDDATA/interp_nomoto06_z0.npy",autointerpolate=True):
        self.name = "Nomoto 06 Yields (Interpolated)"
        self.shortname = "N06i"
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = self.map_elemnames_to_elemmass(self.elemnames)
        self.numtypes  = 28 #M = 13 to 40 Msun, inclusive
        self.massarr = np.arange(13.,40+1)

        self.datafile = datafile
        if not os.path.exists(datafile) and autointerpolate:
            self.create_interp_nomoto06_yields()
        self.yieldarr = np.load(datafile)
        #self.yieldfn = None #directly index the array, redefining get_yields
        
    def get_yields(self,sn_type):
        sn_type = self.clean_sn_type_input(sn_type).reshape(-1)
        indices = np.array(sn_type)-1
        return self.yieldarr[indices,:]
    
    def create_interp_nomoto06_yields(self):
        print "Interpolating nomoto06 and saving into",self.datafile
        from scipy.interpolate import interp1d
        Minput = [13,15,18,20,25,30,40]
        Moutput = np.arange(13,41,1)
        nomoto06arr = np.zeros((7,6))
        outputarr = np.zeros((len(Moutput),6))
        n06y = nomoto06yields()
        for i in xrange(7):
            nomoto06arr[i,:] = n06y(i+1)
        for j in xrange(6):
            f = interp1d(Minput,nomoto06arr[:,j])
            outputarr[:,j] = f(Moutput)
        np.save(self.datafile,outputarr)
        #import pylab as plt
        #plt.figure(figsize=(8,10))
        #elemarr = ['C','O','Mg','Si','Ca','Fe']
        #for j,elem in enumerate(elemarr):
        #    plt.subplot(3,2,j+1)
        #    plt.plot(Moutput,outputarr[:,j],'o-')
        #    plt.xlabel('M'); plt.ylabel(elemarr[j])
        #plt.savefig("PLOTS/interp_nomoto06_z0.png")

class nomoto06interpyields_Cenhance(nomoto06interpyields):
    def __init__(self,p=0.1,f=100.,
                 datafile="YIELDDATA/interp_nomoto06_z0.npy",autointerpolate=True):
        super(nomoto06interpyields_Cenhance,self).__init__(datafile=datafile,autointerpolate=autointerpolate)
        self.p = p; self.f = float(f)
        self.name += (" Carbon enhanced p=%0.3f f=%0.1f" % (self.p, self.f))
        self.shortname += ("p%0.3ff%0.1f" % (self.p, self.f))
        self.numtypes *= 2
        self.yieldarr = np.concatenate([self.yieldarr,self.yieldarr])
        assert self.yieldarr.shape == (self.numtypes,self.numyields)
        self.yieldarr[0:(self.numtypes/2),0] *= self.f #enhance carbon in the first half
        print "--Cenhance: use modify_sntypepdf(sntypepdf) as the sntypepdf"
        
    def modify_sntypepdf(self,sntypepdf):
        sntypepdf /= np.sum(sntypepdf) #normalize
        newpdf = np.concatenate([self.p*sntypepdf,(1-self.p)*sntypepdf])
        print newpdf,newpdf.shape,np.sum(newpdf)
        assert np.abs(np.sum(newpdf) - 1.0) < 3*np.finfo(float).eps
        return newpdf

#def _simple_yields(sn_type):
#    if sn_type==1:
#        return [0.5,0.6]
#    if sn_type==2:
#        return [0.3,0.9]
#    raise ValueError("sn_type invalid")
#def simple_yields(sn_type):
#    if return_yield_length(sn_type):
#        return 2
#    if return_num_sn_types(sn_type):
#        return 2
#    try:
#        yields = []
#        for this_sn_type in sn_type:
#            yields.append(_simple_yields(this_sn_type))
#        yields = np.array(yields)
#        return yields
#    except TypeError:
#        yields = np.array(_testyield(sn_type))
#        return np.reshape(np.array(yields),(1,-1))
#
#
#def interp_nomoto06_yields(sn_type):
#    """
#    sn_type 1 to 28 <=> mass from 13 to 40 Msun in integer increments
#    """
#    if return_yield_length(sn_type):
#        return 6
#    if return_num_sn_types(sn_type):
#        return 28
#    indices = np.array(sn_type)-1
#    #enforce sn_type:
#    assert np.all(np.logical_and(indices >= 0, indices <= 27))
#    yields = np.load("DATA/interp_nomoto06_z0.npy")
#    return yields[indices,:]
#
#def convert_nomoto06_yields(yieldarr,boxsize,cellres):
#    """
#    Converts a nomoto06 yieldarr into [X/Fe]
#    """
#    #yieldarr is shape N x 6, the 6th column is Fe
#    assert yieldarr.shape[1] == 6
#    molecmass = [12.,16.,24.3,28.1,40.1,55.8]
#    asplund09 = np.array([8.43,8.69,7.60,7.51,6.34,7.50])
#    asplund09fe = asplund09-7.50
#
#    metal_num = yieldarr/molecmass
#    metal_ratios = np.transpose(np.transpose(metal_num)/metal_num[:,5])
#    abundances = np.log10(metal_ratios) - asplund09fe
#    posnres = boxsize/float(cellres)
#    #print posnres
#    # convert iron mass into n_Fe/n_H (particles/cc, assuming n_H = 0.1/cc), then to [Fe/H]
#    #print yieldarr[:,5]
#    abundances[:,5] = np.log10(yieldarr[:,5]/(posnres**3 * 55.8) * 4.08 * 10**-8 * 10.) - (7.5 - 12.0)
#    return abundances

