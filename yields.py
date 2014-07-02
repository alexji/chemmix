import numpy as np
import pylab as plt
import os
import karlsson

import hw10

def map_elemnames_to_elemmass(elemnames):
    elemmap = {'C': 12.0, 'O': 16.0,
               'Mg': 24.3, 'Si': 28.1,
               'Ca': 40.1, 'Fe': 55.8}
    return np.array([elemmap[elem] for elem in elemnames])
def map_elemnames_to_asplund09(elemnames):
    """ Solar log10(nZ/nH) """
    elemmap = {'C': 8.43, 'O': 8.69,
               'Mg': 7.60, 'Si': 7.51,
               'Ca': 6.34, 'Fe': 7.50}
    return np.array([elemmap[elem]-12.0 for elem in elemnames])

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

class randombase(yieldsbase):
    def draw_yields(self,N):
        sn_types = karlsson.draw_from_distr(N,self.sntypearr,self.sntypepdf)
        return self.get_yields(sn_types)

class ryI05N06(randombase):
    def __init__(self,p):
        self.p = p
        self.sntypearr = np.array([1,2])
        self.sntypepdf = np.array([p,1-p])
        self.name = "Iwamoto05 and Nomoto06 Yields (p=%f)" % (p)
        self.shortname = "I05N06p%0.2f" % (p)
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
        self.numtypes  = 2
    def yieldfn(self,sn_type):
        if sn_type==1: #I05
            return [.20, 1.1e-5,.39,3e-4,3.6e-5,2.7e-6]
        if sn_type==2: #N06 M=20
            return [ 0.21100001,  2.11000009,  0.150354,    0.099692,    0.00624737,  0.072287  ]

class ryN06(randombase):
    def __init__(self,imf):
        assert len(imf)==7
        assert np.sum(imf)-1 < 10**-10
        self.sntypepdf = np.array(imf)
        self.sntypearr = np.array([1,2,3,4,5,6,7])
        
        self.name = "Nomoto 06 Yields"
        self.shortname = "N06"
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
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

class I05N06yields(yieldsbase):
    def __init__(self):
        self.name = "Iwamoto05 and Nomoto06 Yields"
        self.shortname = "I05N06"
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
        self.numtypes  = 2
    def yieldfn(self,sn_type):
        if sn_type==1: #I05
            return [.20, 1.1e-5,.39,3e-4,3.6e-5,2.7e-6]
        if sn_type==2: #N06 M=20
            return [ 0.21100001,  2.11000009,  0.150354,    0.099692,    0.00624737,  0.072287  ]

class nomoto06yields(yieldsbase):
    def __init__(self):
        self.name = "Nomoto 06 Yields"
        self.shortname = "N06"
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
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
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
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
    def __init__(self,p,f,
                 datafile="YIELDDATA/interp_nomoto06_z0.npy",autointerpolate=True):
        super(nomoto06interpyields_Cenhance,self).__init__(datafile=datafile,autointerpolate=autointerpolate)
        self.p = p; self.f = float(f); assert p >= 0 and p <= 1
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
        assert np.abs(np.sum(newpdf) - 1.0) < 3*np.finfo(float).eps
        return newpdf

class nomoto06interpyields_Cadd(nomoto06interpyields):
    def __init__(self,p,MC,
                 datafile="YIELDDATA/interp_nomoto06_z0.npy",autointerpolate=True):
        super(nomoto06interpyields_Cadd,self).__init__(datafile=datafile,autointerpolate=autointerpolate)
        self.p = p; self.MC = float(MC); assert p >= 0 and p <= 1
        self.name += (" Carbon enhanced p=%0.3f MC=%f" % (self.p, self.MC))
        self.shortname += ("p%0.3fMC%f" % (self.p, self.MC))
        self.numtypes *= 2
        self.yieldarr = np.concatenate([self.yieldarr,self.yieldarr])
        assert self.yieldarr.shape == (self.numtypes,self.numyields)
        self.yieldarr[0:(self.numtypes/2),0] += self.MC #enhance carbon in the first half
        print "--Cadd: use modify_sntypepdf(sntypepdf) as the sntypepdf"
        
    def modify_sntypepdf(self,sntypepdf):
        sntypepdf /= np.sum(sntypepdf) #normalize
        newpdf = np.concatenate([self.p*sntypepdf,(1-self.p)*sntypepdf])
        assert np.abs(np.sum(newpdf) - 1.0) < 3*np.finfo(float).eps
        return newpdf

class hw10yields(yieldsbase):
    def __init__(self,E,cut,mix):
        self.datafile = hw10.get_hw10_filename(E,cut,mix)
        self.E = E; self.cut = cut; self.mix = mix
        self.name = "HW10 Yields "+hw10.get_hw10_label(E,cut,mix)
        self.shortname = "HW10E"+str(E)+self.cut+"m"+str(mix)
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
        self.numtypes  = 31
        self.massarr = np.arange(10,40+1)
        self.yieldarr = hw10.load_hw10(E,cut,mix)
        
    def get_yields(self,sn_type):
        sn_type = self.clean_sn_type_input(sn_type).reshape(-1)
        indices = np.array(sn_type)-1
        return self.yieldarr[indices,:]

class mixN06HW10yields(nomoto06interpyields,hw10yields):
    def __init__(self,p,E,cut,mix,
                 autointerpolate=True):
        self.p = p; assert p >= 0 and p <= 1
        self.hw10datafile = hw10.get_hw10_filename(E,cut,mix)
        self.E = E; self.cut = cut; self.mix = mix
        self.n06datafile  = "YIELDDATA/interp_nomoto06_z0.npy"
        if not os.path.exists(self.n06datafile) and autointerpolate:
            super(mixN06HW10yields,self).create_interp_nomoto06_yields()
        self.n06yieldarr = np.load(self.n06datafile)
        self.hw10yieldarr= hw10.load_hw10(E,cut,mix)
        self.yieldarr = np.concatenate([self.n06yieldarr,self.hw10yieldarr[3:31,:]])
        
        self.name = "N06i and HW10 combined with p="+str(p)
        self.shortname = "mixN06HW10p"+str(p)
        self.numyields = 6
        self.elemnames = ['C','O','Mg','Si','Ca','Fe']
        self.elemmass = map_elemnames_to_elemmass(self.elemnames)
        self.numtypes  = 56 #M = 13 to 40 Msun, inclusive x 2
        self.massarr = np.arange(13.,40+1)
        print "--mix: p=%f to have N06" % (p)
        print "--mix: use modify_sntypepdf(sntypepdf) as the sntypepdf"
        
    def modify_sntypepdf(self,sntypepdf):
        sntypepdf /= np.sum(sntypepdf) #normalize
        newpdf = np.concatenate([self.p*sntypepdf,(1-self.p)*sntypepdf])
        assert np.abs(np.sum(newpdf) - 1.0) < 3*np.finfo(float).eps
        return newpdf

def plot_yields(yieldobj,subplota,subplotb):
    sntypearr = np.arange(1,yieldobj.numtypes+1)
    yieldarr = yieldobj(sntypearr)
    if yieldobj.numtypes == len(yieldobj.massarr):
        x = yieldobj.massarr
        xlab = 'SN Mass'
    else:
        x = np.arange(1,yieldobj.numtypes+1)
        xlab = 'sn type'

    plt.figure()
    for ploti in xrange(yieldobj.numyields):
        plt.subplot(subplota, subplotb, ploti+1)
        plt.plot(x,yieldarr[:,ploti])
        plt.xlabel(xlab)
        plt.ylabel(yieldobj.elemnames[ploti])
    plt.tight_layout()
    plt.show()

def plot_yield_ratios(yieldobj,Fe_ix,subplota,subplotb):
    sntypearr = np.arange(1,yieldobj.numtypes+1)
    yieldarr = yieldobj(sntypearr)
    
    asplund09fe = map_elemnames_to_asplund09(yieldobj.elemnames) +12-7.50
    metal_num = yieldarr/yieldobj.elemmass
    metal_num = (metal_num.transpose()/metal_num[:,5]).transpose()
    yieldarr = np.log10(metal_num) - asplund09fe
    
    if yieldobj.numtypes == len(yieldobj.massarr):
        x = yieldobj.massarr
        xlab = 'SN Mass'
    else:
        x = np.arange(1,yieldobj.numtypes+1)
        xlab = 'sn type'

    plt.figure()
    for ploti in xrange(yieldobj.numyields):
        plt.subplot(subplota, subplotb, ploti+1)
        plt.plot(x,yieldarr[:,ploti])
        plt.xlabel(xlab)
        plt.ylabel('['+str(yieldobj.elemnames[ploti])+'/Fe]')
    plt.tight_layout()
    plt.show()
