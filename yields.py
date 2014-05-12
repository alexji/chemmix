### IMPORTANT: all yields functions:
### If they receive the argument YIELD_LENGTH,
###  should return the number of yields they output
### If they receive the argument NUM_SN_TYPES,
###  should return the number of SN types

import numpy as np

YIELD_LENGTH = -9999
NUM_SN_TYPES = -99999
def return_yield_length(sn_type):
    try:
        if sn_type == YIELD_LENGTH:
            return True
        else:
            return False
    except ValueError: #sn_type is an array
        return False
def return_num_sn_types(sn_type):
    try:
        if sn_type == NUM_SN_TYPES:
            return True
        else:
            return False
    except ValueError: #sn_type is an array
        return False

def _simple_yields(sn_type):
    if sn_type==1:
        return [0.5,0.6]
    if sn_type==2:
        return [0.3,0.9]
    raise ValueError("sn_type invalid")
def simple_yields(sn_type):
    if return_yield_length(sn_type):
        return 2
    if return_num_sn_types(sn_type):
        return 2
    try:
        yields = []
        for this_sn_type in sn_type:
            yields.append(_simple_yields(this_sn_type))
        yields = np.array(yields)
        return yields
    except TypeError:
        yields = np.array(_testyield(sn_type))
        return np.reshape(np.array(yields),(1,-1))

def _nomoto06_yields(sn_type):
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
def nomoto06_yields(sn_type):
    """
    Type 1: 13 Msun
    Type 2: 15 Msun
    Type 3: 18 Msun
    Type 4: 20 Msun
    Type 5: 25 Msun
    Type 6: 30 Msun
    Type 7: 40 Msun

    Six yields: C, O, Mg, Si, Ca, Fe (in Msun)
    """
    if return_yield_length(sn_type):
        return 6
    if return_num_sn_types(sn_type):
        return 7
    try:
        yields = []
        for this_sn_type in sn_type:
            yields.append(_nomoto06_yields(this_sn_type))
        yields = np.array(yields)
        return yields
    except TypeError:
        yields = np.array(_nomoto06_yields(sn_type))
        return np.reshape(np.array(yields),(1,-1))

def convert_nomoto06_yields(yieldarr,boxsize,cellres):
    """
    Converts a nomoto06 yieldarr into [X/Fe]
    """
    #yieldarr is shape N x 6, the 6th column is Fe
    assert yieldarr.shape[1] == 6
    molecmass = [12.,16.,24.3,28.1,40.1,55.8]
    asplund09 = np.array([8.43,8.69,7.60,7.51,6.34,7.50])
    asplund09fe = asplund09-7.50

    metal_num = yieldarr/molecmass
    metal_ratios = np.transpose(np.transpose(metal_num)/metal_num[:,5])
    abundances = np.log10(metal_ratios) - asplund09fe
    posnres = boxsize/float(cellres)
    #print posnres
    # convert iron mass into n_Fe/n_H (particles/cc, assuming n_H = 0.1/cc), then to [Fe/H]
    #print yieldarr[:,5]
    abundances[:,5] = np.log10(yieldarr[:,5]/(posnres**3 * 55.8) * 4.08 * 10**-8 * 10.) - (7.5 - 12.0)
    return abundances


def interp_nomoto06_yields(sn_type):
    """
    sn_type 1 to 28 <=> mass from 13 to 40 Msun in integer increments
    """
    if return_yield_length(sn_type):
        return 6
    if return_num_sn_types(sn_type):
        return 28
    indices = np.array(sn_type)-1
    #enforce sn_type:
    assert np.all(np.logical_and(indices >= 0, indices <= 27))
    yields = np.load("DATA/interp_nomoto06_z0.npy")
    return yields[indices,:]

def create_interp_nomoto06_yields():
    from scipy.interpolate import interp1d
    Minput = [13,15,18,20,25,30,40]
    Moutput = np.arange(13,41,1)
    nomoto06arr = np.zeros((7,6))
    outputarr = np.zeros((len(Moutput),6))
    for i in xrange(7):
        nomoto06arr[i,:] = _nomoto06_yields(i+1)
    for j in xrange(6):
        f = interp1d(Minput,nomoto06arr[:,j])
        outputarr[:,j] = f(Moutput)
    np.save("DATA/interp_nomoto06_z0",outputarr)
    import pylab as plt
    plt.figure(figsize=(8,10))
    elemarr = ['C','O','Mg','Si','Ca','Fe']
    for j,elem in enumerate(elemarr):
        plt.subplot(3,2,j+1)
        plt.plot(Moutput,outputarr[:,j],'o-')
        plt.xlabel('M'); plt.ylabel(elemarr[j])
    plt.savefig("PLOTS/interp_nomoto06_z0.png")
