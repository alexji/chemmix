import numpy as np
import asciitable
import re
import os
import gc
import pylab as plt

def get_hw10_label(E,cut,mix):
    return "E"+str(E)+" "+cut+" mix"+str(mix)

def get_hw10_filename(E,cut,mix,name='sixelem'):
    assert E in [0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,3,5,10]
    assert cut in ['S4','Ye']
    assert mix in [0,.001,.00158,.00251]
    return "DATA/HW10_"+name+"_E"+str(E)+"_C"+cut+"_mix"+str(mix)+".npy"

def load_hw10(E,cut,mix,name='sixelem'):
    filename = get_hw10_filename(E,cut,mix,name=name)
    if os.path.exists(filename):
        print "Loading "+filename
        return np.load(filename)
    elif name=='sixelem':
        print "Reducing full table "+filename+" (may take a while and/or overflow memory)..."
        print 
        process_hw10_table_sixelem(E,cut,mix)
        gc.collect() #garbage collect large data
        return np.load(filename)
    else:
        raise IOError("No HW10 file available or computable")

def process_hw10_table_sixelem(E,cut,mix):
    filename=get_hw10_filename(E,cut,mix,name='sixelem')
    tab = asciitable.read("DATA/hegerwoosley10.dat")
    #tab.dtype.names = 'Mass','Energy','Cut','Mixing','Isotope','Yield'
    elemnames = ['C','O','Mg','Si','Ca','Fe']
    elemarr = [re.compile(elem+r'[0-9]+\Z') for elem in elemnames]
    nelem = len(elemarr)
    massarr = np.arange(10,41)
    nmass = len(massarr)

    outtab = np.zeros((nmass,nelem))
    for row in tab:
        if row[0] in massarr and row[1]==E and row[2]==cut and row[3]==mix:
            i = int(row[0]-10)
        else:
            continue
        for j,iselem in enumerate(elemarr):
            if iselem.search(row[4]):
                outtab[i,j] += row[5]
                break
    np.save(filename,outtab)
    
def plot_hw10_sixelem(E,cut,mix,ratio=False,**kwargs):
    tab = load_hw10(E,cut,mix)
    Marr = np.arange(10,41)
    elemnames = ['C','O','Mg','Si','Ca','Fe']
    molecmass = np.array([12.,16.,24.3,28.1,40.1,55.8])
    asplund09 = np.array([8.43,8.69,7.60,7.51,6.34,7.50])
    asplund09fe = asplund09-7.50

    if ratio: #plot [X/Fe] ratio
        metal_num = tab/molecmass
        metal_num = (metal_num.transpose()/metal_num[:,5]).transpose() #number ratio over iron
        metal_abundances = np.log10(metal_num) - asplund09fe #[X/Fe]

    for i in range(6):
        plt.subplot(3,2,i+1)
        if ratio:
            plt.plot(Marr,metal_abundances[:,i],**kwargs)
            plt.ylabel('['+elemnames[i]+'/Fe]')
        else:
            plt.plot(Marr,tab[:,i],**kwargs)
            plt.ylabel(elemnames[i])
    return tab

if __name__=="__main__":
    ratioflag = False

    plt.figure()
    nomoto = np.load("DATA/interp_nomoto06_z0.npy")
    if ratioflag:
        molecmass = np.array([12.,16.,24.3,28.1,40.1,55.8])
        asplund09 = np.array([8.43,8.69,7.60,7.51,6.34,7.50])
        asplund09fe = asplund09-7.50
        metal_num = nomoto/molecmass
        metal_num = (metal_num.transpose()/metal_num[:,5]).transpose() #number ratio over iron
        nomoto = np.log10(metal_num) - asplund09fe #[X/Fe]

    elemnames = ['C','O','Mg','Si','Ca','Fe']
    for i in range(6):
        plt.subplot(3,2,i+1)
        plt.plot(np.arange(13,41),nomoto[:,i],'k',lw=2)
        plt.ylabel(elemnames[i])
    plot_hw10_sixelem(0.9,'S4',0,     ratio=ratioflag,color='blue',  label='E0.9')
    plot_hw10_sixelem(0.9,'S4',.00251,ratio=ratioflag,color='blue',  ls='dashed',lw=2)
    plot_hw10_sixelem(1.2,'S4',0,     ratio=ratioflag,color='green', label='E1.2')
    plot_hw10_sixelem(1.2,'S4',.00251,ratio=ratioflag,color='green', ls='dashed',lw=2)
    plot_hw10_sixelem(1.5,'S4',0,     ratio=ratioflag,color='red',   label='E1.5')
    plot_hw10_sixelem(1.5,'S4',.00251,ratio=ratioflag,color='red',   ls='dashed',lw=2)
    plot_hw10_sixelem(3.0,'S4',0,     ratio=ratioflag,color='cyan',  label='E3.0')
    plot_hw10_sixelem(3.0,'S4',.00251,ratio=ratioflag,color='cyan',  ls='dashed',lw=2)
    plot_hw10_sixelem(10, 'S4',0,     ratio=ratioflag,color='purple',label='E10')
    plot_hw10_sixelem(10, 'S4',.00251,ratio=ratioflag,color='purple',ls='dashed',lw=2)
    
    plot_hw10_sixelem(1.2,'Ye',0,     ratio=ratioflag,color='green',ls='dotted')
    plot_hw10_sixelem(10,'Ye',0,      ratio=ratioflag,color='purple',ls='dotted')
    
    leg = plt.legend(loc='best')
    ltext = leg.get_texts(); llines = leg.get_lines()
    plt.setp(ltext,fontsize='xx-small'); plt.setp(llines, lw=1.0)
    if not ratioflag:
        for i in range(6):
            plt.subplot(3,2,i+1)
            plt.gca().set_yscale('log')
            plt.ylim((1e-5,10))
    else:
        for i in range(6):
            plt.subplot(3,2,i+1)
            #plt.ylim((-1,1.5))
    plt.show()
