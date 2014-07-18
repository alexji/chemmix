import numpy as np
import pylab as plt

import karlsson
import util
import yields

def draw_reshape_yields(y,k,N,fM,Mplot):
    if N==0: return np.zeros((0,k,y.numyields))
    m = karlsson.draw_from_distr(k*N,Mplot,fM)
    a = y.draw_yields(k*N)
    return np.reshape((a.T/m).T,(N,k,y.numyields))

if __name__=="__main__":
    envname='atomiccoolinghalo'
    sfrname='fidTS500'
    fMkkp = np.load(karlsson.get_fMkkp_filename(envname,sfrname))
    Mmax = 10**8 * 0.1551
    Mplot = np.load(karlsson.get_Mplotkkp_filename(envname,sfrname))
    keepii = Mplot <= Mmax
    badii  = Mplot > Mmax
    Mplot = np.concatenate((Mplot[keepii],[Mmax]))




    kmax,kpmax,ckkp = util.load_ckkp(envname,sfrname)
    #kmax,kpmax,nbins = fMkkp.shape
    fMarr = []
    ckarr = []
    for ikp in range(1,kpmax):
        ik = ikp-1
        fMarr.append(np.concatenate((fMkkp[ik,ikp,keepii],(np.sum(fMkkp[ik,ikp,badii]),))))
        ckarr.append(ckkp[ik,ikp])
    fMarr = np.array(fMarr)
    ckarr = np.array(ckarr)
    ckarr = ckarr/np.sum(ckarr)
    #print ckarr

    plt.plot(np.arange(len(ckarr))+1,ckarr)
    plt.xlabel('k'); plt.ylabel('ck')
    plt.savefig('PLOTS/simpleck.png',bbox_inches='tight')

    ## compute chemarr
    y1 = yields.ryI05N06(0.0)
    y2 = yields.ryI05N06(1.0)

    XH = 0.75
    masstonum = 1.0/(y1.elemmass * XH)
    
    
    kplot = 5
    N = 10**5
    charr = np.zeros((kplot+1,kplot+1,N))
    feharr = np.zeros((kplot+1,kplot+1,N))
    fehbins = np.arange(-9,-2,.1)
    

    for k1 in range(kplot+1):
        for k2 in range(kplot+1):
            ktot = k1+k2; iktot=ktot-1
            if ktot == 0: continue
            a1 = draw_reshape_yields(y1,k1,N,fMarr[iktot],Mplot)
            a2 = draw_reshape_yields(y2,k2,N,fMarr[iktot],Mplot)
            a = np.concatenate((a1,a2),axis=1)
            a = np.sum(a,axis=1)*masstonum
            a = karlsson.convert_to_solar(y1.elemnames,a,verbose=False)
            charr[k1,k2,:] = a[:,0]
            feharr[k1,k2,:] = a[:,5]

    #cfearr = charr-feharr
            
    fig,axarr = plt.subplots(kplot+1,kplot+1,
                             figsize=(9,9),sharex=True,sharey=True)
    fig.subplots_adjust(hspace=0,wspace=0)
    for k1 in range(kplot+1):
        for k2 in range(kplot+1):
            ax = axarr[k1,k2]
            ktot = k1+k2; iktot=ktot-1
            if k1==0: ax.set_title(r'$k_{loFe}='+str(k2)+'$')
            if k2==0: ax.set_ylabel(r'$k_{CCSN}='+str(k1)+'$')
            if ktot == 0: ax.axis('off'); ax.yaxis.set_visible(True); continue

            feh = feharr[k1,k2,:]
            ch =  charr[k1,k2,:]
            fehplot,cfrac,err = karlsson.compute_cfrac(ch,feh,bins=fehbins)

            ax.hist(feh,normed=True,bins=fehbins,color='#D3D3D3',edgecolor='#D3D3D3')
            ax.plot(fehplot,cfrac,color='red',drawstyle='steps-mid')
            ax.plot([-4,-3],[0.75,0.3],'b-')
            ax.set_xlim((np.min(fehbins),np.max(fehbins)))
            ax.set_ylim((0,1.1))

    plt.savefig('PLOTS/simplegrid.png',bbox_inches='tight')
