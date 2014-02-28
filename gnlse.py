# -*- coding: utf-8 -*-

import numpy as np
from functools import partial as funcpartial
from scipy.misc import factorial
from scipy.integrate import complex_ode
from time import time
import scipy.io as sio


def beta0_curve(omvec, om0, betas):
    """
    calculate the dispersion curve via Taylor coefficients
    """
    bc = np.zeros(len(omvec))
    for i in range(len(betas)):
        bc = bc + betas[i]/factorial(i) * (omvec-om0)**i
    return bc


def ramanrespf(tvec):
    """
    temporal raman response function 
    """  
    tau1 = 12.2e-15
    tau2 = 32.0e-15
    rt = (tau1**2 + tau2**2)/(tau1 * tau2**2) *np.multiply( np.exp(-tvec/tau2), np.sin(tvec/tau1))
    rt = np.multiply( rt, tvec>=0)
    return rt



def prepare_sim_params( alpha,
                        betas ,
                        centerwavelength,                        
                        gamma, 
                        length,
                        N, 
                        tempspread=1.0,
                        fr=0.18,                        
                        raman = False,
                        shock = False,
                        nsteps = 500,
                        reltol = 1e-6,
                        abstol = 1e-9,
                        integratortype = 'vode',
                        zpoints = 256):
    """
    prepare_sim_params
    
    creates a dict containing all information necessary for 
    simulation.
    """
    # -------------------------------------------------------
    # COMPUTATION GRID
    # -------------------------------------------------------
    dt0 = centerwavelength / ( 2.0 * 2.99792458e8) #nyquist
    om0 = 2.0 * np.pi * 2.99792458e8 / centerwavelength    
    dt = tempspread * dt0
    points = 2**N
    tvec = np.arange( -points/2, points/2) * dt
    relomvec = 2 * np.pi * np.arange (-points/2, points/2)/(points * dt)
    omvec = relomvec + om0
    dz = length/zpoints    
    # -------------------------------------------------------
    # LINEAR OPERATOR: dispersion and losses 
    # -------------------------------------------------------
    if len(betas) == points: # dispersion curve as vector
        linop = 1.0j * betas
        bk = betas           # store for later use
    else:                    # dispersion via taylor coefficients
        bk = beta0_curve( relomvec, 0, betas)
        linop = 1.0j * bk
    linop += alpha/2.        # loss 
    linop = np.fft.fftshift(linop)
    # -------------------------------------------------------
    # SHOCK TERM: self-steepening
    # -------------------------------------------------------
    if shock==False:  # self-steepening off
        W = 1.0
    else:             # self-steepening on
        #gamma = gamma/om0   # <- original dudley
        #W = relomvec + om0  # <- original dudley
        #
        # i prefer the version below rather than
        # changing the gamma parameter
        # it should give the same results
        W = omvec / om0
        #
        W = np.fft.fftshift(W)
    # -------------------------------------------------------
    # Raman response function
    # -------------------------------------------------------
    RT = ramanrespf( tvec )
    RW = points*np.fft.ifft(np.fft.fftshift(RT)) 

    # -------------------------------------------------------
    # PREPARE OUTPUT DICTIONARY
    # -------------------------------------------------------
    Retval = {}    
    Retval['dt']=dt
    Retval['points']=points
    Retval['tvec']=tvec
    Retval['relomvec']=relomvec
    Retval['om0'] = om0
    Retval['omvec']=omvec
    Retval['dom'] = omvec[2]-omvec[1]
    Retval['raman']=raman
    Retval['fr']=fr
    Retval['RW'] = RW
    Retval['shock']=shock
    Retval['gamma']=gamma
    Retval['linop']=linop
    Retval['betacurve'] = bk
    Retval['length']=length
    Retval['W'] = W
    Retval['dz']=dz
    Retval['zpoints']=zpoints    
    Retval['reltol']=reltol
    Retval['abstol']=abstol
    Retval['nsteps']=nsteps
    Retval['integratortype']=integratortype
    return Retval



def prepare_integrator(simp, inifield):
    simpsub = dict( (k, simp[k]) for k in ('gamma','raman','linop','W','dz','dt','RW','fr'))
        # the line below creates a new function handle as some the scipy integrator functions
        # seem not to wrap additional parameters (simp in this case) of the RHS function 
        # correctly  (as SCIPY 0.14.0.dev-a3e9c7f)
    GNLSE_RHS2 = funcpartial( GNLSE_RHS, simp=simpsub)    
    integrator = complex_ode(GNLSE_RHS2)
        # available types  dop853   dopri5   lsoda    vode
        # zvode also available, but do not use, as complex ode handling already wrapped above
    integrator.set_integrator(simp['integratortype'], atol=simp['abstol'],rtol=simp['reltol'],nsteps=simp['nsteps'])
    integrator.set_initial_value(np.fft.ifft( inifield))
    return integrator


def GNLSE_RHS( z, AW, simp):
    """
    GNLSE_RHS
    solve the generalized Nonlinear Schroedinger equation
    this is derived from the RK4IP matlab script provided in
    "Supercontinuum Generation in Optical Fibers" edited
    by J. M. Dudley and J. R. Taylor (Cambridge 2010).
    see http://scgbook.info/ for the original script.    
    """    
    AT = np.fft.fft( np.multiply( AW , np.exp( simp['linop'] * z)))
    IT = np.abs(AT)**2  

    if simp['raman'] == True:
        RS = simp['dt']  *  np.fft.fft( np.multiply( np.fft.ifft(IT), simp['RW'] ))
        M = np.fft.ifft( np.multiply( AT, 
                                      ( (1-simp['fr'])*IT +  simp['fr'] *  RS  )
                                     )
                        )      
    else:
        M = np.fft.ifft( np.multiply( AT, IT))

    return  1.0j * simp['gamma'] * np.multiply( simp['W'], np.multiply( M, np.exp( -simp['linop'] * z)) )


def instatus( aktl, slength, t1 ):
    frac =  aktl/slength
    if frac>0.0:
        t2 = time()    
        tel = t2-t1        
        trem = (1-frac)*tel/frac
        print("%.4f m / %.4f m (%.1f%%) | %.0f s | %.0f s (%.2f h)"%(aktl,slength,frac*100, tel,trem,trem/3600.))
    
def perform_simulation( simp, inifield):  
    integr = prepare_integrator( simp, inifield)
    zvec = []
    ferg = []
    ferg2 = []
    t1 = time()
    slength = simp['length']
    zvec.append(0)
    ferg.append(np.fft.ifft( inifield))
    #
    # the fft scalingfactor ensures that the energy is conserved in both domains
    #
    scalefak = np.sqrt( simp['dt'] / simp['dom'] * simp['points'] )
    ferg2.append(np.fft.fftshift(np.fft.ifft( inifield)) *scalefak)
    for i in range(simp['zpoints']):
        instatus( integr.t, slength, t1)
        integr.integrate(integr.t + simp['dz'])
        zvec.append(integr.t)        
        tf = np.multiply ( integr.y , np.exp(simp['linop'] * (integr.t) ))
        ferg.append(tf)
        ferg2.append(np.fft.fftshift(tf) * scalefak)
    terg =np.fft.fft(ferg)
    return terg, np.array( ferg2) ,zvec



#
# input and output
#

def saveoutput(filename, tf,ff,zv, simparams):
    outputdict = {}
    outputdict['tvec'] = simparams['tvec']
    outputdict['omvec']=simparams['omvec']
    outputdict['relomvec']=simparams['relomvec']
    outputdict['om0'] = simparams['om0']
    outputdict['betacurve'] = simparams['betacurve']
    outputdict['length']=simparams['length']
    outputdict['zpoints']=simparams['zpoints']
    outputdict['points']=simparams['points']
    
    outputdict['timefield']=tf
    outputdict['freqfield']=ff
    outputdict['zvec']=zv

    outputdict['tfield1'] = tf[0,:]
    outputdict['ffield1'] = ff[0,:]
    outputdict['tfield2'] = tf[simparams['zpoints'],:]
    outputdict['ffield2'] = ff[simparams['zpoints'],:]
    sio.savemat( filename , outputdict)

def loadoutput(filename):
    d = sio.loadmat(filename)
    for k in d.keys():
        if k not in ('__header__','__globals__','__version__','timefield','freqfield'):
            #print k
            d[k]=d[k][0] #sio reconstructs cascaded arrays of (0,1)-arrays (somhow)    
    return d


def inoutplot(d,zparams={}):  # plot input and output (both domains) as well as evolution in one figure
    if 'fignr' in zparams.keys():
        plt.figure(zparams['fignr'])
    else:
        plt.figure(99)
    plt.subplot(221)
    plt.plot( d['tvec'], np.abs(d['tfield2'])**2)
    plt.plot( d['tvec'], np.abs(d['tfield1'])**2,linewidth=3)
    plt.legend(["out","in"])


    plt.subplot(222)
    plt.plot( d['omvec'], db_abs2( d['ffield2']))
    plt.plot( d['omvec'], db_abs2( d['ffield1']),linewidth=3)
    plt.legend(["out","in"])

    plt.subplot(223)
    plt.imshow( np.abs(d['timefield'])**2,aspect='auto',origin='lower')
    plt.colorbar()

    plt.subplot(224)
    ax=plt.imshow( db_abs2(d['freqfield']),aspect='auto',origin='lower')
    plt.colorbar()
    if 'clim' in zparams.keys():
        ax.set_clim(zparams['clim'])
