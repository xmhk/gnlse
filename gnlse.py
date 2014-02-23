# -*- coding: utf-8 -*-

import numpy as np
from functools import partial as funcpartial
from scipy.misc import factorial
from scipy.integrate import complex_ode
from time import time



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
                        reltol = 1e-9,
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
    Retval['raman']=raman
    Retval['fr']=fr
    Retval['RW'] = RW
    Retval['shock']=shock
    Retval['gamma']=gamma
    Retval['linop']=linop
    Retval['length']=length
    Retval['W'] = W
    Retval['dz']=dz
    Retval['zpoints']=zpoints    
    Retval['reltol']=reltol
    Retval['abstol']=abstol
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
    integrator.set_integrator(simp['integratortype'], atol=simp['abstol'],rtol=simp['reltol'])
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
    for i in range(simp['zpoints']):
        instatus( integr.t, slength, t1)
        integr.integrate(integr.t + simp['dz'])
        zvec.append(integr.t)        
        tf = np.multiply ( integr.y , np.exp(simp['linop'] * (integr.t) ))
        ferg.append(tf)
        ferg2.append(np.fft.fftshift(tf)/simp['dt'])
    terg =np.fft.fft(ferg)
    return terg, ferg2,zvec
