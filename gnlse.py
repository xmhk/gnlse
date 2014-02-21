# -*- coding: utf-8 -*-

import numpy as np
from optictools import beta0_curve
from functools import partial as funcpartial
from scipy.integrate import complex_ode
from time import time

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
                        reltol = 1e-5,
                        abstol = 1e-5,
                        integratortype = 'vode',
                        zpoints = 256):
    ###
    dt0 = centerwavelength / ( 2.0 * 2.99792458e8) #nyquist
    om0 = 2.0 * np.pi * 2.99792458e8 / centerwavelength
    dt = tempspread * dt0
    points = 2**N
    tvec = np.arange( -points/2, points/2) * dt
    relomvec = 2 * np.pi * np.arange (-points/2, points/2)/(points * dt)
    if len(betas) == points: # betas as full dispersion vector
        linop = 1.0j * betas
    else: # betas as taylor coefficients
        bk = beta0_curve( relomvec, 0, betas)
        linop = 1.0j * bk
    if shock==False:   # self steepening - FIXIT
        W = 1.0
    else:              # self steepening - FIXIT
        W = 1.0
    linop += alpha/2.        # loss 
    linop = np.fft.fftshift(linop)
  
    dz = length/zpoints    # z-steps
    ###
    Retval = {}
    Retval['dt']=dt
    Retval['points']=points
    Retval['tvec']=tvec
    Retval['relomvec']=relomvec
    Retval['raman']=raman
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


def GNLSE_RHS( z, AW, simp):
    """
    GNLSE_RHS
    solve the generalized Nonlinear Schroedinger equation
    this is derived from the RK4IP matlab script provided in
    "Supercontinuum Generation in Optical Fibers" Edited by J. M. Dudley and J. R. Taylor (Cambridge 2010).
    see http://scgbook.info/ for the original script.    
    """
    
    AT = np.fft.fft( np.multiply( AW , np.exp( simp['linop'] * z)))
    IT = np.abs(AT)**2
    if simp['raman'] == True:
        print "dummy"
    else:
        M = np.fft.ifft( np.multiply( AT, IT))

#    if len(simp.RT)==1 or np.abs(fr)==0:
    #if np.abs(simp.fr)==0:
    #    M = np.fft.ifft( np.multiply( AT, IT))

    #else:
    #    RS = simp.dt * simp.fr * np.fft.fft(   np.multiply ( np.fft.ifft(IT), simp.RW))
    #    M = np.fft.ifft( np.multiply( AT, (1-self.fr)*IT ) + RS)
    return  1.0j * simp['gamma'] * np.multiply( simp['W'], np.multiply( M, np.exp( -simp['linop'] * z)) )

def prepare_integrator(simp, inifield):
    simpsub = dict( (k, simp[k]) for k in ('gamma','raman','linop','W','dz'))
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

def instatus( aktl, slength, t1 ):
    t2 = time()    
    frac =  aktl/slength
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
        ferg2.append(np.fft.fftshift(tf))
    terg =np.fft.fft(ferg)
    return terg, ferg2,zvec
