
from matplotlib import pyplot as plt
import numpy as np                    
import scipy.io as sio
from gnlse import *

def ramanshift_long():
    beta2 = -10e-27
    betas = [0,0,beta2]
    gamma = 2.100000e-02
    T0 = 0.10e-12
    Nsol = 1.0
    P = Nsol**2  * np.abs(beta2) / T0**2 / gamma
    LD = T0**2/np.abs(beta2)
    Zsol = np.pi/2.0 * LD
    flength = 50.0
    simparams = prepare_sim_params(0.0, 
                               betas ,
                               1064.0e-9,
                               gamma,
                               flength,
                               11,  # Npoints
                               1.0, #tempspread
                               zpoints=100,
                               integratortype = 'dopri5',#'dop853'
                               reltol=1e-3, 
                               abstol=1e-6 ,
                               shock=False,
                               raman = True)
    inifield = np.sqrt(P) * 1.0 / np.cosh( (simparams['tvec']+1e-12)/T0)
    tf,ff,zv = perform_simulation( simparams, inifield)
    saveoutput('raman.demo', tf, ff, zv, simparams)
    #
    # output plot
    #
    d = loadoutput('raman.demo')
    inoutplot(d,zparams={"fignr":1})
    plt.show()



def self_steepening():
    """ self-steepening term: 

    compare to the Chapter 'Higher-Order Nonlinear Effects / Self-Steepening'
    
    in

    Govind P. Agrawal
    Nonlinear Fiber Optics
    Fourth Edition
    ELSEVIER
    """
    beta2 = 0.0 #no dispersion
    betas = [0,0,beta2]
    gamma = 2.100000e-02
    P = 1.0
    flength = 20 * 1./(P * gamma)
    simparams = prepare_sim_params(0.0, 
                               betas ,
                               1064.0e-9,
                               gamma,
                               flength,
                               12,  # Npoints
                               1.0, #tempspread
                               zpoints=200,
                               integratortype = 'dopri5',#'dop853'
                               reltol=1e-6, 
                               abstol=1e-12 ,
                               shock=True,
                               raman = False)
    s = 0.01 #self-steepening factor
    T0 = 1/s/simparams['om0']
    inifield = np.sqrt(P) * np.exp( -0.5* (simparams['tvec']/1.0/T0)**2)  #gaussian field

    tf,ff,zv = perform_simulation( simparams, inifield)
    saveoutput('shock.demo', tf, ff, zv, simparams)
    #
    # output plot
    #
    d = loadoutput('shock.demo')
    inoutplot(d,zparams={"fignr":2})
    plt.show()
    


def higher_order_soliton(Nsol):
    """ a higher order soliton 

    compare to the Chapter 'Fiber Solitons'
    
    in

    Govind P. Agrawal
    Nonlinear Fiber Optics
    Fourth Edition
    ELSEVIER

    """
    beta2 = -10.0e-27
    betas = [0,0,beta2]
    gamma = 2.100000e-02
    T0 = .250e-12
    P = Nsol**2 * np.abs(beta2) / gamma / T0**2

    zsol = np.pi/2.0 * T0**2/np.abs(beta2)
    flength = 3 *  zsol

    simparams = prepare_sim_params(0.0, 
                               betas ,
                               1064.0e-9,
                               gamma,
                               flength,
                               12,  # Npoints
                               1.0, #tempspread
                               zpoints=200,
                               integratortype='dop853',
                               reltol=1e-6, 
                               abstol=1e-12 ,
                               shock=False,
                               raman = False)
    inifield = np.sqrt(P) * 1.0 / np.cosh( simparams['tvec']/T0) #higher order soliton

    tf,ff,zv = perform_simulation( simparams, inifield)
    saveoutput('hos.demo', tf, ff, zv, simparams)
    #
    # output plot
    #
    d = loadoutput('hos.demo')
    inoutplot(d,zparams={"fignr":3})
    plt.show()

def soliton_self_frequency_shift_cancellation():
    """
    soliton self frequency shift cancellation: a soliton gets red-shifted due to the 
    raman effect towards the zero-dispersion wavelength of a fiber. Then the shift is
    cancelled and the soliton strongly couples energy to a dispersive wave in the region of 
    normal dispersion.
    
    See D. V. Skryabin, F. Luan, J. C. Knight, and P. S. J. Russell,
    'Soliton Self-Frequency Shift Cancellation in Photonic Crystal Fibers,'
    Science, vol 301, no 5640, pp. 1705-1708, 2003

    for comparison. I didn't get all the fiber parameters (nonlinearity), so the lengths
    for cancellation differ.

    """
    betas = [0, 0, -7.1679349121808307e-26, -3.6107349971451975e-40, 2.7392471626385486e-54, -1.10309281585831e-68, 2.3170537110678921e-83, 1.1816697954223225e-98, -2.0052740903697599e-112, 3.7468050911407223e-127]
    gamma = 0.3
   
    flength = 2.0
    simparams = prepare_sim_params(0.0, 
                                   betas ,
                                   2.99792458e8/250e12,
                                   gamma,
                                   flength,
                                   12,  # Npoints
                                   1.0, #tempspread
                                   zpoints=200,
                                   integratortype='dop853',
                                   reltol=1e-4, 
                                   abstol=1e-7 ,
                                   shock=True,
                                   raman = True,
                                  # ramantype = 'blowwood',#'hollenbeck',  #or  'blowwood', 'linagrawal'
                                  # fr=0.18
                                   )
    tau = 53e-15
    t0 = tau / 2 / np.arcsinh(1)
    p = 215.
   
    inifield = np.sqrt(p) * 1/np.cosh( (simparams['tvec']+3.2e-12)/t0)
   
    tf,ff,zv = perform_simulation( simparams, inifield)
    saveoutput('ssfs_c.demo', tf, ff, zv, simparams)
    #
    # output plot
    #
    d = loadoutput('ssfs_c.demo')
    [ax1,ax2,ax3,ax4]=inoutplot(d,zparams={"fignr":3,'clim':(-280,-220)})
    ax2.set_ylim([-280,-220])
    ax2.set_xlim([150e12,320e12])
    ax4.set_xlim([150e12,320e12])  
    plt.show()

def supercontinuumgeneration():
    """
    example of Supercontinuum generation in a PCF

    see: J. M. Dudley, G. Genty, and S. Coen, 
    'Supercontinuum generation in photonic crystal fiber,'
    Rev. Mod. Phys., vol. 78, no. 4, p. 1135, Oktober 2006.
    """

    betas = [0,0,-11.830e-3*1e-24, 8.1038e-5*1e-36, -9.5205e-8*1e-48, 2.0737e-10*1e-60,
         -5.3943e-13*1e-72, 1.3486e-15*1e-84, -2.5495e-18*1e-96, 3.0524e-21*1e-108,
         -1.7140e-24*1e-120];
    gamma = 0.1
    flength = 0.15
    simparams = prepare_sim_params(0.0, 
                               betas ,
                               835e-9,
                               gamma,
                               flength,
                               13,  # Npoints
                               1.0, #tempspread
                               zpoints=200,      
                               integratortype='dop853', 
                               reltol=1e-3, 
                               abstol=1e-6 ,
                               shock=True,
                               raman = True,
                               ramantype = 'blowwood',#'hollenbeck',  #or  'blowwood', 'linagrawal'
                                fr=0.18  )
    t0 = 28.4e-15
    p = 10e3
    inifield = np.sqrt(p) * 1./np.cosh(simparams['tvec']/t0)    
    tf,ff,zv = perform_simulation( simparams, inifield)
    saveoutput('scg.demo', tf, ff, zv, simparams)
    #
    # output plot
    #
    d = loadoutput('scg.demo')
    inoutplot(d,zparams={"fignr":3, "clim":(-360,-220),'fylim':(-360,-220)})
    plt.show()
          

def compare_raman_response_functions():
    """
    compare the three available raman response functions.
    plots the temporal and the spectral response.

    see e.g.  Q. Lin and G. P. Agrawal, 
    Raman response function for silica fibers
    OL 31, no 21 pp 3086-3088, nov 2006
    """
    simparams = prepare_sim_params(0.0, 
                               [0 ,0,0,0],
                               800e-9,
                               0.0,
                               1.0,
                               11,  # Npoints
                               1.0, #tempspread 
                                   )
    r_bw = raman_blowwood(simparams['tvec'])
    r_la = raman_linagrawal(simparams['tvec'])
    r_hc = raman_hollenbeck(simparams['tvec'])

    ft_bw = np.fft.fftshift(np.fft.ifft(np.fft.fftshift( r_bw)))
    ft_la = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(r_la)))
    ft_hc = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(r_hc)))

    plt.figure(1)
    plt.title("temporal response function R(t)")
    plt.plot( simparams['tvec']/1e-12, r_bw, color="#ff0000")
    plt.plot( simparams['tvec']/1e-12, r_la, color="#000000")
    plt.plot( simparams['tvec']/1e-12, r_hc, color="#0000ff")
    plt.xlabel("time / ps")
    plt.ylabel("response / a.u.")
    plt.legend(["Blow, Wood","Lin, Agrawal","Hollenbeck, Cantrell"])
    plt.xlim([-.2, np.max(simparams['tvec']/1e-12)])

    plt.figure(2)
    plt.subplot(211)
    plt.title("fourier transform R(omega) IMAG part")
    plt.plot( simparams['relomvec']/2e12/np.pi, np.imag(ft_bw), color="#ff0000")
    plt.plot( simparams['relomvec']/2e12/np.pi, np.imag(ft_la), color="#000000")
    plt.plot( simparams['relomvec']/2e12/np.pi, np.imag(ft_hc), color="#0000ff")
    plt.axhline(y=0, color="#999999")
    plt.legend(["Blow, Wood","Lin, Agrawal","Hollenbeck, Cantrell"])
    plt.xlim([-5,45])
    plt.xlabel("frequency / THz")
    plt.ylabel("imag part of R(omega)")

    plt.subplot(212)
    plt.title("fourier transform R(omega) REAL part")
    plt.plot( simparams['relomvec']/2e12/np.pi, np.real(ft_bw), color="#ff0000")
    plt.plot( simparams['relomvec']/2e12/np.pi, np.real(ft_la), color="#000000")
    plt.plot( simparams['relomvec']/2e12/np.pi, np.real(ft_hc), color="#0000ff")
    plt.legend(["Blow, Wood","Lin, Agrawal","Hollenbeck, Cantrell"])
    plt.axhline(y=0, color="#999999")
    plt.xlabel("frequency / THz")
    plt.ylabel("real part of R(omega)")
    plt.xlim([-5,45])
    
    plt.show()


# -------------------------------------------------------
# choose between different demos
# -------------------------------------------------------


if __name__=='__main__':
    #ramanshift_long()
    #self_steepening()
    #higher_order_soliton(4.0)
    #soliton_self_frequency_shift_cancellation()
    #supercontinuumgeneration()
    compare_raman_response_functions()
    
