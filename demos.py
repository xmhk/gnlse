
from matplotlib import pyplot as plt
import numpy as np                    
import scipy.io as sio
from gnlse import *

def helperpolot( tf,ff,simparams):
    plt.figure(1)
    plt.imshow(np.abs(tf)**2,aspect='auto')    
    plt.colorbar()
#    plt.savefig("3dtime.png")

    plt.figure(2)
    plt.imshow(10*np.log10(np.abs(ff)**2),aspect='auto',cmap='hot')#,clim=(-90,20))
    plt.colorbar()
#    plt.savefig("3dfreq.png")

    plt.figure(3)
    plt.subplot(211)
    plt.plot(simparams['tvec'],np.abs( tf[0,:])**2)

    plt.plot(simparams['tvec'],np.abs( tf[simparams['zpoints'],:]**2))
    plt.legend(["in",'out'])
    plt.subplot(212)
    plt.plot(np.angle(tf[0,:]))
    plt.plot(np.angle(tf[simparams['zpoints'],:]))
    plt.legend(['phase in','phase out'])

    plt.figure(4)
    plt.subplot(211)
    ff = np.array(ff)
    plt.plot(simparams['omvec']/2e12/np.pi, np.log( np.abs( ff[0,:])**2))
    plt.plot(simparams['omvec']/2e12/np.pi, np.log( np.abs( ff[simparams['zpoints']-1,:])**2))    
    plt.legend(["spek in ","spek out"])
    plt.subplot(212)
    plt.plot(simparams['omvec']/2e12/np.pi,( np.abs( ff[0,:])**2),"-")
    plt.plot(simparams['omvec']/2e12/np.pi, ( np.abs( ff[simparams['zpoints']-1,:])**2),"-")        
    plt.legend(["spek in ","spek out"])
    plt.show()



def ramanshift_long():
    beta2 = -5e-27
    betas = [0,0,beta2]
    gamma = 2.100000e-02
    T0 = 0.50e-12
    Nsol = 1.0
    P = Nsol**2  * np.abs(beta2) / T0**2 / gamma
    LD = T0**2/np.abs(beta2)
    Zsol = np.pi/2.0 * LD
    flength = 1500.0
    simparams = prepare_sim_params(0.0, 
                               betas ,
                               1064.0e-9,
                               gamma,
                               flength,
                               12,  # Npoints
                               1.0, #tempspread
                               zpoints=1000,
                               integratortype = 'dopri5',#'dop853'
                               reltol=1e-6, 
                               abstol=1e-12 ,
                               shock=False,
                               raman = True)
    inifield = np.sqrt(P) * 1.0 / np.cosh( simparams['tvec']/T0)
    tf,ff,zv = perform_simulation( simparams, inifield)
    helperpolot( tf,ff,simparams)


def self_steepening():
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
    # see agrawal chapter
    s = 0.01 #self-steepening factor
    T0 = 1/s/simparams['om0']
    inifield = np.sqrt(P) * np.exp( -0.5* (simparams['tvec']/1.0/T0)**2)  #gaussian field

    tf,ff,zv = perform_simulation( simparams, inifield)
    helperpolot( tf,ff,simparams)


def higher_order_soliton(Nsol):
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
                             #  integratortype = 'dopri5',#'dop853'
                               integratortype='dop853',
 #                              integratortype='lsoda',
                               reltol=1e-6, 
                               abstol=1e-12 ,
                               shock=False,
                               raman = False)
    inifield = np.sqrt(P) * 1.0 / np.cosh( simparams['tvec']/T0) #higher order soliton

    tf,ff,zv = perform_simulation( simparams, inifield)
    helperpolot( tf,ff,simparams)

# -------------------------------------------------------
# choose between different demos
# -------------------------------------------------------


if __name__=='__main__':
    #ramanshift_long()
    #self_steepening()
    higher_order_soliton(4.0)
