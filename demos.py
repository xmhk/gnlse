
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
    d = loadoutput('raman.demo')
    inoutplot(d,zparams={"fignr":1})
    plt.show()



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
    saveoutput('shock.demo', tf, ff, zv, simparams)
    d = loadoutput('shock.demo')
    inoutplot(d,zparams={"fignr":2})
    plt.show()
    




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
    saveoutput('testoutput.mat', tf, ff, zv, simparams)
    saveoutput('hos.demo', tf, ff, zv, simparams)
    d = loadoutput('hos.demo')
    inoutplot(d,zparams={"fignr":3})
    plt.show()

# -------------------------------------------------------
# choose between different demos
# -------------------------------------------------------


if __name__=='__main__':
#    ramanshift_long()
    #self_steepening()
    higher_order_soliton(4.0)
