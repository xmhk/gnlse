gnlse
=====
Rev 10, 22.03.2014

 * a python script to simulate the propagation of pulses in optical fibers
 * the generalized Nonlinear Schroedinger Equation (gNLSE) is modeled 

![Alt text](gnlse.png "gnlse")

 * integration is done via SCIPYs ode solvers (adaptive stepsize) 
 * it is derived from the RK4IP matlab script written by J.C.Travers, H. Frosz and J.M. Dudley that is provided in  **"Supercontinuum Generation in Optical  Fibers"**,  edited by J. M. Dudley and J. R. Taylor (Cambridge 2010).
 * see [http://scgbook.info/](http://scgbook.info/) for the original script.   

![Alt text](scg.png "supercontinuum generation example")


## available Demos (see **demos.py**):

* Raman shift (soliton self frequency shift)
* self-steepening
* higher-order soliton
* soliton self-frequency shift cancellation
* supercontinuum generation
* compare the different Raman response Models available (not actually a simulations but a temporal and spectral plot)
* losses (simple, independent of frequency) 
* frequency dependent losses (N=1 soliton center frequency changes)

## required python packages:

you will need recent **(!)** versions of 

* numpy
* scipy
* matplotlib (which is **not necessarily needed for simulations**, but required for function 'inoutplot')
* optictools (like above, only required by 'inoutplot');  it can be found on [(github/xmhk)](https://github.com/xmhk/optictools)

have a look at [http://www.scipy.org/](http://www.scipy.org/) and grab the lastest stable versions


## basic usage:

* please have a look at: **demos.py**

1. use **prepare\_sim\_params()** to prepare a time and frequency grid; declare things like the fiber's dispersion, losses, the nonlinear effects that have to be considered,  etc 
2. calculate your input field (in the time domain) 
2. use **perform\_simulation** to propagate your field
2. **saveoutput** can save the simulated data in a matlab-style file
3. **loadoutput** does the reverse thing
3. **inoutplot** can give you a quick overview (temporal and spectral in/output, false-color plot of temporal and spectral evolution)


# function reference (still incomplete)

## Functions you need to access from outside 

### prepare\_sim\_params(args, optargs)

* arguments (give in this order):
  * alpha - loss coefficient (1/m) or vector
  * betas - list or vector [beta0, beta1, beta2]
  * center wavlength (m)
  * gamma (1/(W*m))
  * fiber length (m)
  * N - number of time/frequency discretization poins (factors of 2 recommended)
  * tempspread (numerical factor, how to spread the time window, 1.0 recommended)
	    
* optional arguments
  * zpoints : number of output z-steps
  * raman True/False : include the raman effect? Standard is False
  * ramantype: choose different raman response functions:
      * 'blowwood'   Blow and D. Wood, IEEE J. of Quant. Elec., vol. 25, no. 12, pp. 2665–2673, Dec. 1989.
      * 'linagrawal' Lin and Agrawal, Opt. Lett., vol. 31, no. 21,  pp. 3086–3088, Nov. 2006.
      * 'hollenbeck' (Standard)  Hollenbeck and Cantrell J. Opt. Soc. Am. B / Vol. 19, No. 12 / December 2002
  * fr - fraction of non-instantaneous raman response. From the literature we find:
      * blowwood usually fr = 0.18 
      * linagrawal fr = 0.245
      * hollenbeck models the experimental data quite accureate. fr should be around 0.2
  * shock term True/False - include shock term? Default False
  * nsteps - maximum number of steps for one ode-integration (standard 500)
  * reltol - ode-solver relative accuracy parameter (standard 1e-6)
  * abstol - ode-solver absolute accuracy parameter (standard 1e-9)
  * integratortype: choose one of
      * 'dop853'
      * 'dopri5'
      * 'lsoda'
      * 'vode'
    
      see [corresponding section in the scipy docs](http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html#scipy.integrate.ode). Standard is dopri5

### perform\_simulation

### saveoutput	

### loadoutput

### inoutplot

## Internal functions (you usually don't have to call by yourself)