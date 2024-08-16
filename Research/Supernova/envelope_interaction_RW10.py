#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 18:56:23 2022

@author: hhchoi1022
"""
def planck(wave, temp):
    import numpy as np
    #if len(wave) > 1 :
    #    print('Syntax - bbflux = planck( wave, temp)')
    #    return 0
    #if len(temp) != 1 :
    #    input('Enter a blackbody temperature : ')
    # Gives the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/a
    w = wave/1.e8 # angstroms to cm
    # constants appropriate to cgs units.
    c1 = np.float128(3.7417749e-5)          # =2*!DPI*h*c*c
    c2 = np.float128(1.4387687)             # =h*c/k
    val = c2/w/np.float128(temp)
    bbflux = c1/( (w**5)*(np.exp(val)-1.))
    return bbflux*1.e-8 # convert to ergs cm-2 s-1 A-1

def fearly2_rw10(td, rstar, band) :
    """
    This function produces the shock-heated emssion light curve. Written for the SN 2015F work based on the Rabinak & Waxman (2011) model [2015. M. Im]. Currently, the explosion energy is set to 10^51 erg which is appropriate for Type Ia SN.
    This value may need to be changed for other types of SNe.
    Also, at the end of the program, you will encounter the line
    fearly2_kasen=interpol(mbbflux,xw,6580.)
    Note that the number "6580." is the effective wavelength of the filter in Angstrom.
    In this case, it is set to R-band. Please change the number if you want to plot the light curve in different bands. [2018-05-03, added by M. Im]
    Slightly modified at 2018-05-03 to add comments. [2018-05-03, M. Im].
    """
    import numpy as np
    rstar = np.float128(rstar)
    td = np.float64(td)
    # Rsun = 6.955 * 10**10 cm
    r10 = np.float128(rstar*6.955) # rstar in Rsun unit
    r13 = np.float128(rstar*6.955e-3)
    # rstar is the radius of progenitor or companion in solar radius unit.
    # In K10 model, the emission is from an interaction btw the companion and the ejected materials so r13 is that of companion.
    # In RW11 model, the emission is from the progenitor itself. r13 is that of progenitor
    # Progenitor radius in R/10^10 cm
    Mc = np.float128(1.0/1.40)
    # Mass in chandrasekhar mass unit
    # eff  = 1.0
	# efficiency of conversion of mass to light
    Msun = 1.988e33 # g
    c    = 2.9979e10 # cm/s
    # mc2  = np.log10(Msun) + 2.*np.log10(c)
    # Energy in Msun*c**2
    eff = 1.0
    dm = np.log10(Msun) + 2. * np.log10(2.9979e10) - 51. # Difference btw E51 and Msun*c^2
    loge51 = np.log10(eff*Mc) + np.log10(Msun) + 2.*np.log10(c) - 51. - dm # explosion energy in log unit (10**51 erg)
    #ke  = np.float128(0.2) # Opacity in K10, ke = 0.2 cm**2/g -> k02 = 1, appropriate for e- scattering in fully ionized A/Z=2 elements
    k02 = np.float128(1.0) # Opacity in k/0.2cm**2/g
    #fp = 0.05 # Form factor 0.031 - 0.13 (RW11)
    fp = 0.05 # form factor 0.031 - 0.13
    #v9 = np.float128(1.)   # velocity in 10**9 cm/s unit
    logLt = 40. + np.log10(1.2) + np.log10(r10) + 0.85*loge51 - 0.69*np.log10(Mc) -0.85*np.log10(k02) -0.16*np.log10(fp) - 0.35*np.log10(td) # total luminosity : early UV/optical luminosity emitted from ejecta diffusion including radioactive + cooling emission
    logTeff = np.log10(4105.) + 0.25*np.log10(r10) + 0.016*loge51 + 0.03*np.log10(Mc) - 0.27*np.log10(k02) - 0.022*np.log10(fp) -0.47*np.log10(td)  # Effective temperature of the early emission
    c = np.float(2.9979e8) # speed of light in m/s
    # wavelengh in angstrom, 1A = 10e-10 m
    sigb = np.float128(5.6704e-5) # Stephan-Boltzman constant in CGS [egs*cm**-2 * s**-1 * ]
    logFbb = np.log10(sigb) + 4.*logTeff # Blackbody flux of early emission
    d10pc  = np.float128(10. * 3.0856e18) # 10pc in cm to convert Luminosity to Flux
    fluxm    = (logLt - 2.*np.log10(d10pc) - np.log10(4.*np.pi)) - logFbb # Flux of shock heated emission = Total flux - blackbody component of early emission
    #rph24pi = logLt - logFbb
    abzero = np.float128(-48.600) # in erg*s**-1*cm**-2*Hz**-1
    teff = 10.**(logTeff)
    xw   = 1000. + 80.*np.arange(100)
    xw = np.array(xw, dtype='float128')
    bbflux = planck(xw, teff)
    ff   = -2.5 * (2.*np.log10(xw) -10. -np.log10(c) ) # Angstrom to Hertz conversion factor
    mbbflux = -2.5*np.log10(bbflux) + ff + abzero -2.5*fluxm # convert flux to AB magnitude (Find ref in wiki).
    x_data  = xw
    y_data  = mbbflux
    from scipy.interpolate import UnivariateSpline
    spl     = UnivariateSpline(x_data, y_data, s=0.2, k=5)
    # From 1000A to 8920A, we fit BB spectrum to obtain AB magnitude in specific band from continuous value. ex) What is magnitude in 6580A? 6580 is not generated in the array wx so we fit mbbflux!
    # The last input parameter in the above function (e.g., 4770., and 6580.) are the effective wavelength of the filter in Angstrom. Please modify the value, depending on which filter you use.
    # 4770  : B band, 6580. : R band
    x_array = np.linspace(np.min(x_data), np.max(x_data), num= int((np.max(x_data) -np.min(x_data))/80.))
    if band == 'U':
        fearly2 = np.float128(spl(3656.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'B' :
        #print('Band : '+band+', eff w = 4770.')
        #fearly2 = np.float128(spl(4770)) # LSGT
        fearly2 = np.float128(spl(4353.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'V' :
        fearly2 = np.float128(spl(5477.))
    elif band == 'R' :
        #print('Band : '+band+', eff w = 6580.')
        #fearly2 = np.float128(spl(6580)) # LSGT
        fearly2 = np.float128(spl(6349.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'I' :
        #fearly2 = np.float128(spl(8175.6))
        fearly2 = np.float128(spl(8797.))
    elif band == 'g' :
        fearly2 = np.float128(spl(4770.))
    elif band == 'r' :
        fearly2 = np.float128(spl(6231.))
    return fearly2