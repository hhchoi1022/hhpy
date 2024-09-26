#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 18:49:13 2022

@author: hhchoi1022
"""
#%%
import numpy as np
from typing import Optional
import pyphot
from scipy.interpolate import UnivariateSpline
from astropy.table import Table
from HHsupport_analysis import load_filt_keys
import matplotlib.pyplot as plt
#%%
class CompanionInteraction:
    """
    ### Generate Companion interaction model for Type Ia supernova (Kasen 2010)
    ##### Written by prof.Myungshin Im &&& modified by Hyeonho Choi & Gu Lim
    =======
    Parameters
    =======
    1. rstar : float = radius of the companion star in solar radius unit
    2. m_wd : float = mass of the white dwarf in solar mass unit 
    3. commonangle : bool = commonangle or optimal angle
            If false, calculate the optimal luminosity of the model (Optimal angle)
            If true, calculate the averaged angle that can be detected 
            
    =======
    Methods
    =======
    1. calc_spectrum : calculate the expected spectrum for given phase
    2. calc_magnitude : calculate the expected magnitude for given phase & filterset
    """
    
    def __init__(self,
                 rstar : float,
                 m_wd : float,
                 commonangle : bool = False
                 ):
        
        self.rstar = rstar
        self.wdmass = m_wd
        self.commonangle = commonangle

    def __repr__(self) -> str:
        rstar = self.rstar
        wdmass = self.wdmass
        commonangle = self.commonangle
        
        return (f'CompanionInteractionModel(Companion radius:{rstar}, '
                f'WD mass:{wdmass}, '
                f'CommonAngle:{commonangle})'
        )
    
        
    def _calc_planck(self, wave : Optional[np.array], temp : Optional[np.array]):
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
    
    def _get_filt_lpivot(self, filterset : str):
        lib = pyphot.get_library()
        filter_key = dict(u = 'SDSS_u',
                        g = 'SDSS_g',
                        r = 'SDSS_r',
                        i = 'SDSS_i',
                        U = 'GROUND_JOHNSON_U',
                        B = 'GROUND_JOHNSON_B',
                        V = 'GROUND_JOHNSON_V',
                        R = 'GROUND_COUSINS_R')
        lpivotlist = []
        for band in filterset:
            if band not in 'ugriUBVR':
                print(f'{band} is not registered.')
            cen_wl = lib[filter_key[band]].lpivot.value
            lpivotlist.append(cen_wl)
        return np.array(lpivotlist)
    
    def calc_spectrum(self, td):
        # rstar = 1.0; td    = 0.5 ; band  = 'B'
        rstar = np.float128(self.rstar)
        td = np.float64(td)
        # Rsun = 6.955 * 10**10 cm
        r10 = np.float128(rstar*6.955) # rstar in Rsun unit
        r13 = np.float128(rstar*6.955e-3)

        # rstar is the radius of progenitor or companion in solar radius unit.
        # In K10 model, the emission is from an interaction btw the companion and the ejected materials so r13 is that of companion.
        # In RW11 model, the emission is from the progenitor itself. r13 is that of progenitor
        # Progenitor radius in R/10^10 cm
        Mc = np.float128(self.wdmass/1.40)
        # Mass in chandrasekhar mass unit
        # eff  = 1.0
        # efficiency of conversion of mass to light
        # Msun = 1.988e33 # g
        # c    = 2.9979e10 # cm/s
        # mc2  = np.log10(Msun) + 2.*np.log10(c)
        # Energy in Msun*c**2
        ke  = np.float128(1.0) # Opacity in K10, ke = 0.2 cm**2/g -> k02 = 1, appropriate for e- scattering in fully ionized A/Z=2 elements
        k02 = np.float128(1/5.0) # Opacity in k/0.2cm**2/g
        #fp = 0.05 # Form factor 0.031 - 0.13 (RW11)
        v9 = np.float128(1.)   # velocity in 10**9 cm/s unit
        logLt = 43. + np.log10(2*r13) + 0.25*np.log10(Mc) + (7./4.)*np.log10(v9) + (-0.75)*np.log10(ke) + (-0.5)*np.log10(td)
        if self.commonangle :
            logLt = logLt - 1.
        # Comoving isotropic collisional bolometric luminosity (L_c,iso = Ejecta + L_cooling; BB)
        logTeff = np.log10(2.5) + 4. + 0.25*np.log10(2.*r13) - (35./36.)*np.log10(ke) - (37./72.)*np.log10(td)
        # Effective temperature
        c = np.float(2.9979e8) # speed of light in m/s
        # wavelengh in angstrom, 1A = 10e-10 m
        sigb   = np.float128(5.6704e-5)      # Stephan-Boltzman constant in CGS [egs*cm**-2 * s**-1 * ]
        logFbb = np.log10(sigb) + 4.*logTeff # flux of early emission (cooling) assuming a blackbody
        d10pc  = np.float128(10. * 3.0856e18)# 10pc in cm to convert Luminosity to Flux
        fluxm  = (logLt - 2.*np.log10(d10pc) - np.log10(4.*np.pi)) - logFbb # Cooling = Lc,iso - Ejecta BB
        # What is "?"-band flux from the bolometric luminosity?
        xw   = 1000. + 40.*np.arange(200)
        xw   = np.array(xw, dtype='float128')
        teff    = 10.**(logTeff)   # Blackbody temperature for bolometric collisional luminosity
        bbflux  = self._calc_planck(xw, teff) # Black body flux of total collisional luminosity [ergs cm-2 s-1 A-1]
        ff      = -2.5 * (2.*np.log10(xw) -10. -np.log10(c) ) # Angstrom to Hertz conversion factor
        abzero  = np.float128(-48.600)                        # in erg*s**-1*cm**-2*Hz**-1
        mbbflux = abzero -2.5*(np.log10(bbflux) + fluxm) + ff # (Flux to AB mag of Lc,iso + cooling)
        #mbbflux = abzero -2.5*(np.log10(bbflux) ) + ff # (Flux to AB mag of Lc,iso + cooling)
        #if CommonAngle :
        #    mbbflux = mbbflux + 2.5 # In log scale
        x_data  = xw
        y_data  = mbbflux
        spl     = UnivariateSpline(x_data, y_data, s=0.2, k=5)
        return spl

    def calc_magnitude(self, td : Optional[np.array], filterset : str, show : bool = False):
        tbl_names = ['phase'] + list(filterset)
        mag_tbl = Table(names = tbl_names)
        for day in td:
            spline_model = self.calc_spectrum(td = day)
            lpivotlist = self._get_filt_lpivot(filterset)
            filt_mag = list(spline_model(lpivotlist))
            mag_tbl.add_row([day] + filt_mag)
        if show:
            self._visualize(mag_tbl, dpi = 200)
        return mag_tbl
    
    def _visualize(self, mag_tbl,
                   dpi : int = 300):
        color_key, offset_key, _, _, label_key = load_filt_keys()
        plt.figure(dpi = dpi)
        plt.gca().invert_yaxis()
        for filter_ in mag_tbl.keys()[1:]:
            clr = color_key[filter_]
            offset = offset_key[filter_]
            label = label_key[filter_]
            plt.plot(mag_tbl['phase'], mag_tbl[filter_] + offset, color = clr, label = label)
        plt.legend()
        plt.xlabel('Phase[days]')
        plt.ylabel('Magnitude[AB]')
        plt.show()
# %%
solar_r = 6.955e10
td = np.arange(0.0, 15, 0.1)
A = CompanionInteraction(rstar = 2e12/2/solar_r, m_wd = 1.2, commonangle= False)
A.calc_magnitude(td, filterset = 'UBV', show = True)


#%%
from HHsupport_analysis import mag_to_fnu
from HHsupport_analysis import fnu_to_mag
# %%
mag_to_fnu(-16)
# %%
flux = 10**((-16+48.6)/(-2.5))
# %%
-2.5*np.log10(1e43 /4 / np.pi / (np.float128(10. * 3.0856e18)**2))
# %%
10**((-16 + 48.6)/(-2.5)) * (4 * np.pi * (10 * 3e18)**2)
# %%
