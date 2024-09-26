#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
import numpy as np
from scipy.interpolate import interp1d

def err_weighted_combine(all_wave, all_flux, all_error, dispersion=None,
                         **kwargs):
    """
    Weigthed Combine
    - from costool procedure
      COADDITION: flux elements are coadded pixel by pixel according to several
      methods set by the keyword method. In our extensive testing, modified
      exposure time weighting with flagging seems to produce the best
      coadditions with the least spurious features. Thus we have made method=1
      the default.
        method = -1 - simple mean of pixel values, gives too much weight to
                      exposures
        method =  1 - modified exposure weighting: exposure time is modified at
                      each pixel location by flanging and wire shadows
                      (if selected).
        method =  2 - err^-2 weighting, allows error array tuning, but tends to
                      weight toward lower-flux pixels.
        method =  3 - (S/N)^2 weighting, allows for error array tuning, but
                      tends to weight toward higher-flux pixels.

    This code corresponds to method=2, the error weighted combine.
    
    This code is based on the repo below.
    https://github.com/hamogu/spectrum/blob/master/spectrum/coadd.py

    Parameters
    ----------
    all_wave : numpy array
        stacked wavelength array.
    all_flux : numpy array
        stacked flux array.
    all_error : numpy array
        stacked error array.
    dispersion : TYPE, optional
        dispersion (wavelength array) for result spectra. The default is None,
        which selects first wavelength array of given spectra.
    
    **kwargs : kwargs for scipy.interpolation.interp1d

    the number of spectra to be combined should be identical for all_wave,
    all_flux, all_error.

    Returns
    -------
    error-weighted spectrum (wavelength, flux, error)

    """
    n_spec = len(all_flux)
    
    if dispersion is None:
        dispersion = all_wave[0]
    
    spec_shape = (n_spec,len(dispersion))
    fluxes, errors = np.ma.zeros(spec_shape), np.ma.zeros(spec_shape)
    
    for i in range(n_spec):
        f_interp = interp1d(all_wave[i], all_flux[i], **kwargs)
        e_interp = interp1d(all_wave[i], all_error[i], **kwargs)
        f_new, e_new = f_interp(dispersion), e_interp(dispersion)
        fluxes[i,:] = f_new
        errors[i,:] = e_new
        
    # First, make sure there is no flux defined if there is no error.
    errors = np.ma.fix_invalid(errors)
    if np.ma.is_masked(errors):
        fluxes[errors.mask] = np.ma.masked
    # This can be simplified considerably as soon as masked quantities exist.
    fluxes = np.ma.fix_invalid(fluxes)
    # There are no masked quantities yet, so make sure they are filled here.
    flux = np.ma.average(fluxes, axis=0, weights = 1./errors**2.).filled(np.nan)
    error = np.sqrt(1. / np.ma.sum(1./errors**2., axis=0).filled(np.nan))
    
    return dispersion, flux, error

#%%
import os, glob
#filekey = '/home/hhchoi1022/Desktop/2022-11-02/pHR153-00*_w_spec.csv'
#filekey = '/home/hhchoi1022/Desktop/2022-11-04/pBD+75d325-00*_w_spec.csv'
filekey = '/home/hhchoi1022/Desktop/2022-11-02/pSN2022xkq-000*_wf_spec.csv'
#filekey = '/home/hhchoi1022/Desktop/2022-11-04/pNGC2392-000*_wf_spec.csv'
filelist = glob.glob(filekey)
wave_array = []
flux_array = []
fluxerr_array = []
# %%
from astropy.io import ascii
for file_ in filelist:
    data = ascii.read(file_)
    #wave, flux, fluxerr = data['wave'], data['inten'], data['std']
    wave, flux, fluxerr = data['wave'], data['flux'], data['error']
    wave_array.append(wave)
    flux_array.append(flux)
    fluxerr_array.append(fluxerr)
# %%
from astropy.table import Table
com_wave, com_flux, com_fluxerr = err_weighted_combine(all_wave = wave_array, all_flux = flux_array, all_error = fluxerr_array)
combined_tbl = Table()
combined_tbl['wavelength'] = com_wave
combined_tbl['flux'] = com_flux
combined_tbl['fluxerr'] = com_fluxerr
combined_tbl = combined_tbl[combined_tbl['wavelength'] > 4000]
combined_tbl = combined_tbl[combined_tbl['wavelength'] < 8000]
combined_tbl.write('/home/hhchoi1022/Desktop/2022-11-02/pSN2022xkq_w_spec.csv', format = 'ascii')
#combined_tbl.write('/home/hhchoi1022/Desktop/2022-11-04/pNGC2392_w_spec.csv', format = 'ascii')
#combined_tbl.write('/home/hhchoi1022/Desktop/2022-11-04/pBD+75d325_w_spec.csv', format = 'ascii')
#combined_tbl.write('/home/hhchoi1022/Desktop/2022-11-04/pNGC2392_w_spec.csv', format = 'ascii')
# %% Confirm
import matplotlib.pyplot as plt
for wave, flux, fluxerr in zip(wave_array, flux_array, fluxerr_array):
    plt.errorbar(wave, flux, yerr = fluxerr, c = 'k')
plt.errorbar(com_wave, com_flux, com_fluxerr, c = 'k', label ='Single')
plt.errorbar(com_wave, com_flux, com_fluxerr, c = 'r', label = 'Combined')
#plt.ylim(0, 2e-12)
plt.legend()

# %%
