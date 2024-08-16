#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 13:05:03 2022

@author: hhchoi1022
"""

#%% Library
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from skimage.transform import downscale_local_mean
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Gaussian1D
from astropy.convolution import Gaussian2DKernel
from scipy.signal import fftconvolve
#%% Constant
c_ums = 3e14                  # c in um/s
c = 3e8                       # m/s
h = 6.626e-34                 # Planck constant   [J/Hz]
k = 1.38e-23                  # Boltzman constant [J/K]
rad2arcsec = (180/np.pi*3600) # 206265 arcsec
arcsec2rad = 1/rad2arcsec
#%% Parameters for the telescope
# Mirror
D = 50.5  #########################CHANGABLE#############################
D_obs = 29.8  #########################CHANGABLE#############################
EFL = 1537.3    # [mm]   #########################CHANGABLE#############################
D_eff = np.sqrt(D**2 - D_obs**2)
theta_pixel = 0.517
dQ_RN = 3.51 
I_dark = 0.005 
nxpix = 9576
nypix = 6388
#%%
"""
Make simulated image for 7DT
1. high-resolution FITS image (Oversampling)
    7DT_highres_iamge(ra, dec, oversample, pixel_size, nx, ny, savefile)  !!!!! DONE
2. Load sourecs from the ELCOSMOS catalog : spectrum and position
    2.1. Load catalog 
    2.2. Load SED of each source
3. Add sources in the image
    3.1. Derive pixel coordinates in the image
    3.2. Inject sources to the image    
        3.2.1. Filter construction
        3.2.2. Synthetic photometry with the spectrum
        3.2.3. Conversion between flux and photon rate 
        3.2.4. Inject sources
        3.2.5. Convolution of the image

4. Downsample high-resolution FITS image to 7DT image
5. Add background 
6. Simulate instrumental effect + noises

"""
#%% 1. High resolution FITS image (oversampled)
def generate_highres_image(ra, dec, oversample=2, pixel_size=theta_pixel, nx=nxpix, ny=nypix, savefile=None):

    # Create a new WCS object.
    # WCS = World Coordinate System
    w = wcs.WCS(naxis=2)         # number of axes
    w.wcs.crval = [ra, dec]      # reference point (ra, dec)
    w.wcs.crpix = [nx*oversample/2, ny*oversample/2]         # reference pixel (xp, yp)
    w.wcs.cdelt = np.array([-pixel_size/oversample/3600.0,   # pixel scale
                             pixel_size/oversample/3600.0])

    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]  # WCS projection type
    w.wcs.crota = [0.0, 0.0]                # zero roll angle / rotation

    # Create FITS image 
    header = w.to_header()                         # astropy.io.fits.Header
    im = np.zeros((ny*oversample, nx*oversample))  # note (ny, nx) instead of (nx, ny)
    hdu = fits.PrimaryHDU(im, header=header)

    # Save to FITS file
    if savefile is not None:
        hdu.writeto(savefile, overwrite=True)
    
    return im, w

#%% 2. Load catalog & get sources
    # 2.1. Load catalog
def load_catalog(filename):
    T = Table.read(filename)
    return T
    # 2.2. Load SED
def get_sed(sedfile):
    # To remove duplicate entries in the EL COSMOS SED table. 
    from spherex_helper import remove_duplicate_in_spec
    path = '../data/ELCOSMOS/'
    
    T = Table.read(path+sedfile)
    
    wl = T['wavelength'] # angstrom
    f_lambda = T['flux'] # erg/s/cm2/A

    # conversion to f_nu
    f_nu = f_lambda * wl * (wl / 2.99792e18) / (1e-23 * 1e-6)  # micro Jansky
    wl = wl / 10000      # micron

    # Fix the input SED table
    remove_duplicate_in_spec(wl, f_nu)
    return wl, f_nu

#%% 3. Add sources in the image
    # 3.1. Find pixel coordinates
def get_pixel_coordinates(ra, dec, highres_wcs):
    """
    Retruns the detector (x,y) position for a (ra,dec) & WCS
    """
    sky = SkyCoord(ra*u.deg, dec*u.deg)
    xpos, ypos = highres_wcs.world_to_pixel(sky)
    return xpos, ypos
    
    # 3.2. Make source profile    
        # 3.2.1. Filter construction
from synphot_7DT import synthphot_7DT
response = synthphot_7DT(get_response = True)
def get_response_wl(response, wl_aa):
    idx = np.where(response['cwl'] == wl_aa)
    return response[f'{idx[0][0]}']

        # 3.2.2. Synthetic photometry with the spectrum 
def synth_phot(wave, flux, wave_lvf, resp_lvf, tol=1e-3, return_photonrate = False):
    """
    Quick synthetic photometry routine.

    Parameters
    ----------
    wave : `numpy.ndarray`
        wavelength of input spectrum.
    flux : `numpy.ndarray`
        flux density of input spectrum in f_nu unit
        if `return_countrate` = True, erg/s/cm2/Hz is assumed
    wave_lvf : `numpy.ndarray`
        wavelength of the response function
    resp_lvf : `numpy.ndarray`
        response function. assume that this is a QE.
    tol : float, optional
        Consider only wavelength range above this tolerence (peak * tol).
        The default is 1e-3.

    Returns
    -------
    synthethic flux density in the input unit
        if return_photonrate = True, photon rates [ph/s/cm2]

    """
    index_filt, = np.where(resp_lvf > resp_lvf.max()*tol)

    index_flux, = np.where(np.logical_and( wave > wave_lvf[index_filt].min(),
                                           wave < wave_lvf[index_filt].max() ))

    wave_resamp = np.concatenate( (wave[index_flux], wave_lvf[index_filt]) )
    wave_resamp.sort()
    wave_resamp = np.unique(wave_resamp)
    flux_resamp = np.interp(wave_resamp, wave, flux)
    resp_resamp = np.interp(wave_resamp, wave_lvf, resp_lvf)

    if return_photonrate:
        h_planck = 6.626e-27 # erg/Hz
        return trapezoid(resp_resamp / wave_resamp * flux_resamp, wave_resamp) / h_planck

    return trapezoid(resp_resamp / wave_resamp * flux_resamp, wave_resamp) \
         / trapezoid(resp_resamp / wave_resamp, wave_resamp)


    # 3.2.5. Convolution of the image
def get_psf(fwhm = 1.5, oversample=2):
    pixscale_oversampled = theta_pixel/oversample
    fwhm_pixel = fwhm / pixscale_oversampled
    psf = Gaussian2DKernel(fwhm_pixel / 2.35)
    return psf

#%% 4. Downsampling 
def downsample(highres_image_conv, highres_wcs, oversample=2):
    # image
    lowres_image = downscale_local_mean(highres_image_conv.data, (oversample, oversample))
    lowres_image *= oversample**2

    # wcs
    lowres_wcs = highres_wcs.copy()
    lowres_wcs.wcs.crpix = [x / oversample for x in highres_wcs.wcs.crpix]
    lowres_wcs.wcs.cdelt = [x * oversample for x in highres_wcs.wcs.cdelt]
    
    return lowres_image, lowres_wcs    

#%% Add noise
def get_noise_map(cnt, texp, dQ_RN, I_dark):
    # Read noise
    sigma_read = dQ_RN
    
    # Dark_current - ignore in this example
    sigma_dark = I_dark * texp

    # Use Gaussian instead of full Poisson noise
    varimg = cnt
    noise_map = ( np.random.normal(size=cnt.shape) * np.sqrt(varimg) +
                  np.random.normal(size=cnt.shape) * sigma_read +
                  np.random.normal(size=cnt.shape) * np.sqrt(sigma_dark)
                  )
    return noise_map

#%%
"""
Make simulated image for 7DT
1. high-resolution FITS image (Oversampling)
    7DT_highres_iamge(ra, dec, oversample, pixel_size, nx, ny, savefile)
2. Load sourecs from the ELCOSMOS catalog : spectrum and position
    2.1. Load catalog 
    2.2. Load SED of each source
3. Add sources in the image
    3.1. Derive pixel coordinates in the image
    3.2. Inject sources to the image    
        3.2.1. Filter construction
        3.2.2. Synthetic photometry with the spectrum
        3.2.3. Conversion between flux and photon rate 
        3.2.4. Inject sources
        3.2.5. Convolution of the image

4. Downsample high-resolution FITS image to 7DT image
5. Add background 
6. Simulate instrumental effect + noises

"""
def main(ra_center, 
         dec_center, 
         t_exp = 180,
         wl_cen = 8750,
         oversample = 2, 
         seeing = 1.5,
         catfile = '../data/ELCOSMOS_v1.part.fits', 
         atmfile = '/home/hhchoi1022/Desktop/Gitrepo/7DT/data/transmission_atm_45',
         savefile = '7DS_sim_image.fits'
         ):
    # Oversampling 
    highres_image, highres_wcs = generate_highres_image(ra_center, dec_center, oversample=oversample, pixel_size=theta_pixel, nx=nxpix, ny=nypix, savefile=None)
    # Getting catalog 
    T = load_catalog(catfile)
    # Filter construction
    #response = get_response()
    wave_filt = response['wave']
    resp_filt = get_response_wl(response, wl_cen)
    
    # Atmosphere
    atm = Table.read(atmfile)
    wl_AA_sky = atm['lam']*1e1
    nu_sky = 3e18 / wl_AA_sky
    I_lambda = atm['flux']
    f_lambda = I_lambda * (h * nu_sky) / (1e2**2) / (1e4) * 250   # [erg/s/cm^2/arcsec^2] 250 for equivalent width of the filter
    f_nu_sky = f_lambda  * wl_AA_sky * wl_AA_sky*1e-8 / c # [erg/s/cm^2/arcsec^2]
    photon_rate_sky = synth_phot(wl_AA_sky*1e-4, f_nu_sky, wave_filt, resp_filt, return_photonrate= True)
    I_photo_sky = photon_rate_sky * (np.pi/4*D_eff**2) * (theta_pixel)**2 * (u.electron / u.s)
    
    # Injecting sources
    for target in T:
        # In the case of the sources positioning out of the field, pass
        try: 
            ra, dec = target['ra'], target['dec']
            sedfile = target['sedfile']
            # Getting sed of the source
            wave, flux = get_sed(sedfile)
            # Finding source pixel coordinates
            xpos, ypos = get_pixel_coordinates(ra, dec, highres_wcs)
            
            # Calculation of the photon rate [e/s] of the source
            flux *= 1e-6*1e-23 * 250 #for equivalent width of the filter
            photon_rate_obs = synth_phot(wave, flux, wave_filt, resp_filt, tol=1e-3, return_photonrate=True)
            I_photo_obs = photon_rate_obs * (np.pi/4*D_eff**2) * (u.electron / u.s)
            
            # Inject the source to the image
            xpos_int, ypos_int = map(int, (xpos, ypos))
            highres_image <<= u.electron / u.s
            highres_image[ypos_int, xpos_int] = I_photo_obs
        except:
            pass
    
    
    # Convolution 
    kernel = get_psf(seeing, oversample = 2)
    highres_image_conv = fftconvolve(highres_image, kernel, mode='same')
    
    # Downsampling 
    final_image, final_wcs = downsample(highres_image_conv, highres_wcs, oversample = oversample)
    
    # Add sky background, the result final image [e/s]
    final_image += I_photo_sky.value
    
    # Apply actual image with the exposure time
    cnt = final_image * t_exp
    
    # Derive noise map (Readout noise, )
    noise_map = get_noise_map(cnt, t_exp, dQ_RN, I_dark)
    
    # Add photon
    cnt += noise_map
    
    # Save
    # ELCOSMOS_v1.part.reg
    hdu = fits.PrimaryHDU(cnt, header=final_wcs.to_header())
    hdu.writeto(f'{savefile}', overwrite=True)
    

#%%
#if __name__ == 'main':
#    main()
texp_3yr = 180 * 0.5 * 365 / 7 * 3
ra_center, dec_center = 150.06890, 2.20573
for cwl in response['cwl']:
    main(ra_center, dec_center, t_exp = 180, wl_cen = cwl, oversample = 2, savefile = f'7DT_sim_{int(cwl)}.fits')



# %%


#%%