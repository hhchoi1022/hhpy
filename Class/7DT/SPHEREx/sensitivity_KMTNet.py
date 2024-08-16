# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

###### Sensitivity of the telescope
###### Author : Hyeonho Choi 
###### Date : 22/06/15

#%% Library
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.table import Table
homedir = '/home/hhchoi1022/Desktop/Gitrepo/Class/7DT/SPHEREx'
os.chdir(homedir) 
from spherex_helper import tophat_trans
from scipy.ndimage import gaussian_filter
from scipy.integrate import trapezoid
check = False

#%% Constant
h = 6.626e-27
c = 3e10
rad2arcsec = (180/np.pi*3600)
arcsec2rad = 1/rad2arcsec

##### Filter info : filter_set
##### Total response info : response
##### Parameters
"""
    fwhm    <<< fwhm for aperture
    texp    <<< exposure time
    lambda_7ds    <<< filterset construction
    wave_lvf    <<< [um] lamb_start, lamb_end
    eff_mirrors    <<< efficiency of the mirror
    eff_optics    <<< efficiency of the optics
    eff_filter    <<< efficiency of the filter
    and many configuration for the telescope & detector
"""
#%% Filter
from astropy.io import ascii
B = ascii.read('../data/B.csv')
V = ascii.read('../data/V.csv')
R = ascii.read('../data/R.csv')
I = ascii.read('../data/I.csv')

plt.plot(B['wavelength'],B['QE'])
plt.plot(V['wavelength'],V['QE'])
plt.plot(R['wavelength'],R['QE'])
plt.plot(I['wavelength'],I['QE'])

#%%
filter_name = ['B','V','R','I']
filter_response = [B,V,R,I]
color = {'B':'b',
     'V':'g',
     'R':'r',
     'I':'y'}
lambda_KMTNet = [4353, 5477, 6349, 8797]   #########################CHANGABLE#############################
wave_lvf = np.linspace(0.1, 1.1, 1001)  #########################CHANGABLE#############################V
# Create filter_set definition
filter_set = {'cwl': lambda_KMTNet,
              'wave': wave_lvf}


for filt_name, filt_res, wl_cen in zip(filter_name,filter_response,lambda_KMTNet):
    resp_lvf = np.interp(wave_lvf*1e3, filt_res['wavelength'],filt_res['QE'])
    filter_set.update({f'{filt_name}': resp_lvf})


plt.figure(figsize=(10, 4), dpi = 300)
plt.title('KMTNet filter transmission')
for filt_name in filter_name:
    resp = filter_set[f'{filt_name}']
    plt.plot(wave_lvf, resp, label = f'{filt_name}', c = color[f'{filt_name}'])
plt.legend()
plt.xlim(0.35, 1.1)
plt.ylim(0.00, 1.0)
#%% Telescope 
# Mirror
D = 160  #[cm]#########################CHANGABLE#############################
D_obs = np.sqrt(200/np.pi)  #[cm]#########################CHANGABLE#############################
EFL = 5160    # [mm]   #########################CHANGABLE#############################
D_eff = np.sqrt(D**2 - D_obs**2)

# Detector
dQ_RN = 10           # [e], readout noise   #########################CHANGABLE#############################
I_dark = 0.00111        # [e/s], dark current 4e-/pixel/hour #########################CHANGABLE#############################
pixel_size = 10    # [um], "pitch"  #########################CHANGABLE#############################
theta_pixel = 0.40  # [arcsec], pixel scale  #########################CHANGABLE############################# 
nxpix, nypix = 9216, 9232  # [pixels], detector format, approx. 18k x 18k  #########################CHANGABLE#############################
# QE for the detector(IMX455)
T_qe = Table.read('/home/hhchoi1022/Desktop/Gitrepo/Class/7DT/data/e2VDD.csv', format = 'ascii', names =('wavelength', 'QE'))
T_qe['wavelength'] = T_qe['wavelength'].astype(float)*1e-3
T_qe['wavelength'].unit = u.um

plt.figure(figsize=(8,4), dpi =300)
plt.plot(T_qe['wavelength'], T_qe['QE'], '-', c = 'k')
plt.xlabel('wavelength [$\mu m$]')
plt.ylabel('QE')
plt.title('e2V DD')

#%% Atmospheric transmission
# For sky model, more info here: https://www.eso.org/sci/software/pipelines/skytools/skycalc
# Model parameter
    # 4000~10000AA
    # altitude = 45degree
    # PMW = 2.5mm
atm = Table.read('/home/hhchoi1022/Desktop/Gitrepo/Class/7DT/data/transmission_atm_45')
wl_AA_sky = atm['lam']*1e1
nu_sky = 3e18 / wl_AA_sky
I_lambda = atm['flux']
f_lambda = I_lambda * (h * nu_sky) / (1e2**2) / (1e4)   # [erg/s/cm^2/AA/arcsec^2]
f_nu_sky = f_lambda  * wl_AA_sky * wl_AA_sky*1e-8 / c # [erg/s/cm^2/Hz/arcsec^2]

plt.figure(figsize=(8,4), dpi = 300)
plt.plot(atm['lam']/1e3, atm['trans'], alpha=0.5)
plt.xlabel('wavelength [$\mu m$]')
plt.ylabel('Transmission')
plt.xlim(0.3,1)
trans_smooth = gaussian_filter(atm['trans'], 10)
plt.plot(atm['lam']/1e3, trans_smooth)
#%% Total response curve (Telescope + Atmosphere + Filter )
eff_optics = 0.95   #########################CHANGABLE#############################
response = {'cwl': lambda_KMTNet,
            'wave': wave_lvf}

intp_tel = np.interp(wave_lvf, T_qe['wavelength'], T_qe['QE']  * eff_optics)
intp_atm = np.interp(wave_lvf, atm['lam']*1e-3, trans_smooth)
intp_tot = intp_tel * intp_atm


plt.figure(figsize = (10,6), dpi = 300)
plt.plot(wave_lvf, intp_tel, label = 'telescope', c='k')
plt.plot(wave_lvf, intp_tel / eff_optics, alpha = 0.1, c='k')
plt.plot(wave_lvf, intp_atm, label = 'atmosphere')
#plt.text(0.8,0.35,f'eff_optics = {eff_optics}')
for filt_ in filter_name:
    resp_tot = intp_tot * filter_set[f'{filt_}']
    response.update({f'{filt_}':resp_tot})
    plt.plot(wave_lvf, resp_tot, label = f'{filt_}', c = color[f'{filt_}'])
    plt.plot(wave_lvf, filter_set[f'{filt_}'], c = color[f'{filt_}'], alpha = 0.2)
plt.xlim(0.3,1)
plt.xlabel(r'wavelength[$\mu$m]')
plt.ylabel('response')
plt.legend(loc = 2)
plt.show()


#%% Photometry
    #%% Synthetic photometry

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
        if return_photonrate = True, photon rates [ph/s/cm2/arcsec^2]
        else return = [erg/s/cm^2/arcsec^2]

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
        
    return trapezoid(resp_resamp / wave_resamp * flux_resamp, wave_resamp) / trapezoid(resp_resamp / wave_resamp, wave_resamp)

def make_noise_source(mag, SN):
    magerr = 2.5*np.log10(1+1/SN)
    noised_mag = mag + np.random.normal()*magerr
    return noised_mag, magerr
    
    #%% Target photometry 
    # Target 
seeing = 2
npix_ps = np.pi*(seeing/theta_pixel)**2     # number of pixels of point source
npix_1s = np.pi*(1/theta_pixel)**2     # number of pixels of 1''
texp = 240
spectrum = Table.read('../data/ELCOSMOS/sed_434152.fits')   
wl_AA_obs = spectrum['wavelength']
f_nu_obs = spectrum['flux']  * wl_AA_obs * wl_AA_obs*1e-8 / c # [erg/s/cm^2/Hz]


#%% Photometry
def synphot_KMTNet(wl_AA_obs, f_nu_obs, texp = 240):
    T_sens = (Table( 
                 names=('band', 'wavelength', 'I_photo_sky', 'mag_pxl', 'mag', 'magerr', 'UL5_pxl', 'UL5_pts','SN'),
                 dtype=(np.int16,float,float,float,float,float,float,float,float,) )
             )
    for key in T_sens.colnames:
        T_sens[key].info.format = '.4g'
    for i, cwl in enumerate(response['cwl']): 
        filt = filter_name[i]
        photon_rate = synth_phot(wl_AA_obs*1e-4, f_nu_obs, response['wave'], response[f'{filt}'], return_photonrate= True)
        SB_photo = synth_phot(wl_AA_obs*1e-4, f_nu_obs, response['wave'], response[f'{filt}'], return_photonrate= False)
        photon_rate_sky = synth_phot(wl_AA_sky*1e-4, f_nu_sky, response['wave'], response[f'{filt}'], return_photonrate= True)
        SB_sky = synth_phot(wl_AA_sky*1e-4, f_nu_sky, response['wave'], response[f'{filt}'], return_photonrate= False)
        
        # photo-current or count rate
        I_photo = photon_rate * (np.pi/4*D_eff**2)
        sky_photo =  photon_rate_sky * (np.pi/4*D_eff**2) * (theta_pixel**2)
        
        # noise in count per obs [e]. 
        Q_photo = (I_photo+I_dark)*texp
        dQ_photo = np.sqrt(Q_photo)
        Q_photo_sky = (sky_photo+I_dark)*texp
        dQ_photo_sky = np.sqrt(Q_photo_sky)
        
        # noise in count rate [e/s]
        # read-noise (indistinguishable from signal) should be added 
        dI_photo = np.sqrt(dQ_photo**2 + dQ_RN**2)/texp
        dI_photo_sky = np.sqrt(dQ_photo_sky**2 + dQ_RN**2)/texp
        SN = I_photo*texp / np.sqrt(I_photo*texp + npix_ps*sky_photo*texp + npix_ps*I_dark*texp + npix_ps*dQ_RN**2)
        
        # surface brightness [per arcsec]
        dSB_photo = (dI_photo/I_photo)*SB_photo
        dSB_photo_sky = (dI_photo_sky/sky_photo)*SB_sky
        mag_pxl = -2.5*np.log10(SB_photo*(theta_pixel)**2) - 48.60
        mag_pxl_lim = -2.5*np.log10(5*dSB_photo_sky) - 48.60
        
        # mag [per arcsec]
        #dFnu = np.sqrt(npix_ps) * dSB_photo*(theta_pixel)**2
        dFnu_sky = np.sqrt(npix_ps) * dSB_photo_sky*(theta_pixel)**2
        Fnu = SB_photo
        mag_pts = -2.5*np.log10(Fnu) - 48.60
        mag_pts_lim = -2.5*np.log10(5*dFnu_sky) - 48.60
        mag_pts_noised, magerr = make_noise_source(mag_pts, SN)

        # Add data to the table
        T_sens.add_row([i, cwl, I_photo, mag_pxl, mag_pts_noised, magerr, mag_pxl_lim, mag_pts_lim, SN]) 
    return T_sens

#%%
'''
T_sens = synphot_7DT(wl_AA_obs, f_nu_obs)
T_sens
T_sens
T_sens
2.5*np.log10(1.794/5)
plt.figure(figsize = (12,6), dpi = 300)

plt.plot(wl_AA_obs, -2.5*np.log10(f_nu_obs)-48.6, alpha =0.4)
plt.scatter(T_sens['wavelength'], T_sens['mag'], c= 'r')
plt.errorbar(T_sens['wavelength'],T_sens['mag'], yerr = T_sens['magerr'], capsize = 5)
plt.xlim(0,10000)
#plt.ylim(25,18)
#plt.plot(np.array([4250, 6750, 8250]), [20.2, 19.8, 18.9], 'o', c='orange')

T_sens
'''