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
from spherex_helper import tophat_trans
from scipy.ndimage import gaussian_filter
from scipy.integrate import trapezoid

#%% Constant
h = 6.626e-27
c = 3e10
rad2arcsec = (180/np.pi*3600)
arcsec2rad = 1/rad2arcsec

#%%
def synphot_7DT(spec_AA = None,
                  spec_flamb = None,
                  spec_fnu = None,
                  seeing = 1.5,
                  texp = 180,
                  atmfile = 'transmission_atm_45',
                  QEfile = 'IMX455_QE.csv',
                  show = True,
                  get_response = False
                  ):
    '''
    parameters 
    ----------
    1. spec_AA : list or np.array 
            wavelength array in Angstrom of a spectrum
    2. spec_flamb : list or np.array 
            flux array in f_lambda [erg/s/cm^2/Angstrom] with the same length of spec_AA
    3. spec_fnu : list or np.array 
            flux array in f_nu [erg/s/cm^2/Hz] with the same length of spec_AA
    4. seeing : float
            estimated seeing for photometry and sensitivity 
    5. texp : float
            total integrated exposure time [seconds]
    6. atmfile : str
            absolute path of the atmosphere file
    7. QEfile : str
            absolute path of the telescope QE file 
    8. show : bool
            if true, show response curve
    return
    ----------
    1. result_tbl : astropy.table
            synthetic photometry result
    '''

    if spec_AA != None:
        wl_AA_obs = spec_AA
        if spec_fnu == None:
            f_nu_obs = spec_flamb  * wl_AA_obs * wl_AA_obs*1e-8 / c 
        else:
            f_nu_obs = spec_fnu
    atm = Table.read(atmfile)
    T_qe = Table.read(QEfile, format = 'ascii', names =('wavelength', 'QE'))
    #%% Filter
    fwhm = 250 # [Angstrom] FWHM #########################CHANGABLE#############################
    lambda_7ds = np.arange(4000., 9000., 250)   #########################CHANGABLE#############################
    wave_lvf = np.linspace(0.1, 1.0, 1001)  #########################CHANGABLE#############################

    # Create filter_set definition
    filter_set = {'cwl': lambda_7ds,
                'wave': wave_lvf}
    if show:
        plt.figure(figsize=(10, 4), dpi = 300)
        plt.title('7DS filter transmission')
    for ii, wl_cen in enumerate(lambda_7ds):
        resp_lvf = tophat_trans(wave_lvf, center=wl_cen/1e4, fwhm=fwhm/1e4)
        filter_set.update({f'{ii}': resp_lvf})
        if show:
            plt.plot(wave_lvf, resp_lvf)
    if show:
        plt.text(0.37, 1.1, f'Bandwidth = {fwhm}, Nfilter = {len(lambda_7ds)}')
        plt.xlim(0.35, 0.95)
        plt.ylim(0.00, 1.20)

    #%% Telescope 
    # Mirror
    D = 50.5  #########################CHANGABLE#############################
    D_obs = 29.8  #########################CHANGABLE#############################
    EFL = 1537.3    # [mm]   #########################CHANGABLE#############################
    D_eff = np.sqrt(D**2 - D_obs**2)
    # Detector
    dQ_RN = 3.51           # 
    I_dark = 0.005        # [e/s], dark current  #########################CHANGABLE#############################
    pixel_size = 3.76    # [um], "pitch"  #########################CHANGABLE#############################
    theta_pixel = 0.517  # [arcsec], pixel scale  #########################CHANGABLE############################# 
    nxpix, nypix = 9576, 6388  # [pixels], detector format, approx. 9k x 6k  #########################CHANGABLE#############################
    npix_ps = np.pi*(seeing/theta_pixel)**2     # number of pixels of point source
    # QE for the detector(IMX455)
    T_qe['wavelength'] = T_qe['wavelength'].astype(float)
    T_qe['wavelength'].unit = u.um
    if show:
        plt.figure(figsize=(8,4))
        plt.plot(T_qe['wavelength'], T_qe['QE'], '-')
        plt.xlabel('wavelength [$\mu m$]')      
        plt.ylabel('QE')
        plt.title('IMX455')
    
    #%% Atmospheric transmission
    # For sky model, more info here: https://www.eso.org/sci/software/pipelines/skytools/skycalc
    # Model parameter
        # 4000~10000AA
        # altitude = 45degree
            # PMW = 2.5mm
    wl_AA_sky = atm['lam']*1e1
    nu_sky = 3e18 / wl_AA_sky
    I_lambda = atm['flux']
    f_lambda = I_lambda * (h * nu_sky) / (1e2**2) / (1e4)   # [erg/s/cm^2/AA/arcsec^2]
    f_nu_sky = f_lambda  * wl_AA_sky * wl_AA_sky*1e-8 / c # [erg/s/cm^2/Hz/arcsec^2]
    if show:
        plt.figure(figsize=(8,4), dpi = 300)
        plt.plot(atm['lam']/1e3, atm['trans'], alpha=0.5)
        plt.xlabel('wavelength [$\mu m$]')
        plt.ylabel('Transmission')
        plt.xlim(0.3,1)
    trans_smooth = gaussian_filter(atm['trans'], 10)
    if show:
        plt.plot(atm['lam']/1e3, trans_smooth)
    #%% Sky background
    if show:
        plt.figure(dpi = 300)
        plt.title('Sky brightness')
        plt.plot(atm['lam']*1e-3, f_nu_sky, linewidth  = 0.5, c= 'k')
        plt.xlabel('wavelength [$\mu m$]')
        plt.ylabel(r'$SB_\nu$  $[erg/s/cm^2/Hz/arcsec^2]$')

    #%% Total response curve (Telescope + Atmosphere + Filter )
    eff_mirrors = round(0.92 **2,2) #########################CHANGABLE#############################
    eff_optics = 0.95   #########################CHANGABLE#############################
    eff_filter = 0.95   #########################CHANGABLE#############################
    response = {'cwl': lambda_7ds,
                'wave': wave_lvf}
    intp_tel = np.interp(wave_lvf, T_qe['wavelength'], T_qe['QE'] * eff_mirrors * eff_optics)
    intp_atm = np.interp(wave_lvf, atm['lam']*1e-3, trans_smooth)
    intp_tot = intp_tel * intp_atm
    if show:
        plt.figure(figsize = (10,6), dpi = 300)
        plt.plot(wave_lvf, intp_tel, label = 'telescope', c='k')
        plt.plot(wave_lvf, intp_tel / eff_mirrors / eff_optics, alpha = 0.1, c='k')
        plt.plot(wave_lvf, intp_atm, label = 'atmosphere')
        plt.plot(wave_lvf, intp_tot, label = 'total')
        plt.text(0.8, 0.45, fr'Bandwidth = {fwhm}$\AA$, Nfilter = {len(lambda_7ds)}')
        plt.text(0.8,0.4, fr'eff_mirror = {eff_mirrors}')
        plt.text(0.8,0.35,fr'eff_optics = {eff_optics}')
        plt.text(0.8,0.3, fr'eff_filter   = {eff_filter}')
    for i in range(len(lambda_7ds)):
        resp_tot = intp_tot * filter_set[f'{i}'] * eff_filter
        response.update({f'{i}':resp_tot})
        if show:
            plt.plot(wave_lvf, resp_tot)
    if show:
        plt.xlim(0.3,1)
        plt.xlabel(r'wavelength[$\mu$m]')
        plt.ylabel('response')
        plt.legend(loc = 2)
        plt.show()

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

    #%% Photometry
    def synphot_7DT(wl_AA_obs, f_nu_obs, texp = texp):
        T_sens = (Table( 
                    names=('band', 'wavelength', 'mag', 'magerr', 'UL5_pts','SNR'),
                    dtype=(np.int16,float,float,float,float,float,) )
                )
        for key in T_sens.colnames:
            T_sens[key].info.format = '.4g'
        for i, cwl in enumerate(response['cwl']): 
            photon_rate = synth_phot(wl_AA_obs*1e-4, f_nu_obs, response['wave'], response[f'{i}'], return_photonrate= True)
            SB_photo = synth_phot(wl_AA_obs*1e-4, f_nu_obs, response['wave'], response[f'{i}'], return_photonrate= False)
            photon_rate_sky = synth_phot(wl_AA_sky*1e-4, f_nu_sky, response['wave'], response[f'{i}'], return_photonrate= True)
            SB_sky = synth_phot(wl_AA_sky*1e-4, f_nu_sky, response['wave'], response[f'{i}'], return_photonrate= False)
            
            # photo-current or count rate
            I_photo = photon_rate * (np.pi/4*D_eff**2) #* (theta_pixel**2)
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
            T_sens.add_row([i, cwl, mag_pts_noised, magerr, mag_pts_lim, SN]) 
        return T_sens
    if get_response:
        return response
    else:
        photresult = synphot_7DT(wl_AA_obs, f_nu_obs)
        plt.figure(dpi = 300)
        plt.plot(wl_AA_obs*1e-4, -2.5*np.log10(f_nu_obs) -48.60, linewidth = 1, c= 'k', alpha = 0.3)
        plt.scatter(photresult['wavelength']*1e-4,photresult['mag'], s= 10, marker = '*', c= 'r')
        plt.errorbar(photresult['wavelength']*1e-4,photresult['mag'], yerr = photresult['magerr'], fmt = 'none', elinewidth = 1, capsize = 3, c = 'r')
        #plt.ylim(10,18)
        plt.xlim(0.35,1)
        plt.ylim(17,25)
        plt.xscale('log')
        plt.show()
        return photresult

#%% Sample run 
if __name__ == '__main__':
    spectrum = Table.read('../data/ELCOSMOS/sed_434152.fits')
    photresult, response = synthphot_7DT(spec_AA = spectrum['wavelength'], spec_flamb = spectrum['flux'], texp = 180, show = False)

# %%
