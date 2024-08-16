import numpy as np
import matplotlib.pyplot as plt

from numpy import sqrt, exp, log10, pi
from scipy.integrate import trapezoid


# Constants
c_ums = 3e14                  # c in um/s
c = 3e8                       # m/s
h = 6.626e-34                 # Planck constant   [J/Hz]
k = 1.38e-23                  # Boltzman constant [J/K]
rad2arcsec = (180/np.pi*3600) # 206265 arcsec
arcsec2rad = 1/rad2arcsec


def tophat_trans(x, center=0, fwhm=1, smoothness=0.2):

    from scipy.special import erf, erfc

    t_left  = erfc(+((2*(x-center)/fwhm)-1)/smoothness)/2
    t_right = erfc(-((2*(x-center)/fwhm)+1)/smoothness)/2

    return (t_left*t_right)


def nuInu_ZL(lambda_um, f_ZL=1.7):
    # very rough approximation for ZL
    # nuInu(sky): fit for zodiacal light [nW/m2/sr]
    # f_ZL = a fudge factor for margin
    A_scat = 3800.
    T_scat = 5500.
    b_scat = 0.4
    A_therm = 5000.
    T_therm = 270.
    nuInu = f_ZL * ( A_scat*(lambda_um**b_scat)*((lambda_um)**(-4))/(exp(h*c_ums/(k*T_scat *lambda_um))-1)
                    +A_therm*1000000           *((lambda_um)**(-4))/(exp(h*c_ums/(k*T_therm*lambda_um))-1) )
    return nuInu

def plot_SPHEREx_limit(lambda_um, fnu_limit, ax=None, label=None, **kwarg):
    lambda_min = np.array([0.75, 1.11, 1.64, 2.42, 3.82, 4.42])  # hard-coded
    lambda_max = np.array([1.11, 1.64, 2.42, 3.82, 4.42, 5.00])  # hard-coded
    color = ['darkblue', 'blue', 'green', 'gold', 'red', 'darkred']

    if ax is None:
        fig, ax = plt.subplots()

    for lmin, lmax, c in zip(lambda_min, lambda_max, color):
        gd, = np.where(np.logical_and(lambda_um >= lmin,
                                      lambda_um < lmax))
        ax.plot(lambda_um[gd], fnu_limit[gd], color=c, label=label, **kwarg)

# To remove duplicate entries in the EL COSMOS SED table. 
def remove_duplicate_in_spec(wl, ff):
    uq, uq_ind = np.unique(wl, return_index=True, return_inverse=False)
    wl = wl[uq_ind]
    ff = ff[uq_ind]


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
