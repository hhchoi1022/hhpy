

#%%
from astropy import constants as const
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
#%%
h = const.h.cgs
c = const.c.cgs
k = const.k_B.cgs
pi = np.pi
pc = const.pc.cgs
m_hydrogen = 1.67e-24 * u.g
m_sun = const.M_sun.cgs
r_sun = const.R_sun.cgs
m_u = const.u.cgs
# %% Prob 2
flux = 20 * u.Jy * u.km / u.s
flux_ergs = flux.to(u.erg/u.s/u.cm**2/u.Hz*u.cm/u.s)
nu0 = 1.42 * 1e9 * u.Hz
d = 15 * 1e6 * pc
A21 = 2.9e-15 / u.s
N2 = (flux_ergs*((1)/(4*pi*d**2)*(h*nu0*A21*c/nu0))**-1).value
NH = flux_ergs * (16*pi*d**2)/(3)/(A21*h*c)
#NH = 4/3*N2
M_H = m_hydrogen * NH
M_H_sun = M_H/m_sun

#%% Prob 3 
# prob31
ans31 = r_sun/c
# prob32
crosssection_thomson = 6.65e-25 * u.cm**2
h_ratio = 0.73
he_ratio = 0.25
N_h_e = h_ratio * m_sun / m_hydrogen
N_he_e = 2* he_ratio * m_sun / (4*m_hydrogen)
N_e = N_h_e + N_he_e
#N_tot = m_sun / (h_ratio * m_u + 4* he_ratio * m_u)
#N_e = h_ratio * N_tot + 2 * he_ratio * N_tot
V_sun = (4/3)*pi*r_sun**3
n_e = N_e/V_sun
ans32 = (n_e * crosssection_thomson)**-1
# prob33
distance_move = r_sun**2 / ans32
ans33_sec = (distance_move / c).value * u.s 
ans33_yr = (ans33_sec / 86400/365).value * u.yr
# prob34
#%%
N_interaction = r_sun**2 / ans32 **2
time_interaction = 1e-8 * u.s
ans34_sec = (N_interaction * time_interaction) + ans33_sec
ans34_yr = (ans34_sec/86400/365).value * u.yr

#%% Prob 4
import math

T_star = 4000 # K
R_star = 2.5 # solar radius
ww = np.append(1000+300*np.arange(330), (100000 + 5000*np.arange(30000))) * u.AA
nulist = (c/(ww*1e-8)).value*u.Hz
nulist = (1e10 + 1e10 * np.arange(1e5)) * u.Hz

def Temperature_disk(r,
                     T_star : float = 4000,
                     R_star : float = 2.5):
    """
    Calc disk temperature at distance r with the central star (T_star & R_star)

    Args:
        r : distance from the center of the star in solar radius 
        T_star (float, optional): Central star effective temperature[K]. Defaults to 4000.
        R_star (float, optional): Central star radius in solar radius. Defaults to 2.5.

    Returns:
        temp_disk : disk temperature at the distance r
    """
    term1 = T_star
    term2 = math.asin(R_star/r)
    term3 = R_star/r
    term4 = (1-(R_star/r)**2)**(0.5)
    temp_disk = term1 * (1/pi*(term2 - term3 * term4)) **(0.25)
    return temp_disk
#%%

def planck(temp,
            wave = None,
            nu = None):
    """Calculate planck function for given temperature(K) 

    Args:
        temp : Effective temperature(K)
        wave (optional): wavelength in Angstron unit
        nu (optional): Frequency in Hz unit

    Returns:
        dict: planck functiond in f_nu, f_lamb
    """
    #temp = Temperature_disk(r = r, T_star= T_star, R_star= R_star) * u.K
    if wave is not None:
        w = (wave/1.e8).value*u.cm  # angstroms to cm
        nu = c / w
    else:
        w = (c / nu)
    # constants appropriate to cgs units.
    
    fnu_term1 = 2 * np.pi * h * nu**3 / c**2
    fnu_term2 = np.exp((h*nu)/(k*temp))
    fnu = (fnu_term1 * (1/(fnu_term2 - 1))).value * u.erg / u.s / u.Hz / u.cm**2

    return np.array(fnu)

def planck_disk(r,
                wave = None,
                nu = None,
                T_star : float = 4000,
                R_star : float = 2.5):
    """Calculate planck function for given temperature(K) 

    Args:
        r : distance of the disk from the center (solar radius)
        temp : Effective temperature(K)
        wave (optional): wavelength in Angstron unit
        nu (optional): Frequency in Hz unit

    Returns:
        fnu: planck functiond in f_nu [erg/s/Hz/cm2]
    """
    temp = Temperature_disk(r = r, T_star= T_star, R_star= R_star) * u.K
    if wave is not None:
        w = (wave/1.e8).value*u.cm  # angstroms to cm
        nu = c / w
    else:
        w = (c / nu)
    # constants appropriate to cgs units.
    
    fnu_term1 = 2 * np.pi * h * nu**3 / c**2
    fnu_term2 = np.exp((h*nu)/(k*temp))
    fnu = (fnu_term1 * (1/(fnu_term2 - 1))).value * u.erg / u.s / u.Hz / u.cm**2

    return np.array(fnu)
def luminosity_disk(r,
                    wave = None,
                    nu = None,
                    T_star : float = 4000,
                    R_star : float = 2.5):
    """
    Get luminosity of the disk at the certain distance r (Solar radius unit)

    Args:
        r (_type_): distance of the disk from the center (solar radius)
        wave (_type_, optional): wavelength range for calculation
        nu (_type_, optional): Frequency range for calculation
        T_star (float, optional): Central star temperature(effective)
        R_star (float, optional): Central star radius(solar radius)

    Returns:
        _type_: _description_
    """
    term1 = (2 * pi**2 * r_sun**2).value
    term2 = planck_disk(r = r, wave = wave, nu = nu, T_star = T_star, R_star= R_star)
    term3 = r
    return term1 * term2 * term3 

#%%
radius_range_min = 12.5 # solar radius unit
radius_range_max = 5e4 # solar radius unit
radius_bin = 1
flux = []
flux_tot = np.zeros(len(nulist))
for radius in np.arange(radius_range_min, radius_range_max, radius_bin):
    f = luminosity_disk(r = radius, nu = nulist)
    flux_tot += f
    print(radius)

#%% 
# Prob 4.2
logx = np.log10(nulist.value)
logy = np.log10(nulist.value*flux_tot/const.L_sun.cgs.value)
plt.figure(figsize = (8, 6))
plt.plot(logx, logy, c ='k')
plt.xlim(10, 15)
plt.ylim(-12, 2)
plt.ylabel(r'$\nu L_\nu/L_\odot$', fontsize = 15)
plt.xlabel(r'$\log{\nu}$', fontsize = 15)

#%% 
# Prob 4.3
idx_alpha = (logx < 12.5) & (logx > 11.5)
A = logx[idx_alpha]
B = logy[idx_alpha]
from scipy.optimize import curve_fit
def func(x, a, b):
    return a*x + b
tt = curve_fit(func, A, B)
param = tt[0]
plt.figure(figsize = (8, 6))
plt.plot(logx, logy, c ='k')
plt.xlim(10, 15)
plt.ylim(-12, 2)
plt.ylabel(r'$\nu L_\nu/L_\odot$', fontsize = 15)
plt.xlabel(r'$\log{\nu}$', fontsize = 15)
plt.plot(logx, func(logx, param[0], param[1]), linestyle = '--', c='k', label = r'$\nu L_\nu \propto \nu^{1.326}$')
plt.legend(loc = 2, fontsize = 15)
#%%
# Prob 4.4
flux_star = (4*pi*(R_star*r_sun)**2*planck(temp = 4000*u.K, nu = nulist)).value
logy_star = np.log10(nulist.value*flux_star/const.L_sun.cgs.value)
logx = np.log10(nulist.value)
logy = np.log10(nulist.value*flux_tot/const.L_sun.cgs.value)
logy_both = np.log10((nulist.value*flux_tot+nulist.value*flux_star)/const.L_sun.cgs.value)

plt.figure(figsize = (8, 6))
plt.plot(logx, logy, c ='k', linestyle = '--', label = 'Disk')
plt.plot(logx, logy_star, c='r', linestyle = '--', label = 'Central star')
plt.plot(logx, logy_both, c='r', label = 'Total flux (Central star + Disk')
plt.xlim(10.5, 14.7)
plt.ylim(-6, 2)
plt.ylabel(r'$\nu L_\nu/L_\odot$', fontsize = 15)
plt.xlabel(r'$\log{\nu}$', fontsize = 15)
plt.plot(logx, func(logx, param[0], param[1]), linestyle = '--', c='b', label = r'$\nu L_\nu \propto \nu^{1.326}$')
plt.legend()
#%%