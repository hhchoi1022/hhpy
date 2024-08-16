#%%
import numpy as np
from astropy import constants as const
import astropy.units as u
from scipy import integrate
#%%
c = const.c.cgs
h = const.h.cgs
k = const.k_B.cgs
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
    
    fnu_term1 = 2 * h * nu**3 / c**2
    fnu_term2 = np.exp((h*nu)/(k*temp))
    fnu = (fnu_term1 * (1/(fnu_term2 - 1))).value * u.erg / u.s / u.Hz / u.cm**2

    return np.array(fnu)

# %%


#%% Problem 1
# 1.1
nu = 2e2 * 1e6 * u.Hz # 200MHz
temp = 1e4 * u.K # 10000K
distance = 500 * u.pc # 500pc
BB_intensity = planck(temp = temp, nu = nu) #erg/s/cm2/Hz
ff_flux_density = (1e2* u.Jy.to(u.erg/u.s/u.Hz/u.cm**2)) #erg/s/cm2/Hz
radius_500 = np.sqrt(ff_flux_density/(np.pi*BB_intensity))
radius = distance * radius_500 # radius in pc
# 1.2
factor1 = 3.37e-7 * (1) * (1) * (1) * (1)
tau = 1
EM = tau / factor1 * u.pc / u.cm**6
n_e = np.sqrt(EM / radius / (4/3))
# 1.3
radius_cm = radius.to(u.cm)
comp1 = 4 * np.pi * radius_cm**3
comp2 = k * temp * n_e
E_thermal = comp1  * comp2
comp1 = (4/3) * np.pi * radius_cm**3
comp2 = 1.42e-27
comp3 = n_e**2 * (temp)**0.5 * 1.3
E_ff = (comp1 * comp2 * comp3).value * u.erg / u.s
time_cooling = E_thermal /E_ff
time_cooling_yr = time_cooling / 86400 / 365

#%% Problem 3
mass_sun = const.M_sun.cgs
radius_sun = const.R_sun.cgs
mass_h = 1.67e-24 * u.g

mass_loss_rate_wind = (7e-11 * mass_sun / u.yr).to(u.g/u.s)
velocity_wind = (650 * u.km / u.s).to(u.cm/u.s)
radius_star = 0.75 * radius_sun
distance = (3.2 * u.pc).to(u.cm)
temperature_wind = 1e4 * u.K
nu_range = np.logspace(8, 12, num = 200) / u.s
p_min = radius_star
p_max = 300* radius_star

#%%
def g_ff(nu):
    nu9 = nu/1e9
    temp4 = temperature_wind/1e4
    return (6.155 * (nu9)**-0.118 * temp4**0.177).value
    
def p_c(nu):
    comp1 = 1.85e8
    comp2 = temperature_wind**-0.5
    comp3 = g_ff(nu)
    comp4 = nu**-3
    comp5 = (1 - np.exp(-h*nu/k/temperature_wind))
    comp6 = np.pi
    comp7 = mass_loss_rate_wind / (4 * np.pi * mass_h * velocity_wind)
    p_c = (comp1 * comp2 * comp3 * comp4 * comp5* comp6 * comp7**2)**(1/3)
    return p_c.value * u.cm

def tau_nu_out(nu, p):
    comp1 = 3.69e8
    comp2 = temperature_wind**(-0.5)
    comp3 = g_ff(nu)
    comp4 = nu**-3
    comp5 = (1 - np.exp(-h*nu/k/temperature_wind))
    comp6 = np.pi / 2 / (p)**3
    comp8 = mass_loss_rate_wind / (4 * np.pi * mass_h * velocity_wind)
    tau = (comp1 * comp2 * comp3 * comp4 * comp5* comp6 * comp8**2).value
    return tau

def flux_nu_out(nu):
    comp1 = planck(temp = temperature_wind, nu = nu)
    def flux_monochrome(p, nu):
        p *= u.cm
        nu *= u.Hz
        comp2 = 1 - np.exp(-tau_nu_out(nu = nu, p = p))
        comp3 = 2 * np.pi * p
        return (comp2 * comp3).value
    integrated_value = integrate.quad(flux_monochrome, p_min.value, p_max.value, args = (nu.value))
    flux_nu_integrated = comp1 * integrated_value
    return (flux_nu_integrated / distance**2).value

def tau_nu_in(nu, p):
    comp1 = 3.69e8
    comp2 = temperature_wind**(-0.5)
    comp3 = g_ff(nu)
    comp4 = nu**-3
    comp5 = (1 - np.exp(-h*nu/k/temperature_wind))
    comp6 = np.pi / 2 / (p)**3
    comp8 = mass_loss_rate_wind / (4 * np.pi * mass_h * velocity_wind)
    comp9 = np.pi /2
    comp10 = p * np.sqrt(radius_star**2 - p**2)/radius_star**2
    comp11 = np.arctan(np.sqrt(radius_star**2 - p**2)/p)
    tau = (comp1 * comp2 * comp3 * comp4 * comp5* comp6 * comp8**2 * (comp9 - comp10.value - comp11.value)).value
    return tau

def flux_nu_in(nu):
    comp1 = planck(temp = temperature_wind, nu = nu)
    def flux_monochrome(p, nu):
        p *= u.cm
        nu *= u.Hz
        comp2 = 1 - np.exp(-tau_nu_in(nu = nu, p = p))
        comp3 = 2 * np.pi * p
        return (comp2 * comp3).value
    integrated_value = integrate.quad(flux_monochrome, 0, p_min.value,  args = (nu.value))
    flux_nu_integrated = comp1 * integrated_value
    return (flux_nu_integrated / distance**2).value
#%%
# Figure 1
def function(nu, a, b):
    return b* nu ** a
fit_result_pc = curve_fit(function, nu_range, p_c(nu_range)/radius_star)

plt.figure(dpi = 300, figsize = (4,3))
plt.plot(nu_range, p_c(nu_range)/radius_star, c= 'k')
plt.scatter(8.64e10, function(8.64e10, fit_result_pc[0][0], fit_result_pc[0][1]), c='r', marker = '*')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\nu \ [Hz]$')
plt.ylabel(r'$p_c \ [R_\star]$')

#%%
function(8.64e10, fit_result_pc[0][0], fit_result_pc[0][1])
#%%
# Calculation of the flux
nu_range = np.logspace(8, 12, num = 1000) / u.s
flux_out_all = []
flux_in_all = []
flux_both_all = []
for nu in nu_range:
    flux_out = flux_nu_out(nu)[0]
    flux_in = flux_nu_in(nu)[0]
    flux_both = flux_out + flux_in
    flux_out_all.append(flux_out)
    flux_in_all.append(flux_in)
    flux_both_all.append(flux_both)
#%%
# Plot
plt.figure(dpi = 300, figsize = (4,3))
plt.plot(nu_range, np.array(flux_out_all) * 1e23, c='k', linestyle = '--', label = r'$F_{out}$')
plt.plot(nu_range, np.array(flux_in_all) * 1e23, c='k', linestyle = ':', label = r'$F_{in}$')
plt.plot(nu_range, np.array(flux_both_all) * 1e23, c ='r', label = r'$F_{both}$')
#plt.plot(nu_range, function(nu_range, x0 = 10**5.3))
plt.xscale('log')
plt.yscale('log')
plt.ylim(5e-6, 1e-3)
plt.xlabel(r'$\nu \ [Hz]$')
plt.ylabel(r'$F \ [Jy]$')
plt.legend()

#%%
from scipy.optimize import curve_fit
#%%
def function(nu, a, b):
    return b* nu ** a

#%%
# Fitting
min_idx =500
max_idx =550
flux_both_all = np.array(flux_both_all)
fit_result = curve_fit(function, nu_range[min_idx:max_idx], flux_both_all[min_idx:max_idx]* 1e23)
# Plot
plt.figure(dpi = 300, figsize = (4,3))
#plt.plot(nu_range, np.array(flux_out_all) * 1e23, c='k', linestyle = '--', label = r'$F_{out}$')
#plt.plot(nu_range, np.array(flux_in_all) * 1e23, c='k', linestyle = ':', label = r'$F_{in}$')
plt.plot(nu_range, np.array(flux_both_all) * 1e23, c ='k', label = r'$F_{both}$')
plt.plot(nu_range[min_idx:max_idx], flux_both_all[min_idx:max_idx]*1e23, c= 'r', label = 'Fit region')
plt.plot(nu_range, function(nu_range, a = fit_result[0][0], b = fit_result[0][1]), linestyle = ':', c='k', label = rf'$F_\nu \ \propto$ {np.round(fit_result[0][0], 3)}')
plt.axvline(8.64e10, c ='r', linestyle = '--')

#plt.plot(nu_range, function(nu_range, x0 = 10**5.3))
plt.xscale('log')
plt.yscale('log')
plt.ylim(5e-6, 1e-3)
plt.xlabel(r'$\nu \ [Hz]$')
plt.ylabel(r'$F \ [Jy]$')
plt.legend()
#%% Problem 4
# 4.2
mass_e = 9.11e-28 # g
sigma_thom = (6.65e-29 * u.m**2 ).to(u.cm**2)
c = const.c.cgs
e_charge = 4.80e-10 #cm(3/2)g(1/2)s(-1)
nu = 1e9
obs_flux = 2720e-23
distance = 3.4e3*u.pc.to(u.cm)
radius = 2.47*u.pc.to(u.cm)
# Observed flux
comp1 = 2/9
comp2 = radius
comp3 = distance
comp4 = sigma_thom
comp5 = c
comp6 = 8*np.pi
comp7 = (2*np.pi*mass_e*c)
comp8 = e_charge
comp9 = nu
A = comp1 * comp2**3 / comp3**2 * comp4 * comp5 / comp6 * comp7**0.25 / comp8**0.25 * comp9**-0.75
C = obs_flux/A
# Calculate B 
comp1 = 7*np.pi*1/5*C
comp2 = mass_e
comp3 = c**2
B = (comp1 * comp2 * comp3)**(4/15)
#%%
# 4.3
comp1 = 15/7
comp2 = B
comp3 = 8*np.pi
comp4 = 4*np.pi/3
comp5 = radius
E_tot = comp1*comp2**2/comp3*comp4*comp5**3
# %%
nu_max = 1e12
comp1 = obs_flux * 1e9 / 0.25 
comp2 = (nu_max/1e9)
comp3 = 4*np.pi*distance**2
luminosity = comp1*comp2**0.25*comp3
t_cool = E_tot/luminosity
#%%
t_cool/86400/365
