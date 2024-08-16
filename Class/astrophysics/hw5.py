#%%
import astropy.constants as const
import numpy as np
import astropy.units as u

# %%
k = const.k_B.cgs
c = const.c.cgs
h = const.h.cgs
e_charge = 4.80e-10 #cm(3/2)g(1/2)s(-1)
mass_e = 9.11e-28 # g
sigma_thom = (6.65e-29 * u.m**2 ).to(u.cm**2)
a = 8 * np.pi**5 * k**4 / 15 / h**3 / c**3
# %% Prob1
# 1.2
nu_12 = 1.43 * 1e9 * u.Hz 
flux_12 = (72 * 1e-3 * u.Jy).to(u.erg/u.cm**2/u.Hz/u.s)
radius_1 = (0.0123 * u.pc).to(u.cm)
distance_1 = (3.6 * 1e6 * u.pc).to(u.cm)
comp1 = c**2
comp2 = 2 * np.pi * k * nu_12**2
comp3 = flux_12
comp4 = (radius_1/distance_1)**-2
Tb_temp1 = comp1 / comp2 * comp3 * comp4
# 1.3
flux_13 = (72 * 1e-3 * u.Jy).to(u.erg/u.cm**2/u.Hz/u.s)
nu_13 = 1.43 * 1e9 * u.Hz 
comp1 = np.pi
comp2 = (radius_1 / distance_1)**2
comp3 = 2 * nu_13**2.5
comp4 = mass_e
comp5 = (2*np.pi*mass_e*c)
comp6 = flux_13 * e_charge**0.5
B = (comp1 * comp2 * comp3 * comp4 * comp5**0.5 / comp6)**2
# 1.4
flux_14 = (25 * 1e-3 * u.Jy).to(u.erg/u.cm**2/u.Hz/u.s)
nu_14 = 23 * 1e9 * u.Hz 
comp1 = (1/4/np.pi)
comp2 = 2/3
comp3 = sigma_thom
comp4 = c
comp5 = (B**2)/8/np.pi
comp6 = nu_14**-1
comp7 = flux_14
comp8 = (4/3*np.pi*radius_1**3)
comp9 = distance_1**2
n0 = (comp1 * comp2 * comp3 * comp4 * comp5 * comp6)**-1 * comp7 / comp8 * comp9
# %% Prob 2
# 2.2
comp1 = 324 * k**5
comp2 = e_charge**2
comp3 = mass_e**6
comp4 = c**13
comp5 = 1e12**5 * 1e9
C = comp1 * comp2 / comp3 / comp4 * comp5
# 2.3
flux_23 = (2.19e4 * u.Jy).to(u.erg/u.cm**2/u.Hz/u.s)
nu_23 = 12.6e6 * u.Hz
comp1 = c**2
comp2 = flux_23
comp3 = 2*np.pi*k*nu_23**2
comp4 = 10 / 206265
T = comp1 * comp2 / comp3 / comp4**2
ratio = 0.415 * 0.61**5 * 0.0126
# %% Problem 3
temperature_3 = 2.7
U_gamma_3 = a * temperature_3**4
comp1 = 3*mass_e*c**2
comp2 = 4 * sigma_thom * c * U_gamma_3
factor_3 = comp1 / comp2/1/86400/365
#%% Problem 4
n_e = 1e-3
r = (1e6 * u.pc).to(u.cm)
temperature_4 = 1e8
tau = n_e * sigma_thom * 2 * r
y = 4*k*temperature_4/mass_e/c**2 * tau
delT = 2.8e-4 /2 * 2.73
# %%
delT
# %%
