#%%
#from astropy.constants import constant
from astropy import constants as const
import numpy as np
from astropy import units as u
#%%
R_sun = const.R_sun.cgs
L_sun = const.L_sun.cgs

R_earth = const.R_earth.cgs
AU = const.au.cgs 

G = const.G.cgs
c = const.c.cgs
pi = np.pi * u.sr
unit_flux = u.erg / u.s / u.cm**2
# %% Prob 1
#(1)
earth_solid_angle = pi*R_earth**2 / AU**2
earth_fraction = earth_solid_angle/(4*pi)
ans_11 = R_earth**2 / 4 / AU**2
#(2)
F_sun = 6.33e10 * unit_flux
ans_12 = F_sun/pi
#(3)
D_mars = 246*1e6*1e3*1e2 * u.cm
sun_solid_angle_from_mars = (pi * R_sun**2 / D_mars**2)
sun_solid_angle_from_earth = (pi * R_sun**2 / AU**2)
ans_13 = (sun_solid_angle_from_mars * ans_12 / c)
#%% Prob 2
#(3)
pc = const.pc.cgs
H0 = 70*u.km.to(u.cm)*u.cm / u.s /(1e6*pc)
rho_L = 2.3*1e8*L_sun/(1e6*pc)**3 
I_sky = rho_L * c / (4*pi*H0) / (4.25e10) *u.sr / u.arcsec**2
nu = (c/(5000*u.AA.to(u.cm)*u.cm))
#7.415 * 10**-16 / (3 * 10**10 / (5000 * 10**-8))
I_nu_sky = I_sky/nu
m_sky = -2.5*np.log10(I_nu_sky.value) - 48.6
# %% Prob 3
#(1)
ans_31 = 4e-15*10**(7.7)/0.23*((2e10)**0.23 -(2e7)**0.23) * u.erg / u.s / u.cm**2 / u.sr
#(2)
ans_32 = pi * ans_31 * (2/3.4e3)**2
#(3)
ans_33 = 4 * pi * (2*pc)**2 * pi * ans_31 / u.sr
# %%
