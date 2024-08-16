


#%%
from astropy import constants as const
from astropy import units as u
import numpy as np
m_e = 9.1094e-28 * u.g # g
e = 4.8032e-10 * u.cm**(3/2) * u.g**(1/2) * u.s**(-1)
#%%

#(3)
p = 3.3503e-2 *u.s
dpdt = 4.21e-13 *u.s / u.s
mass = 1.4 * const.M_sun.cgs
radius = (12 * u.km).to(u.cm)
omega = 2*np.pi/p
domegadt = -omega**2/(2*np.pi)*dpdt
T = -omega/domegadt
tc = -omega/domegadt/2
tc_year = tc/86400/365

comp1 = 12/5
comp2 = mass * const.c.cgs**3
comp3 = radius**4 * omega**2 * T
alpha = comp1 * comp2 / comp3
Bp_alpha = (alpha)**0.5

comp1 = 2/5 * mass * radius **2
dEdt = comp1 * omega * domegadt

# %% Problem 4
#%% Prob3
f12 = 0.4162
lamb_lyalpha_AA = 1215.67 * u.AA# AA
lamb_lyalpha_cm = lamb_lyalpha_AA.to(u.cm)
# (1)
A21 = (0.6670e16 * 2/8 * 0.4162 / lamb_lyalpha_AA**2).value / u.s

# (2)
comp1 = 8*np.pi**2*e**2
comp2 = 3*m_e*const.c.cgs*(lamb_lyalpha_cm)**2
gamma = comp1/comp2
cen_pi = 4/gamma
comp3 = np.pi * e**2
comp4 = m_e*const.c.cgs
ans_32 = comp3/comp4*f12*cen_pi

# %%
