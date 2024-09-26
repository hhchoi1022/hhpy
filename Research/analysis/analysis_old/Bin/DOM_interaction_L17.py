
#%%
import numpy as np
#%%

def planck(wave, temp):
    #if len(wave) > 1 :
    #    print('Syntax - bbflux = planck( wave, temp)')
    #    return 0
    #if len(temp) != 1 :
    #    input('Enter a blackbody temperature : ')
    # Gives the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/a
    w = wave/1.e8 # angstroms to cm
    # constants appropriate to cgs units.
    c1 = np.float128(3.7417749e-5)          # =2*!DPI*h*c*c
    c2 = np.float128(1.4387687)             # =h*c/k
    val = c2/w/np.float128(temp)
    bbflux = c1/( (w**5)*(np.exp(val)-1.))
    return bbflux*1.e-8 # convert to ergs cm-2 s-1 A-1
# %%
from scipy.optimize import fsolve
factor_compression = 1.5
factor_DOM = 0.1
mass_DOM = 0.03
kappa = 0.03
velocity_ejecta = 12000
velocity_interaction = 10600
velocity_DOM = 5000
velocity_diffusion = 10000
time_interaction = 1e4

# constants
c = 2.998e10 # cm/s
msun = 1.988e33 # g
rsun = 6.953e10 # cm
fraction_56_low = 0.5
fraction_56_high = 0.1
# variables
energy_explosion = 1.0e51 # erg
mass_ejecta = 1.4 * msun # g
mass_DOM = 0.05 * msun # g
kappa = 0.02 #cm^2g^-1
velocity_DOM = 5e8 # cm/s
fraction_DOM = 0.1

velocity_e = (energy_explosion / (6*mass_ejecta))**0.5 # cgs
A = mass_ejecta / (8*np.pi*velocity_e**3)
def velocity_diffusion(t):
    def func_velocity(x):
        return t**2 * (4*np.pi*c*velocity_e) / (3*kappa*mass_ejecta) - np.exp(-x) * (x+1)
    return velocity_e * fsolve(func_velocity, (10))

def mass_diffusion(t):
    component1 = (4*np.pi*c*velocity_e)/(3*kappa*mass_ejecta)*t**2
    component2 = (velocity_diffusion(t))/(2*velocity_e)*np.exp((-velocity_diffusion(t))/(velocity_e))
    return mass_ejecta*(component1 + component2)

def func_velocity_shock(x):
    component1 = (mass_DOM/fraction_DOM)*(x-velocity_DOM/velocity_e)
    component2 = (4*np.pi*velocity_e**3*np.exp(-x))
    component3 = (x**2+4*x+6)
    return component1 - component2*component3

velocity_shock = velocity_e * fsolve(func_velocity_shock, (10))
def shock_profile():
    component1 = A * mass_shock / mass_ejecta
    
#%%
velocity_diffusion(1e3)
#%%
func1 = lambda x: (4*np.pi*c*velocity_e) / (3*kappa*mass_ejecta) / np.exp(-x) / (1+x)


#%%
logLt = np.log10(2*np.pi) + np.log10(c) - np.log10(kappa) - np.log10(velocity_ejecta) - np.log10(factor_compression) + 2*np.log10(velocity_interaction-velocity_DOM) + np.log10(time_interaction) + 2*np.log10(velocity_diffusion)
logRph = np.log10(velocity_ejecta) + np.log10(factor_compression) + np.log10(np.log(3*))
# %%
