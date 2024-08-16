

#%%
from companion_interaction_K10 import CompanionInteractionK10
from DOM_interaction_L17 import DOMInteractionL17
from astropy import constants as const
import numpy as np
from convert_AB_Vega import ABVegaMagnitude
#%%
radius = (1e12/const.R_sun.cgs).value
comp = CompanionInteractionK10(rstar = radius, kappa = 0.2, m_wd = 1.4)
dom = DOMInteractionL17(E_exp = 1.5, M_ej = 1.4, kappa = 0.05, M_dom = 0.1, V_dom = 5000, f_dom = 0.1, t_delay = 5e3, f_comp = 1.5 )
# %%
td = 0.0 + np.arange(0, 10, 0.1)
#%%
result_dom = dom.calc_magnitude(td = td, filterset = 'UBVRI')
#%%
result_comp = comp.calc_magnitude(td = td, filterset = 'UBVRI')
result_comp_U_vega = ABVegaMagnitude(result_comp['U'], magsys = 'AB', filter_ = 'U').vega
result_comp_B_vega = ABVegaMagnitude(result_comp['B'], magsys = 'AB', filter_ = 'B').vega
result_comp_V_vega = ABVegaMagnitude(result_comp['V'], magsys = 'AB', filter_ = 'V').vega
# %%
import matplotlib.pyplot as plt
# %%
plt.plot(result_dom['phase'], result_dom['Temperature_eff'], c ='k')
plt.plot(result_comp['phase'], result_comp['Temperature_eff'])

# %%
mej_range = np.arange(1, 1.41, 0.05)
eexp_range = np.arange(1, 1.51, 0.05)
mej_grid, eexp_grid = np.meshgrid(mej_range, eexp_range)
# %%
kappa = 0.05
V_dom = 5000
result_DOM_U_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
result_DOM_B_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
result_DOM_V_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
for i, mej in enumerate(mej_range):
    print(i)
    for j, eexp in enumerate(eexp_range):
        DOM = DOMInteractionL17(E_exp=eexp, M_ej=mej, kappa = kappa, M_dom = 0.1, V_dom = V_dom, t_delay = 5e3, f_dom=0.1, f_comp = 1.5)
        result = DOM.calc_magnitude(td = td, filterset = 'UBV')
        U_mag = ABVegaMagnitude(result['U'], magsys ='AB', filter_ ='U').vega
        B_mag = ABVegaMagnitude(result['B'], magsys ='AB', filter_ ='B').vega
        V_mag = ABVegaMagnitude(result['V'], magsys ='AB', filter_ ='V').vega
        result_DOM_U_vega[i, j] = np.array(U_mag)
        result_DOM_B_vega[i, j] = np.array(B_mag)
        result_DOM_V_vega[i, j] = np.array(V_mag)
        #plt.plot(td, result['U']+2, c ='purple', alpha = 0.5)
        #plt.plot(td, result['B'], c ='b', alpha = 0.5)
        #plt.plot(td, result['V']-2, c ='g', alpha = 0.5)
        #plt.plot(td, result_comp['U']+2, linestyle ='--', c ='purple')
        #plt.plot(td, result_comp['B'], linestyle ='--', c ='b')
        #plt.plot(td, result_comp['V']-2, linestyle ='--', c ='g')
#%%
U_max = np.max(np.max(result_DOM_U_vega, axis=1),axis= 0)
B_max = np.max(np.max(result_DOM_B_vega, axis=1),axis= 0)
V_max = np.max(np.max(result_DOM_V_vega, axis=1),axis= 0)

U_min = np.min(np.min(result_DOM_U_vega, axis=1),axis= 0)
B_min = np.min(np.min(result_DOM_B_vega, axis=1),axis= 0)
V_min = np.min(np.min(result_DOM_V_vega, axis=1),axis= 0)
#%% Light curve
plt.figure(dpi = 300)
alpha = 0.3
plt.fill_between(td, y1 = U_min+2, y2 = U_max+2, color ='purple', alpha = alpha)
plt.fill_between(td, y1 = B_min, y2 = B_max, color ='b', alpha = alpha)
plt.fill_between(td, y1 = V_min-2, y2 = V_max-2, color ='g', alpha = alpha)
plt.plot(td, result_comp_U_vega+2, linestyle ='--', c ='purple')
plt.plot(td, result_comp_B_vega, linestyle ='--', c ='b')
plt.plot(td, result_comp_V_vega-2, linestyle ='--', c ='g')
plt.xlim(0, 10)
plt.ylim(-13, -21)
plt.grid()
plt.ylabel('Absolute magnitude [Vega]')
plt.xlabel('Time from explosion [days]')
# %% Color
plt.figure(dpi = 300)
plt.plot(td, result_comp_U_vega - result_comp_B_vega)
for result_DOM_B_vega_single_eexp, result_DOM_V_vega_single_eexp in zip(result_DOM_U_vega, result_DOM_B_vega):
    for result_DOM_B_vega_single, result_DOM_V_vega_single in zip(result_DOM_B_vega_single_eexp, result_DOM_V_vega_single_eexp):
        plt.plot(td, result_DOM_B_vega_single - result_DOM_V_vega_single, c='k', alpha = 0.3)

