

#%%
from companion_interaction_K10 import CompanionInteractionK10
from DOM_interaction_L17 import DOMInteractionL17
from astropy import constants as const
import numpy as np
from convert_AB_Vega import ABVegaMagnitude
#%%
radius = 1.4
comp = CompanionInteractionK10(rstar = radius, kappa = 0.2, m_wd = 1.2)
dom = DOMInteractionL17(E_exp = 0.9, M_ej = 1.0, kappa = 0.05, M_dom = 0.09, V_dom = 5000, f_dom = 0.07,  t_delay = 1e3, f_comp = 1.5 )
# %%
td = np.arange(0.1, 15, 0.1)
#%%
result_dom, _, _ = dom.calc_magnitude(td = td, filterset = 'UBVRIugri', visualize = False)
#%%
result_comp, _, _ = comp.calc_magnitude(td = td, filterset = 'UBVRIugri')
result_comp_U_vega = ABVegaMagnitude(result_comp['U'], magsys = 'AB', filter_ = 'U').AB
result_comp_B_vega = ABVegaMagnitude(result_comp['B'], magsys = 'AB', filter_ = 'B').AB
result_comp_V_vega = ABVegaMagnitude(result_comp['V'], magsys = 'AB', filter_ = 'V').AB
result_comp_R_vega = ABVegaMagnitude(result_comp['R'], magsys = 'AB', filter_ = 'R').AB
result_comp_I_vega = ABVegaMagnitude(result_comp['I'], magsys = 'AB', filter_ = 'I').AB
result_comp_g_vega = ABVegaMagnitude(result_comp['I'], magsys = 'AB', filter_ = 'g').AB
result_comp_r_vega = ABVegaMagnitude(result_comp['I'], magsys = 'AB', filter_ = 'r').AB

# %%
import matplotlib.pyplot as plt
# %%
plt.figure(dpi = 300)
plt.plot(result_dom['phase'], result_dom['Temperature_eff'], c ='r', label = 'DEI')
plt.plot(result_comp['phase'], result_comp['Temperature_eff'], c = 'b', label = 'CEI')
plt.ylabel(r'$T_{effective}$', fontsize = 15)
plt.xlabel('Days since explosion', fontsize = 15)

plt.legend()
#%%
from HHsupport_analysis import load_filt_keys
color_key, offset_key, _, _, label_key = load_filt_keys()
plt.figure(dpi = 300)
plt.gca().invert_yaxis()
for filter_ in 'B':
    plt.plot(result_dom['phase'], result_dom[filter_] + offset_key[filter_], c =color_key[filter_], linestyle = '--', label = label_key[filter_])
    plt.plot(result_comp['phase'], result_comp[filter_] + offset_key[filter_], c =color_key[filter_], linestyle = ':')
plt.ylim(-13, -18)
# %%
mej_range = np.arange(1, 1.41, 0.05)
eexp_range = np.arange(1, 1.51, 0.05)
# %%
kappa = 0.03
V_dom = 5000
f_dom = 0.15
result_DOM_U_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
result_DOM_B_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
result_DOM_V_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
result_DOM_R_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
result_DOM_I_vega = np.zeros([len(mej_range), len(eexp_range), len(td)])
tot_lumlist = np.zeros([len(mej_range), len(eexp_range)])
for i, mej in enumerate(mej_range):
    print(i)
    for j, eexp in enumerate(eexp_range):
        DOM = DOMInteractionL17(E_exp=eexp, M_ej=mej, kappa = kappa, M_dom = 0.1, V_dom = V_dom, t_delay = 5e3, f_dom=f_dom, f_comp = 1.5)
        result, _, _ = DOM.calc_magnitude(td = td, filterset = 'UBVRI')
        U_mag = ABVegaMagnitude(result['U'], magsys ='AB', filter_ ='U').vega
        B_mag = ABVegaMagnitude(result['B'], magsys ='AB', filter_ ='B').vega
        V_mag = ABVegaMagnitude(result['V'], magsys ='AB', filter_ ='V').vega
        R_mag = ABVegaMagnitude(result['R'], magsys ='AB', filter_ ='R').vega
        I_mag = ABVegaMagnitude(result['I'], magsys ='AB', filter_ ='I').vega
        tot_lum = np.sum(result['Luminosity_shock'])
        result_DOM_U_vega[i, j] = np.array(U_mag)
        result_DOM_B_vega[i, j] = np.array(B_mag)
        result_DOM_V_vega[i, j] = np.array(V_mag)
        result_DOM_R_vega[i, j] = np.array(R_mag)
        result_DOM_I_vega[i, j] = np.array(I_mag)
        tot_lumlist[i, j] = tot_lum
#%%
U_max = np.max(np.max(result_DOM_U_vega, axis=1),axis= 0)
B_max = np.max(np.max(result_DOM_B_vega, axis=1),axis= 0)
V_max = np.max(np.max(result_DOM_V_vega, axis=1),axis= 0)
R_max = np.max(np.max(result_DOM_R_vega, axis=1),axis= 0)
I_max = np.max(np.max(result_DOM_I_vega, axis=1),axis= 0)

U_min = np.min(np.min(result_DOM_U_vega, axis=1),axis= 0)
B_min = np.min(np.min(result_DOM_B_vega, axis=1),axis= 0)
V_min = np.min(np.min(result_DOM_V_vega, axis=1),axis= 0)
R_min = np.min(np.min(result_DOM_R_vega, axis=1),axis= 0)
I_min = np.min(np.min(result_DOM_I_vega, axis=1),axis= 0)
#%%
U_min = result_DOM_U_vega[0][0]
B_min = result_DOM_B_vega[0][0]
V_min = result_DOM_V_vega[0][0]
U_max = result_DOM_U_vega[-1][-1]
B_max = result_DOM_B_vega[-1][-1]
V_max = result_DOM_V_vega[-1][-1]


#%% Light curve
plt.figure(dpi = 300)
alpha = 0.3


#plt.fill_between(td, y1 = I_min-4, y2 = I_max-4, color ='k', alpha = alpha, label = 'I-4')
#plt.fill_between(td, y1 = R_min-3, y2 = R_max-3, color ='r', alpha = alpha, label = 'R-3')
plt.fill_between(td, y1 = V_min-2, y2 = V_max-2, color ='g', alpha = alpha, label = 'V-2')
plt.fill_between(td, y1 = B_min, y2 = B_max, color ='b', alpha = alpha, label = 'B')
plt.fill_between(td, y1 = U_min+2, y2 = U_max+2, color ='purple', alpha = alpha, label = 'U+2')
plt.plot(td, result_comp_U_vega+2, linestyle ='--', c ='purple')
plt.plot(td, result_comp_B_vega, linestyle ='--', c ='b')
plt.plot(td, result_comp_V_vega-2, linestyle ='--', c ='g')
#plt.plot(td, result_comp_R_vega-3, linestyle ='--', c ='r')
#plt.plot(td, result_comp_I_vega-4, linestyle ='--', c ='k')
plt.xlim(0, 10)
plt.ylim(-13, -21)
plt.grid()
plt.legend(loc = 1)
plt.ylabel('Absolute magnitude [Vega]')
plt.xlabel('Time from explosion [days]')
# %% Color
plt.figure(dpi = 300)
plt.plot(td, result_comp_U_vega - result_comp_B_vega)
for result_DOM_B_vega_single_eexp, result_DOM_V_vega_single_eexp in zip(result_DOM_U_vega, result_DOM_B_vega):
    for result_DOM_B_vega_single, result_DOM_V_vega_single in zip(result_DOM_B_vega_single_eexp, result_DOM_V_vega_single_eexp):
        plt.plot(td, result_DOM_B_vega_single - result_DOM_V_vega_single, c='k', alpha = 0.3)


#%%
from astropy.io import ascii
home_dir = '/home/hhchoi1022/Desktop/Gitrepo/Research/Supernova/DOM_model/'
range_E_exp = [0.5]#[0.5, 0.75, 1.0, 1.25, 1.5]
range_M_ej = [1.0, 1.2, 1.4]
range_kappa = [0.03, 0.05]
range_t_delay = [5e3, 7.5e3, 1e4]
range_f_comp = [1.5]
range_M_dom = [0.01, 0.04, 0.07, 0.10]
range_v_dom = [5e3, 7.5e3, 1e4]
range_f_dom = [0.05, 0.10, 0.15]

range_E_exp = [0.5]
range_M_ej = [1.0]
range_kappa = [0.03, 0.05]
range_t_delay = [1e4]
range_f_comp = [1.5]
range_M_dom = [0.1]
range_v_dom = [1e4]
range_f_dom = [0.10]

td = np.arange(0.1, 15, 0.1)
result_all = dict()
results_U = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_B = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_V = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_R = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_I = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_u = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_g = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_r = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])
results_i = np.zeros([len(range_E_exp), len(range_M_ej), len(range_kappa), len(range_t_delay), len(range_f_comp), len(range_M_dom), len(range_v_dom), len(range_f_dom), len(td) ])

for i1, E_exp in enumerate(range_E_exp):
    for i2, M_ej in enumerate(range_M_ej):
        for i3, kappa in enumerate(range_kappa):
            for i4, t_delay in enumerate(range_t_delay): 
                for i5, f_comp in enumerate(range_f_comp):
                    for i6, M_dom in enumerate(range_M_dom):
                        for i7, v_dom in enumerate(range_v_dom):
                            for i8, f_dom in enumerate(range_f_dom):
                                result = ascii.read(f'{home_dir}{E_exp}_{M_ej}_{kappa}_{t_delay}_{f_comp}_{M_dom}_{v_dom}_{f_dom}.dat', format = 'fixed_width')
                                results_U[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['U'])
                                results_B[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['B'])
                                results_V[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['V'])
                                results_R[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['R'])
                                results_I[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['I'])
                                results_u[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['u'])
                                results_g[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['g'])
                                results_r[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['r'])
                                results_i[i1, i2, i3, i4, i5, i6, i7, i8] = np.array(result['i'])
                                

# %%
max_U = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_U, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_B = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_B, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_V = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_V, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_R = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_R, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_I = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_I, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_u = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_u, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_g = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_g, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_r = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_r, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
max_i = np.max(np.max(np.max(np.max(np.max(np.max(np.max(np.max(results_i, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)

min_U = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_U, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_B = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_B, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_V = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_V, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_R = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_R, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_I = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_I, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_u = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_u, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_g = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_g, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_r = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_r, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
min_i = np.min(np.min(np.min(np.min(np.min(np.min(np.min(np.min(results_i, axis=0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0),axis= 0)
#%%
max_U = ABVegaMagnitude(magnitude = max_U, magsys = 'AB', filter_ = 'U').AB
max_B = ABVegaMagnitude(magnitude = max_B, magsys = 'AB', filter_ = 'B').AB
max_V = ABVegaMagnitude(magnitude = max_V, magsys = 'AB', filter_ = 'V').AB
max_g = ABVegaMagnitude(magnitude = max_g, magsys = 'AB', filter_ = 'g').AB
max_r = ABVegaMagnitude(magnitude = max_r, magsys = 'AB', filter_ = 'r').AB

min_U = ABVegaMagnitude(magnitude = min_U, magsys = 'AB', filter_ = 'U').AB
min_B = ABVegaMagnitude(magnitude = min_B, magsys = 'AB', filter_ = 'B').AB
min_V = ABVegaMagnitude(magnitude = min_V, magsys = 'AB', filter_ = 'V').AB
min_g = ABVegaMagnitude(magnitude = min_g, magsys = 'AB', filter_ = 'g').AB
min_r = ABVegaMagnitude(magnitude = min_r, magsys = 'AB', filter_ = 'r').AB

#%% Light curve
plt.figure(dpi = 300)
alpha = 0.3
"""
for i1, E_exp in enumerate(range_E_exp):
    for i2, M_ej in enumerate(range_M_ej):
        for i3, kappa in enumerate(range_kappa):
            for i4, t_delay in enumerate(range_t_delay): 
                for i5, f_comp in enumerate(range_f_comp):
                    for i6, M_dom in enumerate(range_M_dom):
                        for i7, v_dom in enumerate(range_v_dom):
                            for i8, f_dom in enumerate(range_f_dom):
                                plt.plot(td, results_U[i1, i2, i3, i4, i5, i6, i7, i8], linewidth = 1, alpha = 0.3, c='cyan')
                                plt.plot(td, results_B[i1, i2, i3, i4, i5, i6, i7, i8], linewidth = 1, alpha = 0.3, c='b')
                                plt.plot(td, results_V[i1, i2, i3, i4, i5, i6, i7, i8], linewidth = 1, alpha = 0.3, c='g')
"""
plt.fill_between(td, y1 = min_U-3, y2 = max_U-3, color ='purple', alpha = alpha, label = 'U-3')
plt.fill_between(td, y1 = min_B, y2 = max_B, color ='b', alpha = alpha, label = 'B')
plt.fill_between(td, y1 = min_V+3, y2 = max_V+3, color ='g', alpha = alpha, label = 'V+3')
#plt.fill_between(td, y1 = min_g-2.5, y2 = max_g-2.5, color ='k', alpha = 0.7, label = 'g')
#plt.fill_between(td, y1 = min_r-1.5, y2 = max_r-1.5, color ='r', alpha = 0.7, label = 'r')

plt.plot(td, result_comp_U_vega-3, linestyle ='--', c ='purple')
plt.plot(td, result_comp_B_vega, linestyle ='--', c ='b')
plt.plot(td, result_comp_V_vega+3, linestyle ='--', c ='g')
#plt.plot(td, result_comp_g_vega, linestyle ='--', c ='k')
#plt.plot(td, result_comp_r_vega, linestyle ='--', c ='r')
#plt.plot(td, result_comp_R_vega, linestyle ='--', c ='b')
#plt.plot(td, result_comp_I_vega, linestyle ='--', c ='b')

plt.xlim(0, 10)
plt.ylim(-5, -24)
plt.grid()
plt.ylabel('Absolute magnitude [Vega]')
plt.xlabel('Time from explosion [days]')
plt.legend()
# %%
sample_parameters = dict(E_exp = 1.4,
                         M_ej = 1.0,
                         kappa = 0.05,
                        
                         M_dom = 0.11,
                         V_dom = 5e3,
                         f_dom = 0.14,
                         
                         t_delay = 200,
                         f_comp = 1.5)
model = DOMInteractionL17(**sample_parameters)
t_range = np.arange(0.01, 10, 0.1)
DEI= model.calc_magnitude(td = t_range, filterset = 'UBgVri', visualize = True)
DEI_vel = []
for t in t_range:
    DEI_vel.append(model._velocity_diffusion(t = 86400*t))

#%%
model =  CompanionInteractionK10(rstar = 2.1, m_wd = 1.2, v9 = 1, commonangle = False)
t_range = np.arange(0.01, 10, 0.1)
CEI = model.calc_magnitude(td = t_range, filterset = 'UBgVri', visualize = True)


# How to derive the mininum index of multi-dimensional numpy array?
# %%
import matplotlib.pyplot as plt
#plt.plot(t_range, DEI[0]['Luminosity_shock'], c ='k')
#plt.plot(t_range, CEI[0]['Luminosity_shock'], c= 'r')
plt.plot(t_range, DEI[0]['Temperature_eff'], c ='k', label = 'DEI')
plt.plot(t_range, CEI[0]['Temperature_eff'], c= 'r', label = 'CEI')
plt.ylim(0, 25000)
plt.legend()
#plt.yscale('log')
#plt.gca().invert_yaxis()
# %%
plt.plot(t_range, DEI[0]['Luminosity_shock'], c ='k')
plt.plot(t_range, CEI[0]['Luminosity_shock'], c= 'r')
plt.yscale('log')
# %%

plt.plot(t_range, np.array(DEI_vel)/1e2/1e3)

#%%
from astropy.io import ascii
ashall = ascii.read('../../Data/SN2021aefx/observation/velocity_ejecta/ashall.csv')
ashall['phase'] = [59529.87, 59531.27, 59532.08]
ni_hv = ascii.read('../../Data/SN2021aefx/observation/velocity_ejecta/ni_hv.csv')
ni_pv = ascii.read('../../Data/SN2021aefx/observation/velocity_ejecta/ni_pv.csv')
ni_hv['phase'] = [59529.86, 59530.87, 59531.12, 59531.70, 59532.07, 59532.84, 59533.63, 59533.86, 59534.50, 59535.07, 59536.45, 59536.86, 59538.07, 59538.85]
ni_hv['error'] = np.array([0.11, 0.07, 0.12, 0.08, 0.03, 0.03, 0.04, 0.03, 0.04, 0.03, 0.05, 0.04, 0.03, 0.03]) * 10000
ni_pv['phase'] = [59529.86, 59530.87, 59531.12, 59531.70, 59532.07, 59532.84, 59533.63, 59533.86, 59534.50, 59535.07, 59536.45, 59536.86, 59538.07, 59538.85]
ni_pv['error'] = np.array([0.29, 0.13, 0.25, 0.15, 0.06, 0.06, 0.07, 0.05, 0.05, 0.02, 0.02, 0.01, 0.01, 0.01]) * 10000
# %%
effective_range = 59530.349
t_exp = 59528.53
#plt.scatter(ashall['phase']-t_exp, ashall['v'])
plt.figure(dpi = 300, figsize = (4,3))
plt.plot(t_range, np.array(DEI_vel)/1e2/1e3/1e3, c = 'k', linestyle = '--', label = 'DOM model')
plt.scatter(ni_hv['phase']-t_exp, ni_hv['v']/1e3, edgecolor ='r' , facecolor = 'none', label = 'HVF[Si II]')
plt.errorbar(ni_hv['phase']-t_exp, ni_hv['v']/1e3, yerr = ni_hv['error']/1e3, fmt = 'None', c ='r')
plt.scatter(ni_pv['phase']-t_exp, ni_pv['v']/1e3, edgecolor ='b' , facecolor = 'none', label = 'PVF[Si II]')
plt.errorbar(ni_pv['phase']-t_exp, ni_pv['v']/1e3, yerr = ni_pv['error']/1e3, fmt = 'None', c ='b')
plt.yticks([10, 20, 30, 40, 50, 60], [10, 20, 30, 40, 50, 60])
plt.ylim(5, 55)
plt.xlim(0, 10)
plt.xlabel(r'Days since $t_{exp, DOM}$')
plt.ylabel('Ejecta velocity [km/s]')
plt.legend()
#plt.yscale('log')
# %%
t_range = np.arange(0.001, 10, 0.001) * 86400
DEI_rad = []
DOM_rad = []
DOM_start_rad = 
for t in t_range:
    DEI_rad.append(model._radius_photosphere(t = t))
# %%

DEI_rad
# %%
