#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 19:58:41 2022
\
@author: hhchoi1022
"""
#%%

import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters, minimize, report_fit
from astropy.table import Table, vstack
from astropy.io import ascii
from HHsupport_analysis import mag_to_flux
from HHsupport_analysis import load_filt_keys
from HHsupport_analysis import flux_to_mag

#%% Load data
homedir = '../data/SN2021aefx/'
obsdir = homedir+'observation/lightcurve'

#IMSNG
all_tbl = Table.read(obsdir+'/IMSNG/IMSNG_hostmwextin3.10.dat', format = 'ascii.fixed_width')
filterlist = sorted(list(set(all_tbl['filter'])))
observatorylist = sorted(list(set(all_tbl['observatory'])))
UL_tbl = all_tbl[all_tbl['status'] == 'UL']
all_tbl = all_tbl[all_tbl['status'] == 'detected']
UL_tbl['absmag'] = UL_tbl['UL5_4'] - 31.18
all_tbl['absmag'] = all_tbl['mag'] - 31.18
UL_tbl_g, UL_tbl_r, UL_tbl_i = UL_tbl[UL_tbl['filter'] == filterlist[0]], UL_tbl[UL_tbl['filter'] == filterlist[2]], UL_tbl[UL_tbl['filter'] == filterlist[1]]
all_g, all_r, all_i = all_tbl[all_tbl['filter'] == filterlist[0]], all_tbl[all_tbl['filter'] == filterlist[2]], all_tbl[all_tbl['filter'] == filterlist[1]]
all_KCT, all_LSGT, all_RASA = all_tbl[all_tbl['observatory'] == observatorylist[0]], all_tbl[all_tbl['observatory'] == observatorylist[1]], all_tbl[all_tbl['observatory'] == observatorylist[2]],
KCT_g, LSGT_g = all_g[all_g['observatory'] == observatorylist[0]], all_g[all_g['observatory'] == observatorylist[1]]
KCT_r, LSGT_r, RASA_r = all_r[all_r['observatory'] == observatorylist[0]], all_r[all_r['observatory'] == observatorylist[1]], all_r[all_r['observatory'] == observatorylist[2]]
KCT_i, LSGT_i = all_i[all_i['observatory'] == observatorylist[0]], all_i[all_i['observatory'] == observatorylist[1]]

#Hosseinzadeh 2022
extern_tbl = ascii.read(obsdir+'/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat', format = 'fixed_width')
extern_tbl['absmag'] = extern_tbl['mag'] - 31.18
extern_tbl = extern_tbl[extern_tbl['observatory'] == 'LasCumbres1m']
extern_filterlist = sorted(list(set(extern_tbl['filter'])))
extern_g, extern_r, extern_i = extern_tbl[extern_tbl['filter'] == 'g'], extern_tbl[extern_tbl['filter'] == 'r'], extern_tbl[extern_tbl['filter'] == 'i']
extern_U, extern_B, extern_V = extern_tbl[extern_tbl['filter'] == 'U'], extern_tbl[extern_tbl['filter'] == 'B'], extern_tbl[extern_tbl['filter'] == 'V']

#Combine
all_extern_tbl = vstack([all_tbl, extern_tbl])
# %%

#%% Slicing the data for fitting
early_tbl = all_extern_tbl[(all_extern_tbl['obsdate'] > 59531)&(all_extern_tbl['obsdate'] < 59537)]
early_tbl.sort('obsdate')
early_tbl['flux'] = mag_to_flux(early_tbl['mag'])
early_tbl['e_flux'] = early_tbl['e_mag']*mag_to_flux(early_tbl['mag'])*2.303/2.5
filter_tbls = early_tbl.group_by('filter').groups
filter_key = list(set(early_tbl['filter']))
filter_key.sort()
fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
y_fit = [np.array((fit_table[filter_]['flux'].tolist())) for filter_ in filter_key]
e_y_fit = [np.array((fit_table[filter_]['e_flux'].tolist())) for filter_ in filter_key]

#%% Setting loss function
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux

def fireball_filter(params, time, filter_):
    exptime = params['exptime']
    amp = params[f'amp_{filter_}']
    alpha= params[f'alpha_{filter_}']
    return fireball_model(np.array(time), amp, exptime, alpha)

def calc_chisq(params, x, y, e_y, filter_key):
    tot_chisq = []
    for i, filter_ in enumerate(filter_key):
        obstime = x[i]
        obsflux = y[i]
        obserr = e_y[i]
        modelflux = fireball_filter(params, obstime, filter_)
        chisq = (((obsflux - modelflux)/obserr)**2)
        tot_chisq.append(chisq)
    return np.concatenate(tot_chisq)
#%% Setting parameters
fit_params = Parameters()
fit_params.add('exptime', value = 59529, min = 59525, max = 59535)
for filter_ in list(fit_table.keys()):
    fit_params.add(f'amp_{filter_}', value = 200, min = 0, max = 2000)
    fit_params.add(f'alpha_{filter_}', value = 2, min = 1, max = 10)
#%% Fitting
out = minimize(calc_chisq, fit_params, args = (x_fit, y_fit, e_y_fit, filter_key))
report_fit(out.params)

#%% Visualization
color_key, offset_key, _, _, label_key = load_filt_keys(filter_key)
plt.figure(dpi = 300)
plt.gca().invert_yaxis()
for filter_ in filter_key:
    exptime = out.params['exptime']
    amp = out.params[f'amp_{filter_}']
    alpha= out.params[f'alpha_{filter_}']
    filter_tbl = all_extern_tbl[all_extern_tbl['filter'] == filter_]
    filter_tbl.sort('obsdate')
    flux_model = fireball_filter(out.params, filter_tbl['obsdate'], filter_)
    mag_model = flux_to_mag(flux_model)
    
    plt.text(59531,10,'%.4f'%exptime)
    plt.scatter(filter_tbl['obsdate'], filter_tbl['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_])
    plt.plot(filter_tbl['obsdate'], mag_model + offset_key[filter_], c = color_key[filter_], label = f'[{label_key[filter_]}]alpha = {round(alpha.value,2)}')
plt.xlim(59528, 59540)
plt.ylim(24, 8)
plt.legend(loc = 4)
# %%
