#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 19:58:41 2022
\
@author: hhchoi1022
"""
#%%

#%%
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
#import sys
#sys.path.append('/data7/yunyi/temp_supernova/Gitrepo/')
from Research.analysis.observedphot import ObservedPhot
from lmfit import Parameters, minimize
from scipy.interpolate import UnivariateSpline
from Research.model import DOMInteractionL17
from Research.helper import Helper
import matplotlib
#%matplotlib inline
#%% Observation
helper = Helper()
DM = 31.18
ZP = 25
filepath_all = '/data1/supernova_rawdata/SN2021aefx/photometry/all_phot_MW_dereddening_Host_dereddening.dat'
model_directory = '/data1/supernova_model/FB_model'
fit_filterset = 'UBVugri'
fit_start_mjd : int = 59531
fit_end_mjd : int = 59538

# Query data
obs_tbl = ascii.read(filepath_all, format = 'fixed_width')

# Exclude some observatories with poor photometry
observed_data = ObservedPhot(obs_tbl)
observed_data.exclude_observatory(['LasCumbres0.4m', 'Swift'])
observed_data.exclude_filter(['Unfilte'])

# Get detected data
detected_tbl = observed_data.get_data_detected()

# Construct table for fitting
fit_idx = [filter_ in fit_filterset for filter_ in detected_tbl['filter']]
fit_tbl = detected_tbl[fit_idx]
fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
fit_tbl.sort('obsdate')

# Add systematic error
systematic_error_cut = 0.03
def adjust_errors(errors, systematic_error):
    return np.sqrt(errors**2 + systematic_error**2)
obs_lowerror_tbl = fit_tbl[fit_tbl['e_mag'] < systematic_error_cut]
adjusted_errors = np.round(adjust_errors(obs_lowerror_tbl['e_mag'], systematic_error= systematic_error_cut),3)
obs_lowerror_tbl['e_mag'] = adjusted_errors
fit_tbl[fit_tbl['e_mag'] < systematic_error_cut] = obs_lowerror_tbl
fit_tbl['e_mag'] = np.round(fit_tbl['e_mag'], 3)

# Add flux 
fit_tbl['flux'] = helper.mag_to_flux(fit_tbl['mag'])
fit_tbl['e_flux'] = fit_tbl['e_mag']*helper.mag_to_flux(fit_tbl['mag'])*2.303/2.5
fit_tbl['absmag'] = (fit_tbl['mag'] - DM).round(3)

# Visualize
plt.figure(dpi = 400, figsize = (6,6))
ax1, ax2 = observed_data.show_lightcurve(day_binsize = 5,
                              scatter_linewidth=0.6, 
                              scatter_size=50, 
                              scatter_alpha = 1,
                              errorbar_linewidth=0.7, 
                              errorbar_capsize= 3, 
                              color_UB = True,
                              color_BV = True, 
                              color_gr = True, 
                              UL = True, 
                              UL_alpha = 0.8,
                              label = True, 
                              label_location=0, 
                              )
ax1.fill_betweenx(y = [ 30, 0], x1 = fit_start_mjd, x2 = fit_end_mjd, color = 'gray', alpha = 0.2)
ax1.set_ylim(22, 6)
plt.xlim(59525, 59540)
#%% Slicing the data for fitting
filter_tbls = fit_tbl.group_by('filter').groups
filter_key = list(set(fit_tbl['filter']))
filter_key.sort()
fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
y_fit = [np.array((fit_table[filter_]['flux'].tolist())) for filter_ in filter_key]
e_y_fit = [np.array((fit_table[filter_]['e_flux'].tolist())) for filter_ in filter_key]

#%% Setting loss function
def fireball_model(time, 
                   amplitude, 
                   exptime, 
                   alpha):
    flux = amplitude * (time - exptime )**alpha
    np.nan_to_num(flux, copy=False, nan=0.0001)
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
for filter_ in fit_filterset:
    fit_params.add(f'amp_{filter_}', value = 200, min = 0, max = 2000)
    fit_params.add(f'alpha_{filter_}', value = 2, min = 1, max = 10)
#%% Fitting
out = minimize(calc_chisq, fit_params, args = (x_fit, y_fit, e_y_fit, filter_key))
result_values = out.params.valuesdict()
#%% Visualization
color_key, offset_key, _, _, label_key = helper.load_filt_keys(filter_key)
plt.figure(dpi = 300, figsize = (4.5, 6.5))
plt.gca().invert_yaxis()
tbl_UL = observed_data.get_data_ul()
tbl_obs = observed_data.get_data_detected()
ax1, ax2 = observed_data.show_lightcurve(day_binsize = 5,
                            scatter_linewidth=0.5, 
                            scatter_size=50, 
                            scatter_alpha = 1,
                            errorbar_linewidth=0, 
                            errorbar_capsize=0, 
                            color_UB = True,
                            color_BV = True, 
                            color_gr = True, 
                            UL = True, 
                            UL_alpha = 0.8,
                            label = True, 
                            label_location=4, 
                            )
phase_min_FB = np.max([59526, result_values['exptime']])
phase_range_FB = np.arange(phase_min_FB, 59538, 0.1)

for filter_ in filter_key:
    exptime = out.params['exptime']
    amp = out.params[f'amp_{filter_}']
    alpha= out.params[f'alpha_{filter_}']
    flux_model = fireball_filter(out.params, phase_range_FB, filter_)
    mag_model = helper.flux_to_mag(flux_model)
    ax1.plot(phase_range_FB, mag_model + offset_key[filter_], c = color_key[filter_], label = f'[{label_key[filter_]}]alpha = {round(alpha.value,2)}')
    # For color plot 
    if filter_ == 'U':
        mag_U_model = mag_model
    if filter_ == 'B':
        mag_B_model = mag_model
    if filter_ == 'V':
        mag_V_model = mag_model
    if filter_ == 'g':
        mag_g_model = mag_model
    if filter_ == 'r':
        mag_r_model = mag_model

#ax2.plot(phase_range_FB, mag_U_model - mag_B_model -0.5, c = 'cyan', label = 'U-B', linestyle= '-', linewidth = 1, alpha = 1)
ax2.plot(phase_range_FB, mag_B_model - mag_V_model + 0.5, c = 'b', label = 'B-V', linestyle= '-', linewidth = 1, alpha = 1)
ax2.plot(phase_range_FB, mag_g_model - mag_r_model, c = 'g', label = 'g-r', linestyle= '-', linewidth = 1, alpha = 1)
ax1.set_xlim(59528, 59537)
ax2.set_xlim(59528, 59537)
ax1.set_ylim(22.5, 8)
#%%
ax1.clear()
show_idx = [0,3,7, 9, 10, 12]
obs_spec = ascii.read('/data1/supernova_rawdata/SN2021aefx/photometry/all_spec_MW_dereddening_Host_dereddening.dat', format = 'fixed_width')
obs_spec_phot = ObservedPhot(data_tbl  = obs_spec)
obs_spec_phot.data.sort('obsdate')
filt_spec_tbl = obs_spec_phot.get_filt_data(obs_spec_phot.data)
UB_tbl = helper.match_table(filt_spec_tbl['U'], filt_spec_tbl['B'], key = 'obsdate')
BV_tbl = helper.match_table(filt_spec_tbl['B'], filt_spec_tbl['V'], key = 'obsdate')
gr_tbl = helper.match_table(filt_spec_tbl['g'], filt_spec_tbl['r'], key = 'obsdate')

import glob
import matplotlib.cm as cm  # Import the colormap
from Research.spectroscopy.spectroscopyfile import SpectroscopyFile
from Research.spectroscopy import Spectrum
from Research.spectroscopy import TimeSeriesSpectrum
from Research.model import CompanionInteractionK10

num_files = len(show_idx)  # Determine the number of files
colormap = cm.get_cmap('cool', num_files)  # Choose a colormap and set the number of colors
bb_temp = [11000, 6000, 8000, 8000, 9000, 10000]
CEI_model = CompanionInteractionK10(rstar = 2.35, m_wd = 1.3, v9 = 0.9).get_LC(td = np.arange(0.1, 10, 0.1))
exptime_CEI = 59528.97627959757
CEI_model['phase'] = CEI_model['phase'] + exptime_CEI
spl_temp,_ = helper.interpolate_spline(list(CEI_model['phase']), list(CEI_model['Temperature_eff']), show = False)

for i, idx in enumerate(show_idx):
    
    file_ = UB_tbl[idx]['filename_1']
    obsdate = UB_tbl[idx]['obsdate_1']
    print(obsdate - 59529.3318)
    specfile = SpectroscopyFile(file_)
    flux = specfile.flux
    wl = specfile.wavelength
    spec = Spectrum(wl, flux, flux_unit='flamb')
        
    # Get a color from the colormap
    color = colormap(i)
    Planck = helper.planck
    val = Planck(temperature=spl_temp(obsdate), wl_AA=wl)
    spec_bb = Spectrum(wl, val['flamb'], flux_unit='flamb')
    
    # If spec.show() supports a color parameter
    spec.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.05, normalize_cenwl=7500, color=color, label = specfile.obsdate, offset = -2*i, axis = ax1, linewidth = 1)
    if i < 0:
        spec_bb.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.05, normalize_cenwl=7500, color='black', offset = -2*i, axis = ax1, linestyle = '--', linewidth = 0.5)
    ax2.scatter(BV_tbl['obsdate_1'][idx], BV_tbl['mag_1'][idx] - BV_tbl['mag_2'][idx] + 0.5, facecolor = 'b', edgecolor = color, marker = '*', s = 100, alpha = 1)
    ax2.scatter(gr_tbl['obsdate_1'][idx], gr_tbl['mag_1'][idx] - gr_tbl['mag_2'][idx], facecolor = 'g', edgecolor = color, marker = '*', s = 100, alpha = 1)

ax1.tick_params(axis='x', which='both', direction='in', top=True)
ax2.tick_params(axis='x', which='both', direction='in', top=True)
ax1.set_ylim(-10, 5.5)
ax1.set_xlim(3000, 10000)
ax1.set_xticks(np.arange(3000, 11000, 1000), np.arange(3000, 11000, 1000))
# Move x-ticks to the top for ax1
ax1.xaxis.tick_top()  # This moves the x-tick labels to the top of ax1
ax1.xaxis.set_label_position('top')  # This moves the x-axis label to the top
ax1.set_ylabel(rf'Normalized flux ($F_\lambda$) + offset')
ax2.set_xlim(59528, 59537)
# %%
