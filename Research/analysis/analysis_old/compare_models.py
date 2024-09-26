#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:23:50 2022

@author: hhchoi1022
"""
#%%
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec

from astropy.io import ascii
from astropy.table import unique
from astropy.table import Table 
from astropy.table import vstack

import pyphot
from pyphot import unit

from HHsupport_analysis import load_filt_keys
from HHsupport_analysis import mag_to_flux, flux_to_mag
from HHsupport_analysis import interpolate_spline
from HHsupport_analysis import read_HESMA
from HHsupport_analysis import read_Polin2019
os.chdir('/Users/hhchoi1022/Gitrepo/analysis')
color_key, offset_key, _, _, label_key = load_filt_keys()
distancemodulus = 31.17
#%%  Load data
homedir = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/'
obsdir = homedir+'observation/lightcurve'

#IMSNG
all_tbl = Table.read(obsdir+'/IMSNG/IMSNG_hostmwextin3.10.dat', format = 'ascii.fixed_width')
filterlist = sorted(list(set(all_tbl['filter'])))
observatorylist = sorted(list(set(all_tbl['observatory'])))
UL_tbl = all_tbl[all_tbl['status'] == 'UL']
all_tbl = all_tbl[all_tbl['status'] == 'detected']
UL_tbl['absmag'] = UL_tbl['UL5_4'] - distancemodulus
all_tbl['absmag'] = all_tbl['mag'] - distancemodulus
UL_tbl_g, UL_tbl_r, UL_tbl_i = UL_tbl[UL_tbl['filter'] == filterlist[0]], UL_tbl[UL_tbl['filter'] == filterlist[2]], UL_tbl[UL_tbl['filter'] == filterlist[1]]
all_g, all_r, all_i = all_tbl[all_tbl['filter'] == filterlist[0]], all_tbl[all_tbl['filter'] == filterlist[2]], all_tbl[all_tbl['filter'] == filterlist[1]]
all_KCT, all_LSGT, all_RASA = all_tbl[all_tbl['observatory'] == observatorylist[0]], all_tbl[all_tbl['observatory'] == observatorylist[1]], all_tbl[all_tbl['observatory'] == observatorylist[2]],
KCT_g, LSGT_g = all_g[all_g['observatory'] == observatorylist[0]], all_g[all_g['observatory'] == observatorylist[1]]
KCT_r, LSGT_r, RASA_r = all_r[all_r['observatory'] == observatorylist[0]], all_r[all_r['observatory'] == observatorylist[1]], all_r[all_r['observatory'] == observatorylist[2]]
KCT_i, LSGT_i = all_i[all_i['observatory'] == observatorylist[0]], all_i[all_i['observatory'] == observatorylist[1]]

#Hosseinzadeh 2022
extern_tbl = ascii.read(obsdir+'/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat', format = 'fixed_width')
extern_tbl['absmag'] = extern_tbl['mag'] - distancemodulus
extern_tbl = extern_tbl[extern_tbl['observatory'] == 'LasCumbres1m']
extern_filterlist = sorted(list(set(extern_tbl['filter'])))
extern_g, extern_r, extern_i = extern_tbl[extern_tbl['filter'] == 'g'], extern_tbl[extern_tbl['filter'] == 'r'], extern_tbl[extern_tbl['filter'] == 'i']
extern_U, extern_B, extern_V = extern_tbl[extern_tbl['filter'] == 'U'], extern_tbl[extern_tbl['filter'] == 'B'], extern_tbl[extern_tbl['filter'] == 'V']

#Ashall 2022
#ashall_tbl = ascii.read(obsdir+'/Ashall2022/Ashall2022_hostmwextin3.10.dat', format = 'fixed_width')
#ashall_tbl['absmag'] = ashall_tbl['mag'] - distancemodulus
#ashall_u, ashall_g, ashall_r, ashall_i = ashall_tbl[ashall_tbl['filter'] == 'u'], ashall_tbl[ashall_tbl['filter'] == 'g'], ashall_tbl[ashall_tbl['filter'] == 'r'], ashall_tbl[ashall_tbl['filter'] == 'i']
#ashall_B, ashall_V = ashall_tbl[ashall_tbl['filter'] == 'B'], ashall_tbl[ashall_tbl['filter'] == 'V']

#Combine
all_extern_tbl = vstack([all_tbl, extern_tbl])

#%%
#%% Observation
# noextin

filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_noextin.dat'
filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_noextin.dat'
filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_noextin.dat'

# extin2.36

filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat'
filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostmwextin3.10.dat'
filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostmwextin3.10.dat'

tbl1 = ascii.read(filepath_1, format = 'fixed_width')
tbl2 = ascii.read(filepath_2, format = 'fixed_width')
tbl3 = ascii.read(filepath_3, format = 'fixed_width')

tbl_obs = vstack([tbl1, tbl2])
observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
observed_data.data.sort('obsdate')
#observed_data.data = observed_data.get_data_detected()
err_tbl_key = observed_data.data['e_mag'] < 0.01
observed_data.data['e_mag'][err_tbl_key] = 0.03
observed_data.exclude_observatory(['LasCumbres0.4m', 'Swift'])
show_tbl = observed_data.data
#%% Lightcurve
from HHsupport_analysis import load_marker_keys

def plot_LC(exptime, delmag, show_mjd_range = (-5,50), model_tbl = None, days = None, filter_key = 'UBVgri', stretch = None, model = False, title = 'Empty title', text = None):
    import matplotlib.patches as mpatches
    import matplotlib.ticker as ticker

    color_key, offset_key, _, _, label_row_dict = load_filt_keys(filter_key)
    marker= load_marker_keys(4)
    marker = ['.','P','^','s']
    label_row = list(label_row_dict.values())
    label_column = ['H22','KCT','LSGT','RASA36']
    mksize = 30
    scatterlinewidth = 0.7
    
    plt.subplots(figsize = (6,7), dpi = 1000)
    plt.title(os.path.basename(title))
    #plt.gca().invert_yaxis()
    gs = gridspec.GridSpec(nrows = 2, ncols =1, height_ratios= [10, 3], width_ratios = [6])
    ax1 = plt.subplot(gs[0])
    
    ### OBSERVATION ###
    #UL
    ax1.set_title(os.path.basename(title))
    ax1.scatter(UL_tbl_g['obsdate']-exptime, UL_tbl_g['UL5_4']+offset_key['g'], marker = 'v', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =0.5, alpha = 0.3)
    ax1.scatter(UL_tbl_r['obsdate']-exptime, UL_tbl_r['UL5_4']+offset_key['r'], marker = 'v', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =0.5, alpha = 0.3)
    ax1.scatter(UL_tbl_i['obsdate']-exptime, UL_tbl_i['UL5_4']+offset_key['i'], marker = 'v', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =0.5, alpha = 0.3)
    #KCT
    ax1.scatter(KCT_g['obsdate']-exptime, KCT_g['mag']+offset_key['g'], marker = 'P', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =scatterlinewidth)
    ax1.scatter(KCT_r['obsdate']-exptime, KCT_r['mag']+offset_key['r'], marker = 'P', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.scatter(KCT_i['obsdate']-exptime, KCT_i['mag']+offset_key['i'], marker = 'P', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =scatterlinewidth)
    ax1.errorbar(KCT_g['obsdate']-exptime, KCT_g['mag']+offset_key['g'], KCT_g['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['g'])
    ax1.errorbar(KCT_r['obsdate']-exptime, KCT_r['mag']+offset_key['r'], KCT_r['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['r'])
    ax1.errorbar(KCT_i['obsdate']-exptime, KCT_i['mag']+offset_key['i'], KCT_i['e_mag'], fmt = 'none', elinewidth = 0.4, c = color_key['i'])
    #LSGT
    ax1.scatter(LSGT_g['obsdate']-exptime, LSGT_g['mag']+offset_key['g'], marker = '^', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =scatterlinewidth)
    ax1.scatter(LSGT_r['obsdate']-exptime, LSGT_r['mag']+offset_key['r'], marker = '^', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.scatter(LSGT_i['obsdate']-exptime, LSGT_i['mag']+offset_key['i'], marker = '^', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =scatterlinewidth)
    ax1.errorbar(LSGT_g['obsdate']-exptime, LSGT_g['mag']+offset_key['g'], LSGT_g['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['g'])
    ax1.errorbar(LSGT_r['obsdate']-exptime, LSGT_r['mag']+offset_key['r'], LSGT_r['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['r'])
    ax1.errorbar(LSGT_i['obsdate']-exptime, LSGT_i['mag']+offset_key['i'], LSGT_i['e_mag'], fmt = 'none', elinewidth = 0.4, c = color_key['i'])
    #RASA36
    ax1.scatter(RASA_r['obsdate']-exptime, RASA_r['mag']+offset_key['r'], marker = 's', s=mksize,facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.errorbar(RASA_r['obsdate']-exptime, RASA_r['mag']+offset_key['r'], RASA_r['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['r'])
    #Hosseinzadeh2022
    ax1.scatter(extern_g['obsdate']-exptime, extern_g['mag']+offset_key['g'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_r['obsdate']-exptime, extern_r['mag']+offset_key['r'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_i['obsdate']-exptime, extern_i['mag']+offset_key['i'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =scatterlinewidth)
    ax1.errorbar(extern_g['obsdate']-exptime, extern_g['mag']+offset_key['g'], extern_g['e_mag'] , fmt = 'none', elinewidth = 0.4,  c = color_key['g'])
    ax1.errorbar(extern_r['obsdate']-exptime, extern_r['mag']+offset_key['r'], extern_r['e_mag'] , fmt = 'none', elinewidth = 0.4,  c = color_key['r'])
    ax1.errorbar(extern_i['obsdate']-exptime, extern_i['mag']+offset_key['i'], extern_i['e_mag'], fmt = 'none', elinewidth = 0.4, c = color_key['i'])
    ax1.scatter(extern_U['obsdate']-exptime, extern_U['mag']+offset_key['U'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['U'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_B['obsdate']-exptime, extern_B['mag']+offset_key['B'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['B'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_V['obsdate']-exptime, extern_V['mag']+offset_key['V'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['V'],linewidth  =scatterlinewidth)
    ax1.errorbar(extern_U['obsdate']-exptime, extern_U['mag']+offset_key['U'], extern_U['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['U'])
    ax1.errorbar(extern_B['obsdate']-exptime, extern_B['mag']+offset_key['B'], extern_B['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['B'])
    ax1.errorbar(extern_V['obsdate']-exptime, extern_V['mag']+offset_key['V'], extern_V['e_mag'], fmt = 'none', elinewidth = 0.4, c = color_key['V'])
    #Ashall2022
    #ax1.scatter(ashall_u['obsdate']-exptime, ashall_u['mag']+offset_key['u'], marker = 'D', s=mksize, facecolors = 'none' ,edgecolors = color_key['u'],linewidth  =0.5, label = 'u-1')
    #ax1.scatter(ashall_g['obsdate']-exptime, ashall_g['mag']+offset_key['g'], marker = 'D', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =0.5, label = 'g')
    #ax1.scatter(ashall_r['obsdate']-exptime, ashall_r['mag']+offset_key['r'], marker = 'D', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =0.5, label = 'r+1')
    #ax1.scatter(ashall_i['obsdate']-exptime, ashall_i['mag']+offset_key['i'], marker = 'D', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =0.5, label = 'i+2')
    #ax1.errorbar(ashall_u['obsdate']-exptime, ashall_u['mag']+offset_key['u'], ashall_u['e_mag'] , fmt = 'none', capsize = 1, c = color_key['u'])
    #ax1.errorbar(ashall_g['obsdate']-exptime, ashall_g['mag']+offset_key['g'], ashall_g['e_mag'] , fmt = 'none', capsize = 1, c = color_key['g'])
    #ax1.errorbar(ashall_r['obsdate']-exptime, ashall_r['mag']+offset_key['r'], ashall_r['e_mag'] , fmt = 'none', capsize = 1, c = color_key['r'])
    #ax1.errorbar(ashall_i['obsdate']-exptime, ashall_i['mag']+offset_key['i'], ashall_i['e_mag'], fmt = 'none', capsize = 1, c = color_key['i'])
    #ax1.scatter(ashall_B['obsdate']-exptime, ashall_B['mag']+offset_key['B'], marker = 'D', s=mksize, facecolors = 'none' ,edgecolors = color_key['B'],linewidth  =0.5, label = 'B-0.5')
    #ax1.scatter(ashall_V['obsdate']-exptime, ashall_V['mag']+offset_key['V'], marker = 'D', s=mksize, facecolors = 'none' ,edgecolors = color_key['V'],linewidth  =0.5, label = 'V+0.5')
    #ax1.errorbar(ashall_B['obsdate']-exptime, ashall_B['mag']+offset_key['B'], ashall_B['e_mag'] , fmt = 'none', capsize = 1, c = color_key['B'])
    #ax1.errorbar(ashall_V['obsdate']-exptime, ashall_V['mag']+offset_key['V'], ashall_V['e_mag'], fmt = 'none', capsize = 1, c = color_key['V'])
    #Label for observation
    rows = [mpatches.Patch(color=color_key[clr]) for clr in color_key]
    columns = [plt.plot([],[], marker[i], markerfacecolor ='w', markeredgecolor='k')[0] for i in range(len(marker))]
    ax1.legend(rows + columns, label_row + label_column, loc=4, ncol = 2)
    ### OBSERVATION END ###
    
    ### MODEL ###
    if model:
        if not stretch == None:
            for i in range(len(filter_key)):
                ax1.plot(stretch*np.array(days), np.array(list(model_tbl.loc[filter_key[i]])[1:])+offset_key[filter_key[i]]+delmag+distancemodulus, linewidth = 1, c = color_key[filter_key[i]])
        else:
            for i in range(len(filter_key)):
                ax1.plot(days, np.array(list(model_tbl.loc[filter_key[i]])[1:])+offset_key[filter_key[i]]+delmag+distancemodulus, linewidth = 1, c = color_key[filter_key[i]])
    
    if not text == None:
        for i, txt in enumerate(text):
            ax1.text(8,20+0.9*i,s = txt)
    
    
    ax1.set_xlabel('days since explosion')
    ax1.set_ylabel('Apparent magnitude')
    ax1.set_xlim(show_mjd_range[0],show_mjd_range[1])
    ax1.set_yscale('linear')
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax1.set_ylim(24,5)
    #ax2 = ax1.twinx()
    ##ax2.set_ylabel('Absolute magnitude')
    #ax2.set_ylim(24-distancemodulus,5-distancemodulus)


#%%%%%%%%%%%%%%%%%%%%%%%%%% One-time codes
#%% Check data (color)
'''
Hosse_U = binning_table(extern_U, 'obsdate')
Hosse_B = binning_table(extern_B, 'obsdate')
Hosse_V = binning_table(extern_V, 'obsdate')
color_tbl_BV = match_table(Hosse_B, Hosse_V, 'obsdate')
color_tbl_UB = match_table(Hosse_U, Hosse_B, 'obsdate')
#color_tbl_BV_e = color_tbl_BV
#color_tbl_UB_e = color_tbl_UB

plt.figure(dpi = 300)
plt.grid(axis = 'y')
plt.scatter(color_tbl_BV_e['obsdate_1']-59529.19, color_tbl_BV_e['mag_1']-color_tbl_BV_e['mag_2'], c= 'b', label ='[corrected]B-V')
plt.scatter(color_tbl_UB_e['obsdate_1']-59529.19, color_tbl_UB_e['mag_1']-color_tbl_UB_e['mag_2'], c ='cyan', label ='[corrected]U-B')
plt.scatter(color_tbl_BV['obsdate_1']-59529.19, color_tbl_BV['mag_1']-color_tbl_BV['mag_2'], marker = 'x', c= 'b', label ='[raw]B-V')
plt.scatter(color_tbl_UB['obsdate_1']-59529.19, color_tbl_UB['mag_1']-color_tbl_UB['mag_2'], marker = 'x', c ='cyan', label ='[raw]U-B')
plt.legend()
plt.xlim(4.5e-1, 6.3e1)
plt.ylabel('Color')
plt.xlabel('MJD-59529')
plt.xscale('log')
'''
#%% two trends in i band?
#%%%%%%%%%%%%%%%%%%%%%%%%%% One-time codes END
'''
from support import load_marker_keys
markers = load_marker_keys()
obslist = set(extern_g['observatory'])
color_key, offset_key, _, _, label_row_dict = load_filt_keys('all')
plt.figure(dpi = 300)
plt.title('g band near maximum')
plt.grid(axis ='y')
for i, obs in enumerate(obslist):
    obs_tbl_g = extern_g[extern_g['observatory'] == obs]
    plt.xlim(59535, 59565)
    plt.ylim(13, 11)
    clr = 'b'
    if obs == 'LasCumbres1m':
        clr = 'r'
    plt.scatter(obs_tbl_g['obsdate'],obs_tbl_g['mag'], label ='[H22]'+obs, marker =markers[i+3], facecolor = 'none', edgecolor = clr )
plt.scatter(ashall_g['obsdate'],ashall_g['mag'], label ='A22', marker =markers[i+4], facecolor = 'none', edgecolor = 'r' )
plt.scatter(all_g['obsdate'],all_g['mag'], label ='IMSNG', marker =markers[i+1], facecolor = 'none', edgecolor = 'r' )
plt.legend()
'''

#%%%%%%%%%%%%%%%%%%%%%%%%%% One-time codes END
#%%
def synth_phot_HESMA(filename, filter_key ='UBVRIugriz'):
    wl, days, f_lamb, f_nu, mag, _ = read_HESMA(filename)
    lib = pyphot.get_library()
    _, _, _, pyphot_key, _ = load_filt_keys(filter_key)
    result_phot = Table()
    result_phot['filter'] = list(filter_key)
    result_phot.add_index('filter')
    for day in days:
        magset = []
        for filt_ in filter_key:
            filt_pyphot = lib[pyphot_key[filt_]]
            flux = filt_pyphot.get_flux(wl*unit['AA'],f_lamb[str(day)]*unit['ergs/s/cm**2/AA'], axis = 1)
            mag = -2.5*np.log10(flux.value) - filt_pyphot.AB_zero_mag
            magset.append(mag)
        result_phot[f'{day}'] = magset
    return result_phot, days

#%%
def plot_spectrum_HESMA(filename):
    wl, days, f_lamb, f_nu, mag = read_HESMA(filename)
    
    # daily spectrum
    
    for day in days:
        plt.figure(dpi = 100)
        plt.title(os.path.basename(filename))
        plt.plot(wl, f_lamb[f'{day}'], linewidth = 0.5, label = f'Days = {day}')
        plt.legend(loc = 1)
        plt.ylabel(r'flux [ergs/s/cm$**2$/AA]')
        plt.xlabel('wavelength[AA]')
        plt.show()
    
    # spectrum variation
    norm = matplotlib.colors.Normalize(
                    vmin = np.min(days),
                    vmax = np.max(days))
    c_m = matplotlib.cm.afmhot
    s_m = matplotlib.cm.ScalarMappable(cmap = c_m, norm = norm)
    s_m.set_array([])
    
    plt.figure(dpi = 300)
    for i, day in enumerate(days):
        plt.plot(wl, f_lamb[f'{day}'], linewidth = 0.5, color = s_m.to_rgba(day))
    cbar = plt.colorbar(s_m)
    cbar.set_label('Phase')
    
    plt.ylabel(r'flux [ergs/s/cm$**2$/AA]')
    plt.xlabel('wavelength[AA]')
    plt.show()

#%% Filtercurve
def show_filtercurve(filter_key = 'BVRI'):
    lib = pyphot.get_library() # Filter library 
    _, _, _, pyphot_key, _ = load_filt_keys(filter_key)
    show_filter_key = list(pyphot_key.values())
    filterlist_ = lib[show_filter_key]
    plt.figure(dpi = 300)
    for i, filter_ in enumerate(filterlist_):
        if filter_.name =='GROUND_JOHNSON_B':
            plt.plot(np.array(filter_.wavelength), filter_.transmit/100, label = f'{show_filter_key[i]}')
        else:
            plt.plot(np.array(filter_.wavelength), filter_.transmit, label = f'{show_filter_key[i]}')
    plt.legend(loc = 4)
    plt.xlabel('wavelength[AA]')
    plt.ylabel('transmission')
    plt.xlim(3000,9000)
    plt.show()

#%%
def chisq_likelihood(theta, x, y, yerr, x_model, y_model, show = False):
    explosion, deltamag = theta
    func_model = interpolate_spline(x_model + explosion, y_model)
    y_model = func_model(x) + deltamag
    chisq = -0.5 * np.sum(((y - y_model)/(yerr))**2)
    if show:
        xgrid = np.arange(np.min(x),np.max(x), 0.005)
        plt.figure(dpi =150)
        plt.gca().invert_yaxis()
        plt.scatter(x,y,marker= '+', s = 5, c = 'r', label  ='Observed data')
        plt.errorbar(x,y,yerr, fmt= 'none', c = 'r', elinewidth = 4)
        plt.plot(xgrid, func_model(xgrid)+deltamag, c = 'k', linewidth = 1, label = 'Best-matched func $\chi_{reduced}$ = %.2f'%(-chisq/len(x-2)))
        plt.legend(loc = 0)
    return chisq

#%%
def chisq_multi_likelihood(theta, filterlist, obs_tbl, model_tbl, title = None, show = False):
    color_key, offset_key, _, _, label_row = load_filt_keys(filterlist)
    explosion, deltamag, stretch_param = theta
    sum_chisq = 0
    exptime_exp = 59529
    if show:
        plt.figure(figsize = (8,7))
        plt.title(title)
        plt.ylim(-22,-11.0)
        plt.xlim(exptime_exp, exptime_exp+60)
    for filt_ in filterlist:
        obs = obs_tbl[obs_tbl['filter'] == filt_]
        x,y,yerr =  obs['obsdate'],obs['absmag'],obs['e_mag']
        
        idx = np.where(x < exptime_exp+45)
        x_obs = x[idx]
        y_obs = y[idx]
        yerr_obs = yerr[idx]
        x_model = stretch_param * np.array(list(model_tbl['phase']))
        #x_model = stretch_param*(np.array(list(model_tbl.columns)[1:]).astype(float))
        y_model = model_tbl[filt_]
        #y_model = np.array(list(model_tbl.loc[filt_])[1:]).astype(float)
        func_model, _ = interpolate_spline(x_model + explosion, y_model)
        y_model = func_model(x_obs) + deltamag
        chisq = np.sum(((y_obs - y_model)/(yerr_obs))**2)
        reduced_chisq = -chisq / (len(x_obs)-2)
    
        if show:
            #xgrid = np.arange(np.min(x_model)+explosion,np.max(x_model)+explosion, 0.005)
            xgrid = np.arange(np.min(x_obs)-3,np.max(x_obs)+3, 0.005)
            plt.gca().invert_yaxis()
            plt.ylabel('Absolute magnitude')
            plt.scatter(x_obs,y_obs+offset_key[filt_],marker= '+', s = 5, c = color_key[filt_], label  =f'[Observation] {label_row[filt_]}')
            plt.errorbar(x_obs,y_obs+offset_key[filt_],yerr_obs, fmt= 'none', c = color_key[filt_], elinewidth = 4)
            plt.plot(xgrid, func_model(xgrid)+deltamag+offset_key[filt_], c = color_key[filt_], linewidth = 1, label = '[Model] $\chi_{reduced}$ = %.2f'%-reduced_chisq)
            plt.legend(loc = 0)
        sum_chisq += reduced_chisq
    sum_chisq = sum_chisq/(len(filterlist))
    plt.show()
    return sum_chisq

def log_prior(theta):
    explosion, deltamag = theta
    if 59525 < explosion < 59535 and -1 < deltamag < 1:
        return 0.0
    return -np.inf

def log_probability(theta, filterlist, obs_tbl, model_tbl, title = None, show = False):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + chisq_multi_likelihood(theta, filterlist, obs_tbl, model_tbl, title = title, show = False)


#%% Fitting
def fit_obs2model(initial, filterlist, obs_tbl, model_tbl, bnds = None, title = None, show = True):
    from scipy.optimize import minimize
    nll = lambda *args: -chisq_multi_likelihood(*args, filterlist = filterlist, obs_tbl = obs_tbl, model_tbl = model_tbl, title = title, show = show)
    if bnds != None:
        soln = minimize(nll, initial, tol = 0.01, options = dict(maxiter=10000), bounds = bnds)
    else:
        soln = minimize(nll, initial, tol = 0.01, options = dict(maxiter=10000)) 
    return soln

#%% MCMC
def sampling_MCMC(filterlist, soln, obs_tbl, model_tbl):
    import emcee
    import corner
    pos = soln.x + [1e-1,1e-1] * np.random.randn(32, 2)
    nwalkers, ndim = pos.shape
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability, args = (filterlist, obs_tbl, model_tbl)
    )
    sampler.run_mcmc(pos, 500, progress=True)
    samples = sampler.get_chain()
    labels = ["exptime", "delmag"]
    flat_samples = sampler.get_chain(discard=100, thin=1, flat=True)
    fig = corner.corner(
        flat_samples, labels=labels, truths = soln.x
        )

#%% Model(Polin 2019) ####################################################################################################
polin_key = homedir+'model/lightcurve/Polin/ddet_Polin2019/*/*.mag'
filelist = sorted(glob.glob(polin_key)) 
initial = 59529, 0, 1 # explosion time, deltamag, stretch
bnds = ((59525, 59529.3315), (-3,3), (0.5,3))
bnds = ((59525, 59529.3315), (-3,3), (0.999,1.001))
fit_start_mjd = 59529
fit_end_mjd = 59539
all_extern_tbl = all_extern_tbl[(all_extern_tbl['obsdate']>fit_start_mjd)&(all_extern_tbl['obsdate']<fit_end_mjd)]

fit_filterlist = 'UBgVri'
#%%
import matplotlib.pyplot as plt
result_Polin = Table()
tbl_files = []
tbl_exptimes = []
tbl_delmags = []
tbl_fun = []
tbl_success = []
tbl_stretches =[]
for file in filelist:
    result_phot = read_Polin2019(file)
    #result_phot.add_index('filter')
    soln = fit_obs2model(initial, fit_filterlist, all_extern_tbl, result_phot, bnds = bnds, title = os.path.basename(file), show= False)
    tbl_exptimes.append(soln.x[0])
    tbl_delmags.append(soln.x[1])
    tbl_files.append(file)
    tbl_fun.append(soln.fun)
    tbl_success.append(soln.success)
    tbl_stretches.append(soln.x[2])
result_Polin['status'] = tbl_success
result_Polin['fun'] = tbl_fun
result_Polin['exptime'] = tbl_exptimes
result_Polin['delmags'] = tbl_delmags
result_Polin['file'] = tbl_files
result_Polin['stretch'] = tbl_stretches
result_Polin.sort('fun')
result_Polin.write('Fit_Polin2019', format = 'ascii.fixed_width', overwrite = True)
#%% Good results visualization
if not 'result_Polin' in globals():
    result_Polin = ascii.read('Fit_Polin2019', format = 'fixed_width')
#%%
from observedphot import ObservedPhot
initial = 59529, 0, 1 # explosion time, deltamag, stretch
fit_filterlist = 'UBgVri'
bnds = ((59525, 59529.3315), (-3,3), (0.999,1.001))
obs = ObservedPhot(show_tbl)
plt.figure(figsize = (8,5), dpi = 300)
for file in result_Polin[0:1]['file']:
    
    result_phot = read_Polin2019(file)
    phase = result_phot['phase']
    soln = fit_obs2model(initial, fit_filterlist, all_extern_tbl, result_phot, bnds = bnds, title = os.path.basename(file), show= False)
    exptime, deltamag, stretch_1st = soln.x
    obs.show_lightcurve(color_BV= False, color_gr= False, label = True, label_location= 4)
    for filter_ in fit_filterlist:
        plt.plot(result_phot['phase'] + exptime, result_phot[filter_] + deltamag +  distancemodulus + offset_key[filter_], c = color_key[filter_], label = label_key[filter_])
    plt.title(os.path.basename(file))
    plt.xlim(59528, 59540)
    plt.ylim(25, 5)

#%% Model(HSEMA) ####################################################################################################
initial = 59529, 0, 1 # explosion time, deltamag, stretch
bnds = ((59528, 59529.3317), (-2,2), (0.999, 1.001))
fit_filterlist = 'BVgri'

from hesmaphot import HesmaPhot
hesma = HesmaPhot()
modellist = hesma._authorlist
#%%

result_HESMA = Table()
tbl_model = []
tbl_spec = []
tbl_exptimes = []
tbl_delmags = []
tbl_fun = []
tbl_success = []
#tbl_stretches = []
tbl_early = []
for model in modellist:
    speckey = homedir+f'model/spectrum/HESMA/{model}/*spectra.dat'
    speclist = sorted(glob.glob(speckey))

    for spec in speclist:
        try:
            model_key = os.path.basename(spec).split('spectra')[0]
            result_phot = hesma.synthphot(spec, early = False)
            result_phot = result_phot[result_phot['phase'] > 0.75]
            soln = fit_obs2model(initial, fit_filterlist, all_extern_tbl, result_phot, bnds = bnds, title = os.path.basename(spec), show= False)
            tbl_exptimes.append(soln.x[0])
            tbl_delmags.append(soln.x[1])
            tbl_model.append(os.path.basename(model))
            tbl_spec.append(spec)
            tbl_fun.append(soln.fun)
            tbl_success.append(soln.success)
            earlykey = homedir+f'model/lightcurve/HESMA/{model}/*{model_key}*early*.dat'
            if len(glob.glob(earlykey)) > 0:
                tbl_early.append(True)
            else:
                tbl_early.append(False)
        except:
            pass
#%%
result_HESMA = Table()
result_HESMA['status'] = tbl_success
result_HESMA['model'] = tbl_model
result_HESMA['spec'] = tbl_spec
result_HESMA['fun'] = tbl_fun
result_HESMA['exptime'] = tbl_exptimes
result_HESMA['delmag'] = tbl_delmags
result_HESMA['earlyLC'] = tbl_early 
result_HESMA.sort('fun')
result_HESMA.write('Fit_HESMA', format = 'ascii.fixed_width', overwrite= True)
#%% Good results visualization
result_HESMA = ascii.read('Fit_HESMA', format = 'fixed_width')
result_HESMA.sort('fun')
#%%
from hesmaphot import HesmaPhot
from observedphot import ObservedPhot
hesma = HesmaPhot()
obs = ObservedPhot(all_extern_tbl)
#%%
for file in result_HESMA[:2]['spec']:
    result_phot = hesma.synthphot(file)
    #plt.title()
    phase = result_phot['phase']
    soln = fit_obs2model(initial, fit_filterlist, all_extern_tbl, result_phot, bnds = bnds, title = os.path.basename(file), show= False)
    exptime, deltamag, stretch_1st = soln.x
    
    plt.figure(figsize = (4,6), dpi = 300)
    plt.title(r'Classic detonation [1.06$M_\odot$ WD] (Sim + 2010)')
    
    obs.show_lightcurve(color_BV= False, color_gr= False, label= True, label_location= 0)
    for filter_ in 'UBVgri':
        plt.plot(result_phot['phase'] + exptime, result_phot[filter_] + deltamag +  distancemodulus + offset_key[filter_], c = color_key[filter_], label = label_key[filter_])
    plt.ylabel('Apparent magnitude[AB]')
    plt.xlim(59525, 59580)
    plt.ylim(25, 5)
    

#%% Early light curve visualization ################################################################################
best_HESMA = result_HESMA[0]
modellist = os.listdir(homedir+'model/lightcurve/HESMA')
modellist = [model for model in modellist if not model.startswith('.DS_')]
model = os.path.basename(os.path.dirname(best_HESMA['spec']))
early_lckey = homedir+f'model/lightcurve/HESMA/{model}/*lightcurves_early.dat'
early_lclist = glob.glob(early_lckey)
early_lc_file = early_lclist[0]
early_lc = ascii.read(early_lc_file, format = 'fixed_width')
#%%
#early_lc_file
# Spline interpolation
spline_B,_ = interpolate_spline(extern_B['obsdate'], extern_B['absmag'], smooth = 0.1)
spline_V,_ = interpolate_spline(extern_V['obsdate'], extern_V['absmag'], smooth = 0.1)
spline_U,_ = interpolate_spline(extern_U['obsdate'], extern_U['absmag'], smooth = 0.3)
all_r.sort('obsdate')
#spline_r,_ = interpolate_spline(all_r['obsdate'], all_r['absmag'], smooth = 0.6)

#HESMA

initial = 59529, 0, 1 # explosion time, deltamag, stretch
bnds = ((59528, 59529.3317), (-2,2), (0.999, 1.001))
fit_filterlist = 'BgVri'
#early_lc = read_HESMA(early_lc_file)
early_lc_file = best_HESMA['spec']
result_phot = hesma.synthphot(best_HESMA['spec'])


soln = fit_obs2model(initial, fit_filterlist, all_extern_tbl, result_phot, bnds = bnds, title = os.path.basename(best_HESMA['spec']), show= False)
explosion_1st, deltamag_1st, stretch_1st = soln.x

#####
show_filt_set = ['U','B','V']
filt_set = list(early_lc.columns[1:])
early_lc['U-B'] = early_lc['U']- early_lc['B']
early_lc['B-V'] = early_lc['B']- early_lc['V']
color_key, offset_key, _, _, label_row  = load_filt_keys(filt_set)

plt.figure(figsize = (5,8), dpi = 300)
#plt.gca().invert_yaxis()
plt.subplots_adjust(hspace=0)
#plt.xticks(visible=False)
gs = gridspec.GridSpec(nrows = 2, ncols =1, height_ratios= [10, 3], width_ratios = [6])
#ax0 = plt.subplot()
ax0 = plt.subplot(gs[0])
ax0.set_title('Classic detonation [1.06M$_\odot$ WD] (Sim10)')#f'{os.path.basename(early_lc_file)}')

for filt_ in show_filt_set:
    ax0.plot(early_lc['t'], early_lc[f'{filt_}']+offset_key[f'{filt_}']+ deltamag_1st, c = color_key[f'{filt_}'], linewidth = 1)
ax0.scatter(extern_U['obsdate']-explosion_1st, extern_U['absmag']+offset_key['U'], marker = '.', c = color_key['U'], label = label_row['U'])
ax0.scatter(extern_B['obsdate']-explosion_1st, extern_B['absmag']+offset_key['B'], marker = '.', c = color_key['B'], label = label_row['B'])
ax0.scatter(extern_V['obsdate']-explosion_1st, extern_V['absmag']+offset_key['V'], marker = '.', c = color_key['V'], label = label_row['V'])
ax0.errorbar(extern_U['obsdate']-explosion_1st, extern_U['absmag']+offset_key['U'], extern_U['e_mag'] , fmt = 'none', capsize = 1, c = color_key['U'])
ax0.errorbar(extern_B['obsdate']-explosion_1st, extern_B['absmag']+offset_key['B'], extern_B['e_mag'] , fmt = 'none', capsize = 1, c = color_key['B'])
ax0.errorbar(extern_V['obsdate']-explosion_1st, extern_V['absmag']+offset_key['V'], extern_V['e_mag'], fmt = 'none', capsize = 1, c = color_key['V'])
#ax0.scatter(all_r['obsdate']-explosion_1st, all_r['absmag']+offset_key['R'], marker = 'o', facecolor = 'none', edgecolor = color_key['R'], label = 'r+2')
ax0.plot(extern_B['obsdate']-explosion_1st, spline_B(extern_B['obsdate'])+offset_key['B'], linestyle = '--', linewidth = 0.5, c = color_key['B'])
ax0.plot(extern_V['obsdate']-explosion_1st, spline_V(extern_V['obsdate'])+offset_key['V'], linestyle = '--', linewidth = 0.5, c = color_key['V'])
ax0.plot(extern_U['obsdate']-explosion_1st, spline_U(extern_U['obsdate'])+offset_key['U'], linestyle = '--', linewidth = 0.5, c = color_key['U'])
#ax0.plot(all_r['obsdate']-explosion_1st, spline_r(all_r['obsdate'])+offset_key['R'], linestyle = '--', linewidth = 0.5, c = color_key['R'])

ax0.text(1.5,-8,s = r'$t_{exp} = %.4f$' %explosion_1st)
ax0.text(1.5,-7,s = 'mag_offset = %.2f' %deltamag_1st)
ax0.legend(loc = 4)
ax0.set_xlabel('MJD')
ax0.set_ylabel('Absolute mag')
ax0.set_xlim(-0.2,10.2)
ax0.set_ylim(-6.5,-24.2)
ax1 = plt.subplot(gs[1])
ax1.set_xlim(-0.2, 10.2)
ax1.set_ylim(-1,1.5)
ax1.axhline(0, c = 'k', linewidth=1)
ax1.plot(early_lc['t'][0:], early_lc['U-B'][0:], label = 'U-B[model]', c= 'cyan')
ax1.plot(early_lc['t'][0:], early_lc['B-V'][0:], label = 'B-V[model]', c= 'b')
ax1.scatter(extern_B['obsdate']-explosion_1st, spline_B(extern_B['obsdate'])-spline_V(extern_B['obsdate']), label = 'B-V', marker = 'x', c= 'b')
ax1.scatter(extern_U['obsdate']-explosion_1st, spline_U(extern_U['obsdate'])-spline_B(extern_U['obsdate']), label = 'U-B', marker = 'x', c= 'cyan')
ax1.set_ylabel(r'$Color$')
ax1.set_xlabel('Days from the explosion')
ax1.legend()

#%% Best_matched + companion interaction(K10) ####################################################################################
if not 'result_Polin' in globals():
    result_Polin = ascii.read('Fit_Polin2019', format = 'fixed_width')
def get_chisq_modelcombined(filt_tbl, spl_model, comp_interac_model, explosion, rstar, explosion_interact):
    obs_mjd,obs_mag,obs_magerr =  filt_tbl['obsdate'],filt_tbl['mag'],filt_tbl['e_mag']
    filter_ = filt_tbl['filter'][0]
    obs_flux = mag_to_flux(obs_mag)
    obs_fluxerr = obs_magerr*obs_flux*2.303/2.5
    model_mag = spl_model(obs_mjd-explosion)+distancemodulus
    model_flux = mag_to_flux(model_mag)
    interactionmodel_mag = np.array(comp_interac_model(obs_mjd-explosion_interact, rstar, filter_)) + distancemodulus + deltamag_1st
    interactionmodel_flux = mag_to_flux(np.array(interactionmodel_mag))
    model_flux = model_flux+ interactionmodel_flux
    chisq = np.sum(((obs_flux - model_flux)/(obs_fluxerr))**2)
    reduced_chisq = chisq/len(obs_flux-2)
    return reduced_chisq

def get_multichisq_modelcombined(theta, obs_tbl, model, comp_interac_model, explosion):
    rstar, explosion_interact = theta[0], theta[1]
    sum_chisq = 0
    filtset = set(obs_tbl['filter'])
    for filt_ in filtset:
        spl_model  = model[filt_]
        filt_tbl = obs_tbl[obs_tbl['filter'] == filt_]
        chisq = get_chisq_modelcombined(filt_tbl, spl_model, comp_interac_model, explosion, rstar, explosion_interact)
        sum_chisq += chisq
    return sum_chisq

explosion_model = explosion_1st 
fit_tbl = vstack([extern_U, extern_B, extern_V])
fit_tbl['filter'][fit_tbl['filter']=='r'] = 'R'
fit_maxcut = explosion_model+10
fit_mincut = explosion_model
show_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut)]
fit_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut) & (fit_tbl['obsdate'] > fit_mincut)]
filtset = 'UBV'#set(fit_tbl['filter'])
rstar = 1
explosion_interact = 59529
initial = tuple([rstar,explosion_interact])
bnds = (0.5, 10),(59525, 59529.0416) # rstar, explosion_interact
xrange = np.arange(0.1,10,0.1)
spl_U,_ = interpolate_spline(early_lc['t'], early_lc['U']+deltamag_1st)
spl_B,_ = interpolate_spline(early_lc['t'], early_lc['B']+deltamag_1st)
spl_V,_ = interpolate_spline(early_lc['t'], early_lc['V']+deltamag_1st)
spl_R,_ = interpolate_spline(early_lc['t'], early_lc['R']+deltamag_1st)
model = dict(U = spl_U,
             B = spl_B,
             V = spl_V,
             R = spl_R)
#%%

from scipy.optimize import minimize
from Companion_interaction_K10_old import Compinteract_K10
nll = lambda *args: get_multichisq_modelcombined(*args, obs_tbl = fit_tbl, model = model, comp_interac_model = Compinteract_K10, explosion = explosion_model)
soln = minimize(nll, x0 = initial, method = "SLSQP", tol = 0.001,options = dict(maxiter=30000), bounds = bnds) # tol = step_size 

 #%% Visualization
import matplotlib.patches as mpatches
rstar, explosion_interact = soln.x[0], soln.x[1]

color_key, offset_key, _, _, label_row = load_filt_keys(filtset)

fig = plt.subplots(figsize = (5,8), dpi = 500)
plt.subplots_adjust(hspace=0)
gs = gridspec.GridSpec(nrows = 2, ncols =1, height_ratios= [10, 3], width_ratios = [6])
ax1 = plt.subplot(gs[0])
ax1.set_xlabel([5])
for filt_ in filtset:
    spl_model  = model[filt_]
    filt_tbl = fit_tbl[fit_tbl['filter'] == filt_]
    filt_show_tbl = show_tbl[show_tbl['filter'] == filt_]
    chisq = get_chisq_modelcombined(filt_tbl, spl_model, Compinteract_K10, explosion_model, rstar, explosion_interact)
    xrange_model = np.arange(0.1, 10, 0.1)
    #xrange_K10 = np.arange(explosion_inter+0.01, fit_maxcut, 0.1)
    flux_model = np.nan_to_num(mag_to_flux(spl_model(xrange_model)+distancemodulus),0)
    flux_K10 = np.nan_to_num(mag_to_flux(np.array(Compinteract_K10(xrange_model, rstar, filt_))+distancemodulus),0)
    flux_both = flux_K10 + flux_model
    # Observation
    ax1.scatter(filt_show_tbl['obsdate']-explosion_model, filt_show_tbl['mag']+offset_key[filt_], facecolor = 'none', edgecolor = color_key[filt_],label = label_row[filt_])
    ax1.errorbar(filt_show_tbl['obsdate']-explosion_model, filt_show_tbl['mag']+offset_key[filt_], filt_show_tbl['e_mag'], fmt ='none', color = color_key[filt_], capsize = 3)
    # Model
    ax1.set_title(r'Pure detonation [1.06M$_\odot$ WD] (Sim + 2010)')
    ax1.plot(xrange_model, flux_to_mag(flux_both)+ offset_key[filt_], linewidth = 1.5,  c=color_key[filt_]) 
    ax1.text(7., spl_model(xrange_model[-1])+offset_key[filt_]-0.3+distancemodulus, s = r'$\bar{\chi}^2$ = %.2f'%chisq)
    if filt_ == 'U':
        obs_U = filt_tbl
        obs_U.sort('obsdate')
        spl_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], smooth = 2)
        days_U, model_U = xrange_model, flux_to_mag(flux_both)
    if filt_ == 'B':
        obs_B = filt_tbl
        obs_B.sort('obsdate')
        spl_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'])
        days_B, model_B = xrange_model, flux_to_mag(flux_both)
    if filt_ == 'V':
        obs_V = filt_tbl
        obs_V.sort('obsdate')
        spl_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'])
        days_V, model_V = xrange_model, flux_to_mag(flux_both)
        
for filt_ in filtset:
    spl_model  = model[filt_]
    flux_model = np.nan_to_num(mag_to_flux(spl_model(xrange_model)+distancemodulus),0)
    ax1.plot(xrange_model, flux_to_mag(flux_model)+offset_key[filt_], linestyle = '--', linewidth = 1,  c=color_key[filt_], label = 'Pure detonation model')
for filt_ in filtset:
    flux_K10 = np.nan_to_num(mag_to_flux(np.array(Compinteract_K10(xrange_model, rstar, filt_))+distancemodulus),0)
    ax1.plot(xrange_model, flux_to_mag(flux_K10) + offset_key[filt_], linestyle = 'dotted', linewidth = 1,  c=color_key[filt_])
ax1.text(6, 17.8, s = r'$R_{companion}$ = %.2f $R_{\odot}$'%(round(rstar,2)))
ax1.text(6, 18.5, s = r'$t_{Ni} = %.4f$' %explosion_model)
ax1.text(6, 19.2, s = r'$t_{Comp}$  = %.4f' %explosion_interact)
ax1.grid(axis = 'y')
ax1.set_ylim(22,6)
ax1.set_xlim(-0.2, 10.2)
#ax1.legend(loc = 4, ncol = 2)
ax1.set_ylabel('apparent magnitude[Vega]')

ax2 = plt.subplot(gs[1], sharex=ax1)
ax2.axhline(0, c = 'k', linewidth=1)
ax2.scatter(obs_B['obsdate']-explosion_model, spl_U(obs_B['obsdate'])-spl_B(obs_B['obsdate']), marker = 'x', label = 'U-B', c = 'cyan')
ax2.scatter(obs_B['obsdate']-explosion_model, spl_B(obs_B['obsdate'])-spl_V(obs_B['obsdate']), marker = 'x', label = 'B-V', c = 'b')
ax2.set_xlim(-0.2,10.2)
ax2.set_ylim(-0.5,1)
ax2.plot(days_U, model_U-model_B, linewidth  =1, label = 'U-B[model]', c = 'cyan')
ax2.plot(days_B, model_B-model_V, linewidth  =1, label = 'B-V[model]', c = 'b')
ax2.set_xlabel('days since explosion')
ax2.set_ylabel('color')
ax2.legend(loc=1)

marker = [':','--','solid']
#label_row = list(label_row_dict.values())
label_column = ['Companion interaction','Pure detonation','Combined']
rows = [mpatches.Patch(color=color_key[clr]) for clr in color_key]
columns = [plt.plot([],[], linestyle = marker[i], c='k')[0] for i in range(len(marker))]
ax1.legend( columns, label_column, loc=4, ncol = 1)
#%%





#%% HESMA model + companion interaction(K10) only for early ####################################################################################

def get_chisq_modelcombined(filt_tbl, spl_model, comp_interac_model, explosion, rstar, explosion_interact, deltamag):
    obs_mjd,obs_mag,obs_magerr =  filt_tbl['obsdate'],filt_tbl['mag'],filt_tbl['e_mag']
    filter_ = filt_tbl['filter'][0]
    obs_flux = mag_to_flux(obs_mag)
    obs_fluxerr = obs_magerr*obs_flux*2.303/2.5
    model_mag = spl_model(obs_mjd-explosion)+distancemodulus
    model_flux = mag_to_flux(model_mag)
    interactionmodel_mag = np.array(comp_interac_model(obs_mjd-explosion_interact, rstar, filter_)) + distancemodulus + deltamag
    interactionmodel_flux = mag_to_flux(np.array(interactionmodel_mag))
    model_flux = model_flux+ interactionmodel_flux
    chisq = np.sum(((obs_flux - model_flux)/(obs_fluxerr))**2)
    reduced_chisq = chisq/len(obs_flux-2)
    return reduced_chisq

def get_multichisq_modelcombined(theta, obs_tbl, model, comp_interac_model):
    rstar, explosion_interact, explosion, deltamag = theta[0], theta[1], theta[2], theta[3]
    sum_chisq = 0
    filtset = set(obs_tbl['filter'])
    for filt_ in filtset:
        spl_model  = model[filt_]
        filt_tbl = obs_tbl[obs_tbl['filter'] == filt_]
        chisq = get_chisq_modelcombined(filt_tbl, spl_model, comp_interac_model, explosion, rstar, explosion_interact, deltamag)
        sum_chisq += chisq
    return sum_chisq

best_HESMA = result_HESMA[1]
modellist = os.listdir(homedir+'model/lightcurve/HESMA')
modellist = [model for model in modellist if not model.startswith('.DS_')]
model = os.path.basename(os.path.dirname(best_HESMA['spec']))
early_lckey = homedir+f'model/lightcurve/HESMA/{model}/*lightcurves_early.dat'
early_lclist = glob.glob(early_lckey)
early_lc_file = early_lclist[0]
early_lc = ascii.read(early_lc_file, format = 'fixed_width')

fit_tbl = vstack([extern_U, extern_B, extern_V])
fit_tbl['filter'][fit_tbl['filter']=='r'] = 'R'
fit_mincut = 59529
fit_maxcut = 59539
show_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut)]
fit_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut) & (fit_tbl['obsdate'] > fit_mincut)]
filtset = 'UBV'#set(fit_tbl['filter'])
rstar = 1
explosion_interact = 59529
explosion_model = 59529
magoffset = 0
initial = tuple([rstar,explosion_interact, explosion_model, magoffset])
bnds = (0.5, 10),(59525, 59529.0416), (59525, 59535), (-2, 2) # rstar, explosion_interact
xrange = np.arange(0.1,10,0.1)
spl_U,_ = interpolate_spline(early_lc['t'], early_lc['U'])
spl_B,_ = interpolate_spline(early_lc['t'], early_lc['B'])
spl_V,_ = interpolate_spline(early_lc['t'], early_lc['V'])
spl_R,_ = interpolate_spline(early_lc['t'], early_lc['R'])
model = dict(U = spl_U,
             B = spl_B,
             V = spl_V,
             R = spl_R)
#%%

from scipy.optimize import minimize
from Companion_interaction_K10_old import Compinteract_K10
nll = lambda *args: get_multichisq_modelcombined(*args, obs_tbl = fit_tbl, model = model, comp_interac_model = Compinteract_K10 )
soln = minimize(nll, x0 = initial, method = "SLSQP", tol = 0.001,options = dict(maxiter=30000), bounds = bnds) # tol = step_size 

 #%% Visualization
import matplotlib.patches as mpatches
rstar, explosion_interact, explosion_model, deltamag = soln.x[0], soln.x[1], soln.x[2], soln.x[3]

color_key, offset_key, _, _, label_row = load_filt_keys(filtset)

fig = plt.subplots(figsize = (5,8), dpi = 500)
plt.subplots_adjust(hspace=0)
gs = gridspec.GridSpec(nrows = 2, ncols =1, height_ratios= [10, 3], width_ratios = [6])
ax1 = plt.subplot(gs[0])
ax1.set_xlabel([5])
for filt_ in filtset:
    spl_model  = model[filt_]
    filt_tbl = fit_tbl[fit_tbl['filter'] == filt_]
    filt_show_tbl = show_tbl[show_tbl['filter'] == filt_]
    chisq = get_chisq_modelcombined(filt_tbl, spl_model, Compinteract_K10, explosion_model, rstar, explosion_interact, deltamag = deltamag)
    xrange_model = np.arange(0.1, 10, 0.1)
    #xrange_K10 = np.arange(explosion_inter+0.01, fit_maxcut, 0.1)
    flux_model = np.nan_to_num(mag_to_flux(spl_model(xrange_model)+distancemodulus + deltamag),0)
    flux_K10 = np.nan_to_num(mag_to_flux(np.array(Compinteract_K10(xrange_model, rstar, filt_))+distancemodulus),0)
    flux_both = flux_K10 + flux_model
    # Observation
    ax1.scatter(filt_show_tbl['obsdate']-explosion_model, filt_show_tbl['mag']+offset_key[filt_], facecolor = 'none', edgecolor = color_key[filt_],label = label_row[filt_], marker ='^')
    ax1.errorbar(filt_show_tbl['obsdate']-explosion_model, filt_show_tbl['mag']+offset_key[filt_], filt_show_tbl['e_mag'], fmt ='none', color = color_key[filt_], capsize = 3)
    # Model
    ax1.set_title(r'Classic detonation [1.06M$_\odot$ WD] (Sim + 2010)')
    ax1.plot(xrange_model, flux_to_mag(flux_both)+ offset_key[filt_], linewidth = 1.5,  c=color_key[filt_]) 
    #ax1.text(7., spl_model(xrange_model[-1])+offset_key[filt_]-0.3+distancemodulus, s = r'$\bar{\chi}^2$ = %.2f'%chisq)
    if filt_ == 'U':
        obs_U = filt_tbl
        obs_U.sort('obsdate')
        spl_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], smooth = 2)
        days_U, model_U = xrange_model, flux_to_mag(flux_both)
    if filt_ == 'B':
        obs_B = filt_tbl
        obs_B.sort('obsdate')
        spl_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'])
        days_B, model_B = xrange_model, flux_to_mag(flux_both)
    if filt_ == 'V':
        obs_V = filt_tbl
        obs_V.sort('obsdate')
        spl_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'])
        days_V, model_V = xrange_model, flux_to_mag(flux_both)
        
for filt_ in filtset:
    spl_model  = model[filt_]
    flux_model = np.nan_to_num(mag_to_flux(spl_model(xrange_model)+distancemodulus),0)
    ax1.plot(xrange_model, flux_to_mag(flux_model)+offset_key[filt_], linestyle = '--', linewidth = 1,  c=color_key[filt_], label = 'Classic detonation model')
for filt_ in filtset:
    flux_K10 = np.nan_to_num(mag_to_flux(np.array(Compinteract_K10(xrange_model, rstar, filt_))+distancemodulus),0)
    ax1.plot(xrange_model, flux_to_mag(flux_K10) + offset_key[filt_], linestyle = 'dotted', linewidth = 1,  c=color_key[filt_])
ax1.text(6, 17.8, s = r'$R_{companion}$ = %.2f $R_{\odot}$'%(round(rstar,2)))
ax1.text(6, 18.5, s = r'$t_{Ni} = %.4f$' %explosion_model)
ax1.text(6, 19.2, s = r'$t_{Comp}$  = %.4f' %explosion_interact)
ax1.grid(axis = 'y')
ax1.set_ylim(22,6)
ax1.set_xlim(-0.2, 10.2)
#ax1.legend(loc = 4, ncol = 2)
ax1.set_ylabel('apparent magnitude[Vega]')

ax2 = plt.subplot(gs[1], sharex=ax1)
ax2.axhline(0, c = 'k', linewidth=1)
ax2.scatter(obs_B['obsdate']-explosion_model, spl_U(obs_B['obsdate'])-spl_B(obs_B['obsdate']), marker = 'x', label = 'U-B', c = 'cyan')
ax2.scatter(obs_B['obsdate']-explosion_model, spl_B(obs_B['obsdate'])-spl_V(obs_B['obsdate']), marker = 'x', label = 'B-V', c = 'b')
ax2.set_xlim(-0.2,10.2)
ax2.set_ylim(-0.5,1)
ax2.plot(days_U, model_U-model_B, linewidth  =1,  c = 'cyan')
ax2.plot(days_B, model_B-model_V, linewidth  =1,  c = 'b')
ax2.set_xlabel('days since explosion')
ax2.set_ylabel('color')
ax2.legend(loc=1)

marker = [':','--','solid']
#label_row = list(label_row_dict.values())
label_column = ['Companion interaction','Classic detonation','Combined']
rows = [mpatches.Patch(color=color_key[clr]) for clr in color_key]
columns = [plt.plot([],[], linestyle = marker[i], c='k')[0] for i in range(len(marker))]
ax1.legend( columns, label_column, loc=4, ncol = 1)





#%% Fitting with Fireball 


# Fireball model with alpha as a free parameter
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux


def get_chisq_fireball(filt_tbl, model, explosion, amplitude, alpha):
    obs_mjd,obs_mag,obs_magerr =  filt_tbl['obsdate'],filt_tbl['mag'],filt_tbl['e_mag']
    obs_flux = mag_to_flux(obs_mag, zp = 25)
    obs_fluxerr = obs_magerr*obs_flux*2.303/2.5
    model_flux = fireball_model(obs_mjd, amplitude, explosion, alpha)
    chisq = np.sum(((obs_flux - model_flux)/(obs_fluxerr))**2)
    reduced_chisq = chisq/len(obs_flux-2)
    return reduced_chisq

def get_multichisq_fireball(theta, obs_tbl, model):
    explosion= theta[0]
    ampl_cut = (len(theta)+1)//2
    amplitudelist, alphalist = theta[1:ampl_cut], theta[ampl_cut:]
    sum_chisq = 0
    filtset = set(obs_tbl['filter'])
    if (len(filtset) != len(amplitudelist)) | (len(filtset) != len(alphalist)):
        raise ValueError('Dimension mismatch');
    for filt_, amplitude, alpha in zip(filtset, amplitudelist, alphalist):
        filt_tbl = obs_tbl[obs_tbl['filter'] == filt_]
        chisq = get_chisq_fireball(filt_tbl, model, explosion, amplitude, alpha)
        sum_chisq += chisq
    average_chisq = sum_chisq/len(filtset)
    return average_chisq


#%% 

fit_tbl = vstack([all_extern_tbl])
fit_maxcut = 59538
fit_mincut = 59530
show_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut)]
fit_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut) & (fit_tbl['obsdate'] > fit_mincut)]
filtset = set(fit_tbl['filter'])
amplitudelist = 1500
alphalist = 2 
explosion = 59529
initial = [explosion]+ [amplitudelist] * len(filtset)+ [alphalist] * len(filtset)
bnd_explosion = 1.5
bnd_amplitude = 1500
bnd_alpha = 1.5
bnd = [bnd_explosion]+ [bnd_amplitude] * len(filtset)+ [bnd_alpha] * len(filtset)

bnds = ()
for value, bound in zip(initial, bnd):
    cut = (value-bound, value+bound)
    bnds += (cut,)

from scipy.optimize import minimize
nll = lambda *args: get_multichisq_fireball(*args, obs_tbl = fit_tbl, model = fireball_model)
soln = minimize(nll, x0 = initial, method = "SLSQP", tol = 0.00001,options = dict(maxiter=30000), bounds = bnds) # tol = step_size 


#%%
exptime = soln.x[0]
ampl_cut = (len(soln.x)+1)//2
amplitudelist, alphalist = soln.x[1:ampl_cut], soln.x[ampl_cut:]
color_key, offset_key, _, _, label_row = load_filt_keys(filtset)
from observedphot import ObservedPhot

filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat'
filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostmwextin3.10.dat'
filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostmwextin3.10.dat'

tbl1 = ascii.read(filepath_1, format = 'fixed_width')
tbl2 = ascii.read(filepath_2, format = 'fixed_width')
tbl3 = ascii.read(filepath_3, format = 'fixed_width')

tbl_obs = vstack([tbl1, tbl2])
det_tbl = tbl_obs[tbl_obs['status'] == 'detected']
ul_tbl = tbl_obs[tbl_obs['status'] == 'UL']
#det_observed = ObservedPhot(tbl_obs)
det_observed = ObservedPhot(det_tbl)
ul_observed = ObservedPhot(ul_tbl)
#%%

import matplotlib.pyplot as plt
plt.figure(dpi = 300, figsize = (3, 5.38))
plt.xlim(59525, 59539)
det_observed.show_lightcurve(label =False, color_BV = False, color_gr = False)
ul_observed.show_lightcurve(label =False, color_BV = False, color_gr = False)
plt.xticks([59527.5, 59532.5, 59537.5], [59527.5, 59532.5, 59537.5])

for filt_, amplitude, alpha in zip(filtset, amplitudelist, alphalist):
    filt_tbl = fit_tbl[fit_tbl['filter'] == filt_]
    filt_show_tbl = show_tbl[show_tbl['filter'] == filt_]
    ul_show_tbl = UL_tbl[UL_tbl['filter'] == filt_]
    xrange = np.arange(exptime, fit_maxcut, 0.1)
    plt.ylim(22, 7)
    plt.xlabel('days [mjd]')
    plt.ylabel('Apparent magnitude')

    #plt.scatter(ul_show_tbl['obsdate'], ul_show_tbl['mag']+offset_key[filt_], facecolor = 'none', edgecolor = color_key[filt_], alpha = 0.3)
    #plt.scatter(filt_show_tbl['obsdate'], filt_show_tbl['mag']+offset_key[filt_], facecolor = 'none', edgecolor = color_key[filt_])
    #plt.errorbar(filt_show_tbl['obsdate'], filt_show_tbl['mag']+offset_key[filt_], filt_show_tbl['e_mag'], fmt ='none', color = color_key[filt_], capsize = 3)
    #plt.axvline(fit_mincut, c ='k', linestyle = '-.', linewidth = 0.5)
    #plt.axvline(fit_maxcut, c ='k', linestyle = '-.', linewidth = 0.5)
    #plt.fill_betweenx([7,22], fit_mincut, fit_maxcut, color = 'k', alpha =0.05)
    plt.plot(xrange, flux_to_mag(fireball_model(xrange, amplitude, exptime, alpha))+offset_key[filt_], linestyle = '--', linewidth = 0.7,  c=color_key[filt_])#, label = rf'[{label_row[filt_]}]$\alpha = %.2f$'%alpha)
    #plt.xlim(-2, 10)
    chisq = get_chisq_fireball(filt_tbl, fireball_model, exptime, amplitude, alpha)
    #plt.text(59536.8, flux_to_mag(fireball_model(xrange[-1], amplitude, exptime, alpha))+offset_key[filt_]-0.3, s = rf'$\chi^2$ = {round(chisq,2)}')
    #plt.text(exptime+0.7, 21, s = r'$t_{exp}$ = %.3f'%exptime)
    print('%s=%.2f'%(filt_, alpha))
#plt.legend(loc = 2)

#%% Fireball + companion interaction

# Fireball + Companion interaction(K10) model
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux

def get_chisq_combined(filt_tbl, fireball_model, comp_interac_model, explosion_Ni, explosion_inter, amplitude, alpha, rstar):
    obs_mjd,obs_mag,obs_magerr =  filt_tbl['obsdate'],filt_tbl['mag'],filt_tbl['e_mag']
    filter_ = filt_tbl['filter'][0]
    obs_flux = mag_to_flux(obs_mag)
    obs_fluxerr = obs_magerr*obs_flux*2.303/2.5
    fireballmodel_flux = fireball_model(obs_mjd, amplitude, explosion_Ni, alpha)
    interactionmodel_mag = np.array(comp_interac_model(obs_mjd-explosion_inter, rstar, filter_)) + distancemodulus
    interactionmodel_flux = mag_to_flux(np.array(interactionmodel_mag))
    model_flux = fireballmodel_flux+ interactionmodel_flux
    chisq = np.sum(((obs_flux - model_flux)/(obs_fluxerr))**2)
    reduced_chisq = chisq/len(obs_flux-2)
    return reduced_chisq

def get_multichisq_combined(theta, obs_tbl, fireball_model, comp_interac_model):
    explosion_Ni, rstar, explosion_inter = theta[0], theta[-2], theta[-1]
    ampl_cut = (len(theta))//2
    amplitudelist, alphalist = theta[1:ampl_cut], theta[ampl_cut:-2]
    sum_chisq = 0
    filtset = set(obs_tbl['filter'])
    for filt_, amplitude, alpha in zip(filtset, amplitudelist, alphalist):
        filt_tbl = obs_tbl[obs_tbl['filter'] == filt_]
        chisq = get_chisq_combined(filt_tbl, fireball_model, comp_interac_model, explosion_Ni, explosion_inter, amplitude, alpha, rstar)
        sum_chisq += chisq
    return sum_chisq
fit_tbl = vstack([all_tbl, extern_tbl])
fit_maxcut = 59537
fit_mincut = 59529
show_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut)]
fit_tbl = fit_tbl[(fit_tbl['obsdate'] < fit_maxcut) & (fit_tbl['obsdate'] > fit_mincut)]
filtset = set(fit_tbl['filter'])
amplitudelist = 1500
alphalist = 2 
explosion_Ni = 59528
explosion_inter = 59528.191
rstar = 1
initial = [explosion_Ni]+ [amplitudelist] * len(filtset)+ [alphalist] * len(filtset) + [rstar] + [explosion_inter]
bnd_explosion_Ni = [59527.5,59530]
bnd_explosion_inter = [59527.5,59529.33]
bnd_amplitude = [100,5000]
bnd_alpha = [1.4,3]
bnd_rstar = [0.5,10]
bnd = [bnd_explosion_Ni]+ [bnd_amplitude] * len(filtset)+ [bnd_alpha] * len(filtset) + [bnd_rstar] + [bnd_explosion_inter]
bnds = ()
for bound in bnd:
    bnds += (tuple(bound),)
#%%
from scipy.optimize import minimize
from Companion_interaction_K10_old import Compinteract_K10
nll = lambda *args: get_multichisq_combined(*args, obs_tbl = fit_tbl, fireball_model = fireball_model, comp_interac_model = Compinteract_K10)
soln = minimize(nll, x0 = initial, method = "SLSQP", tol = 0.1,options = dict(maxiter=3000), bounds = bnds) # tol = step_size 
#soln_combined = soln

#%%
explosion_Ni, rstar, explosion_inter = soln.x[0], soln.x[-2], soln.x[-1]
ampl_cut = (len(soln.x))//2
amplitudelist, alphalist = soln.x[1:ampl_cut], soln.x[ampl_cut:-2]
color_key, offset_key, _, _, label_row = load_filt_keys(filtset)
amplitude = amplitudelist[1]
alpha = alphalist[1]
plt.figure(figsize = (3,4.9), dpi = 500)
plt.subplots_adjust(hspace=0)
plt.xticks(visible = False)

scatterlinewidth = 1.5
mksize =100
amplitudelist
alphalist
gs = gridspec.GridSpec(nrows = 2, ncols =1, height_ratios= [10, 3], width_ratios = [6])
ax1 = plt.subplot(gs[0])
ax3 = plt.subplot(gs[1])
ax1.set_xticks([59550])
alphalist_temp = alphalist[0],alphalist[1],alphalist[2],alphalist[5],alphalist[4], alphalist[3]
amplitudelist_temp = amplitudelist[0],amplitudelist[1],amplitudelist[2],amplitudelist[5],amplitudelist[4],amplitudelist[3]

#for filt_, amplitude, alpha in zip(filtset, amplitudelist, alphalist):
for filt_, amplitude, alpha in zip(list(filtset), amplitudelist_temp , alphalist_temp):
    print(filt_, alpha)
    filt_tbl = fit_tbl[fit_tbl['filter'] == filt_]
    filt_show_tbl = show_tbl[show_tbl['filter'] == filt_]
    chisq = get_chisq_combined(filt_tbl, fireball_model, Compinteract_K10, explosion_Ni, explosion_inter, amplitude, alpha, rstar)
    xrange_fireball = np.arange(explosion_Ni+0.01, fit_maxcut, 0.1)
    xrange_K10 = np.arange(explosion_inter+0.01, fit_maxcut, 0.1)
    flux_fireball = np.nan_to_num(fireball_model(xrange_fireball, amplitude, explosion_Ni, alpha))
    flux_K10 = np.nan_to_num(mag_to_flux(np.array(Compinteract_K10(xrange_K10-explosion_inter, rstar, filt_))+distancemodulus),0)
    flux_both = flux_K10 + np.nan_to_num(fireball_model(xrange_K10, amplitude, explosion_Ni, alpha))
    
    # Observation
    #temp#######
    ax1.set_title('Fireball + Companion_interaction')
    ax1.scatter(UL_tbl_g['obsdate'], UL_tbl_g['mag']+offset_key['g'], marker = 'v', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =0.5, alpha = 0.3)
    ax1.scatter(UL_tbl_r['obsdate'], UL_tbl_r['mag']+offset_key['r'], marker = 'v', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =0.5, alpha = 0.3)
    ax1.scatter(UL_tbl_i['obsdate'], UL_tbl_i['mag']+offset_key['i'], marker = 'v', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =0.5, alpha = 0.3)
    #KCT
    ax1.scatter(KCT_g['obsdate'], KCT_g['mag']+offset_key['g'], marker = 'P', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =scatterlinewidth)
    ax1.scatter(KCT_r['obsdate'], KCT_r['mag']+offset_key['r'], marker = 'P', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.scatter(KCT_i['obsdate'], KCT_i['mag']+offset_key['i'], marker = 'P', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =scatterlinewidth)
    ax1.errorbar(KCT_g['obsdate'], KCT_g['mag']+offset_key['g'], KCT_g['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['g'])
    ax1.errorbar(KCT_r['obsdate'], KCT_r['mag']+offset_key['r'], KCT_r['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['r'])
    ax1.errorbar(KCT_i['obsdate'], KCT_i['mag']+offset_key['i'], KCT_i['e_mag'], fmt = 'none', elinewidth = 0.4, c = color_key['i'])
    #LSGT
    ax1.scatter(LSGT_g['obsdate'], LSGT_g['mag']+offset_key['g'], marker = '^', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =scatterlinewidth)
    ax1.scatter(LSGT_r['obsdate'], LSGT_r['mag']+offset_key['r'], marker = '^', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.scatter(LSGT_i['obsdate'], LSGT_i['mag']+offset_key['i'], marker = '^', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =scatterlinewidth)
    ax1.errorbar(LSGT_g['obsdate'], LSGT_g['mag']+offset_key['g'], LSGT_g['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['g'])
    ax1.errorbar(LSGT_r['obsdate'], LSGT_r['mag']+offset_key['r'], LSGT_r['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['r'])
    ax1.errorbar(LSGT_i['obsdate'], LSGT_i['mag']+offset_key['i'], LSGT_i['e_mag'], fmt = 'none', elinewidth = 0.4, c = color_key['i'])
    #RASA36
    ax1.scatter(RASA_r['obsdate'], RASA_r['mag']+offset_key['r'], marker = 's', s=mksize,facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.errorbar(RASA_r['obsdate'], RASA_r['mag']+offset_key['r'], RASA_r['e_mag'] , fmt = 'none', elinewidth = 0.4, c = color_key['r'])
    #H22
    ax1.scatter(extern_g['obsdate'], extern_g['mag']+offset_key['g'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['g'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_r['obsdate'], extern_r['mag']+offset_key['r'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['r'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_i['obsdate'], extern_i['mag']+offset_key['i'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['i'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_U['obsdate'], extern_U['mag']+offset_key['U'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['U'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_B['obsdate'], extern_B['mag']+offset_key['B'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['B'],linewidth  =scatterlinewidth)
    ax1.scatter(extern_V['obsdate'], extern_V['mag']+offset_key['V'], marker = '.', s=mksize, facecolors = 'none' ,edgecolors = color_key['V'],linewidth  =scatterlinewidth)
    #temp_end####
    #ax1.scatter(filt_show_tbl['obsdate'], filt_show_tbl['mag']+offset_key[filt_], facecolor = 'none', edgecolor = color_key[filt_])
    #ax1.errorbar(filt_show_tbl['obsdate'], filt_show_tbl['mag']+offset_key[filt_], filt_show_tbl['e_mag'], fmt ='none', color = color_key[filt_], capsize = 3)
    # Model
    ax1.plot(xrange_fireball, flux_to_mag(flux_fireball)+offset_key[filt_], linestyle = '--', linewidth = 1,  c=color_key[filt_], label = rf'[{label_row[filt_]}]$\alpha = {round(alpha,2)}$')
    ax1.plot(xrange_K10, flux_to_mag(flux_K10) + offset_key[filt_], linestyle = '-.', linewidth = 1,  c=color_key[filt_])
    ax1.plot(xrange_K10, flux_to_mag(flux_both)+ offset_key[filt_], linewidth = 1.5,  c=color_key[filt_]) 
    ax1.text(59535.8, flux_to_mag(fireball_model(xrange_fireball[-1], amplitude, explosion_Ni, alpha))+offset_key[filt_]-0.3, s = r'$\bar{\chi}^2$ = %.2f'%chisq)
    if filt_ == 'U':
        obs_U = filt_tbl
        obs_U.sort('obsdate')
        spl_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], smooth = 2)
        days_U, model_U = xrange_K10, flux_to_mag(flux_both)
    if filt_ == 'B':
        obs_B = filt_tbl
        obs_B.sort('obsdate')
        spl_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'])
        days_B, model_B = xrange_K10, flux_to_mag(flux_both)
    if filt_ == 'V':
        obs_V = filt_tbl
        obs_V.sort('obsdate')
        spl_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'])
        days_V, model_V = xrange_K10, flux_to_mag(flux_both)
    if filt_ == 'g':
        obs_g = filt_tbl
        obs_g.sort('obsdate')
        spl_g,_ = interpolate_spline(obs_g['obsdate'], obs_g['mag'])
        days_g, model_g = xrange_K10, flux_to_mag(flux_both)
    if filt_ == 'r':
        obs_r = filt_tbl
        obs_r.sort('obsdate')
        spl_r,_ = interpolate_spline(obs_r['obsdate'], obs_r['mag'])
        days_r, model_r = xrange_K10, flux_to_mag(flux_both)
    if filt_ == 'i':
        obs_i = filt_tbl
        obs_i.sort('obsdate')
        spl_i,_ = interpolate_spline(obs_i['obsdate'], obs_i['mag'])
        days_i, model_i = xrange_K10, flux_to_mag(flux_both)
ax3.scatter(obs_B['obsdate']-explosion_inter, spl_U(obs_B['obsdate'])-spl_B(obs_B['obsdate']), label = 'U-B', facecolor = 'none', edgecolor ='cyan')
ax3.scatter(obs_B['obsdate']-explosion_inter, spl_B(obs_B['obsdate'])-spl_V(obs_B['obsdate']), label = 'B-V', facecolor = 'none', edgecolor ='g')
ax3.scatter(obs_g['obsdate']-explosion_inter, spl_g(obs_g['obsdate'])-spl_r(obs_g['obsdate']), label = 'g-r', facecolor = 'none', edgecolor ='orange')


ax3.plot(days_U-explosion_inter, model_U-model_B, c = 'cyan', linewidth  =0.5)
ax3.plot(days_B-explosion_inter, model_B-model_V, c = 'g', linewidth  =0.5)
ax3.plot(days_g-explosion_inter, model_g-model_r, c = 'orange', linewidth = 0.5)
ax3.set_xlim(-1, 8.1)
ax3.set_xlabel('days since explosion')
xrange_K10
len(days_U)
len(xrange_K10)

ax3.set_ylabel('Color')
ax3.legend()
#ax3.plot(xtrange_K10,)
ax1.text(59530.5, 23.2, s = r'$t_{Ni}$ = %.4f'%(round(explosion_Ni,4)))
ax1.text(59530.5, 23.9, s = r'$t_{Comp}$ = %.4f'%(round(explosion_inter,4)))
ax1.text(59530.5, 22.5, s = r'$R_{companion}$ = %.2f $R_{\odot}$'%(round(rstar,2)))
ax1.grid()
ax1.set_ylim(25, 7)
ax1.set_xlim(explosion_inter-1, explosion_inter+8.1)
ax1.legend(loc = 4)
ax1.set_ylabel('Apparent mag')
#ax2 = ax1.twinx()
#ax2.set_ylabel('Absolute mag')
#ax2.set_ylim(25-distancemodulus,7-distancemodulus)
#ax2.legend(loc = 4)


#%%
def log_prior(theta, bounds):
    explosion_Ni, rstar, explosion_inter = theta[0], theta[-2], theta[-1]
    ampl_cut = (len(theta))//2
    amplitudelist, alphalist = theta[1:ampl_cut], theta[ampl_cut:-2]
    bnd_explosion_Ni, bnd_rstar, bnd_explosion_inter = bounds[0], bounds[-2], bounds[-1]
    ampl_cut = (len(bounds))//2
    bnd_amplitude = bounds[1:ampl_cut][0]
    bnd_alpha = bounds[ampl_cut:-2][0]
    if all(bnd_amplitude[0]<=amplitude<=bnd_amplitude[1] for amplitude in amplitudelist) and all(bnd_alpha[0]<=alpha<=bnd_alpha[1] for alpha in alphalist):
        if bnd_explosion_Ni[0]<=explosion_Ni<=bnd_explosion_Ni[1] and bnd_rstar[0]<=rstar<=bnd_rstar[1] and bnd_explosion_inter[0] <= explosion_inter <= bnd_explosion_inter[1]:
            return 0.0
    return -np.inf

def log_probability(theta, bounds, obs_tbl, fireball_model, comp_interac_model):
    lp = log_prior(theta, bounds)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_multichisq_combined(theta,  obs_tbl, fireball_model, comp_interac_model)


#%% MCMC
def sampling_MCMC(soln, obs_tbl, bounds, fireball_model, comp_interac_model, discard = 0):
    import emcee
    import corner
    pos = soln.x + ([0.3]+[100]*((len(soln.x)-3)//2)+[0.3]*((len(soln.x)-3)//2)+[0.3]+[0.3]) * np.random.randn(500, 17)
    nwalkers, ndim = pos.shape
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability, args = (obs_tbl, bounds, fireball_model, comp_interac_model)
    )
    sampler.run_mcmc(pos, 5000, progress=True)
    samples = sampler.get_chain()
    samples[discard:][:][:]
    for i, sample in enumerate(samples):
        for j, sample_1 in enumerate(sample):
            idx = [0,15,16] ######################### Chooise parameters here!
            if j == 0:
                cut_sample = sample_1[idx]
            else:
                cut_sample = np.vstack([sample_1[idx], cut_sample])
        if i == 0:
            cut_samples = cut_sample
        else:
            cut_samples = np.vstack([cut_samples, cut_sample])
    flat_samples = cut_samples
    #flat_samples = sampler.get_chain(discard=100, thin=1, flat=True)
    labels = [r"t_{Ni}", r"R_{companion}",r"t_{interact}"]
    cut_samples.shape
    fig = corner.corner(
        cut_samples, labels = labels, truths = [soln.x[0],soln.x[15],soln.x[16]]
        )
sampling_MCMC(soln, all_extern_tbl, bnds, fireball_model, Compinteract_K10, discard = 1000)





#%%





