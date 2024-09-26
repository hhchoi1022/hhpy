#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 17:34:11 2022

@author: hhchoi1022
"""
#%%
from astropy.io import ascii
from astropy.table import vstack
from astropy.table import Table
from HHsupport_analysis import load_filt_keys
from HHsupport_analysis import formatting_SNcosmo
import astropy.units as u
import sncosmo
import numpy as np
import os, glob
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import seaborn as sns
#%%
color_key, offset_key, filter_key_sncosmo, _, name_key = load_filt_keys()
#%% Refister filter
def sncosmo_registerfilter(responsefile,
                           name,
                           force = True):
    tbl = ascii.read(responsefile, format = 'fixed_width')
    band = sncosmo.Bandpass(tbl['wavelength'], tbl['response'], wave_unit = u.AA, name = name)
    sncosmo.register(band, force = force)

list_responsefile = glob.glob('/Users/hhchoi1022/Gitrepo/Config/transmission/KCT_STX16803/STX16803_sdss2_?')
KCT_filterkey = dict()
for responsefile in list_responsefile:

    filter_ = responsefile[-1]
    filename = os.path.basename(responsefile)
    sncosmo_registerfilter(responsefile, filename, force = True)
    KCT_filterkey[filter_] = filename
RASA36_responsefile = glob.glob("/Users/hhchoi1022/Gitrepo/Config/transmission/RASA36_KL4040/KL4040*")[0]
RASA_filterkey = dict()
filter_ = RASA36_responsefile[-1]
filename = os.path.basename(RASA36_responsefile)
sncosmo_registerfilter(RASA36_responsefile, filename, force = True)
RASA_filterkey[filter_] = filename

list_responsefile = glob.glob('/Users/hhchoi1022/Gitrepo/Config/transmission/LasCumbres/Las*')
Las_filterkey = dict()
for responsefile in list_responsefile:
    filter_ = responsefile[-1]
    filename = os.path.basename(responsefile)
    sncosmo_registerfilter(responsefile, filename, force = True)
    Las_filterkey[filter_] = filename
#%% Data 
from HHsupport_analysis import remove_rows_table
from observedphot import ObservedPhot
'''
file_A22 = '../data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_noextin.dat'
file_H22 = '../data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_noextin.dat'
file_IMSNG = '../data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_noextin.dat'
file_A22 = '../data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostmwextin3.10.dat'
file_H22 = '../data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat'
file_IMSNG = '../data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostmwextin3.10.dat'


#tbl_A22 = ascii.read(file_A22, format ='fixed_width')
tbl_H22 = ascii.read(file_H22, format ='fixed_width')
tbl_IMSNG = ascii.read(file_IMSNG, format ='fixed_width')
tbl_tot = vstack([tbl_IMSNG, tbl_H22])
'''
#tot_tbl = ascii.read('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Alldata_H_M_cor310.dat', format = 'fixed_width')
tot_tbl = ascii.read('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Alldata_H_M_cor310.dat', format = 'fixed_width')
#tot_tbl = ascii.read('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostmwextin3.10.dat', format = 'fixed_width')
#tot_tbl = ascii.read('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat', format = 'fixed_width')
#tot_tbl = ascii.read('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostmwextin3.10.dat', format = 'fixed_width')
observed_data = ObservedPhot(tot_tbl, MW_extinction_corrected = True, Host_extinction_corrected = True)
observed_data.exclude_observatory(['DLT40','LasCumbres0.4m','Swift', 'Swope','LasCumbres1m'])
observed_data.data['e_mag'][observed_data.data['e_mag']<0.01] = 0.03

fit_tbl = observed_data.get_data_detected()
filterkeylist = [] 
for observatory, filter_ in zip(fit_tbl['observatory'], fit_tbl['filter']):
    filter_key = filter_key_sncosmo[filter_]
    '''    
    if observatory == 'KCT':
        filter_key = KCT_filterkey[filter_]
    if observatory == 'RASA36':
        filter_key = RASA_filterkey[filter_]
    if (observatory == 'LasCumbres1m') & (filter_ in 'UBVRI'):
        filter_key = Las_filterkey[filter_]
    '''
    filterkeylist.append(filter_key)

fit_tbl['filter_sncosmo'] = filterkeylist
show_tbl = fit_tbl
#%%
phase_min = 59533
phase_max = 59590
remove_filter = ['V']#['U', 'B','V']
fit_tbl = fit_tbl[(fit_tbl['obsdate'] > phase_min) & (fit_tbl['obsdate'] < phase_max)]
remove_rows_table(fit_tbl, column_key='filter', remove_keys= remove_filter )
formatted_fit_tbl = formatting_SNcosmo(fit_tbl['obsdate'], fit_tbl['mag'], fit_tbl['e_mag'], fit_tbl['filter_sncosmo'], magsys = 'ab')

#%%SNcosmo
#chi-square
def calc_chisq(formatted_fit_tbl, filt_):
    chisquare = 0
    data = formatted_fit_tbl[formatted_fit_tbl['band'] == filt_]
    for obsdate in data['mjd']:
        flux = fitted_model.bandflux(filt_,obsdate, zp = 25, zpsys = 'ab')
        obs_flux = data[data['mjd'] == obsdate]['flux'][0]
        obs_fluxerr = data[data['mjd'] == obsdate]['fluxerr'][0]
        delflux = obs_flux - flux
        pull = delflux/obs_fluxerr
        chisquare += pull**2
    reduced_chisq = chisquare/(len(data)-2)
    return reduced_chisq

source = sncosmo.get_source('salt2')
model = sncosmo.Model(source=source)
dust = sncosmo.CCM89Dust()
import sfdmap
dustmap = sfdmap.SFDMap("/Users/hhchoi1022/Gitrepo/Config/sfddata-master")
ebv = dustmap.ebv(64.9708333, -54.9480556)
#model.add_effect(dust, 'mw', 'obs')
#model.set(mwebv = ebv)
#model.set(mwr_v = 3.1)
#model.add_effect(dust, 'host', 'rest')
#model.set(hostebv = 0.0)
#model.set(hostr_v = 3.1)
model.set(z = 0.005017)
result , fitted_model= sncosmo.fit_lc(
    formatted_fit_tbl, model,
    #['t0', 'amplitude'], # hsiao
    ['t0', 'x0', 'x1', 'c'], #salt2 or salt3
    bounds = {}
    )

figtext =  ''
for band in set(formatted_fit_tbl['band']):
    chisq = round(calc_chisq(formatted_fit_tbl, band),2)
    figtext +=f"$[reduced\ \  \chi^2$]{band}={chisq}\n"
sncosmo.plot_lc(formatted_fit_tbl, model=fitted_model, errors=result.errors, figtext = figtext, ncol = 3,  xfigsize = 10, tighten_ylim=True, color = 'black')
#%%stretch parameter & delmag
import numpy as np
from HHsupport_analysis import flux_to_mag
result.parameters
z = result.parameters[0]
t0 = result.parameters[1]
x0 = result.parameters[2]
x1 = result.parameters[3]
c = result.parameters[4]
e_x1 = result.errors['x1']
e_c = result.errors['c']

def saltt3_to_salt2(x1,c):
    x1_salt2 = ((0.985/0.138)*x1 - c - 0.005 - (0.985/0.138))*0.138/0.985/1.028
    c_salt2 = (c - 0.002 * x1_salt2 - 0.013)/0.985
    return x1_salt2, c_salt2
#x1, c = saltt3_to_salt2(x1, c)

param_stretch = 0.98+ 0.091*x1+ 0.003*x1**2- 0.00075*x1**3
e_param_stretch = np.sqrt((0.091*e_x1)**2+(0.003*2*e_x1)**2+(0.00075*3*e_x1)**2)
delmag = 1.09- 0.161*x1+ 0.013*x1**2- 0.00130*x1**3
e_delmag = np.sqrt((0.161*e_x1)**2+(0.013*2*e_x1)**2+(0.00130*3*e_x1)**2)
t_max = result.parameters[1]

t_range = np.arange(t_max - 15, t_max + 30, 0.01)
Bmag_range = fitted_model.bandmag('bessellb', 'ab', t_range)
magB_idx = np.argmin(Bmag_range)
magB_max = Bmag_range[magB_idx]
timeB_max = t_range[magB_idx]
Vmag_range = fitted_model.bandmag('bessellv', 'ab', t_range)
magV_idx = np.argmin(Vmag_range)
magV_max = Vmag_range[magV_idx]
timeV_max = t_range[magV_idx]
Rmag_range = fitted_model.bandmag('bessellr', 'ab', t_range)
magR_idx = np.argmin(Rmag_range)
magR_max = Rmag_range[magR_idx]
timeR_max = t_range[magR_idx]
gmag_range = fitted_model.bandmag('sdssg', 'ab', t_range)
magg_idx = np.argmin(gmag_range)
magg_max = gmag_range[magg_idx]
timeg_max = t_range[magg_idx]
rmag_range = fitted_model.bandmag('sdssr', 'ab', t_range)
magr_idx = np.argmin(rmag_range)
magr_max = rmag_range[magr_idx]
timer_max = t_range[magr_idx]
imag_range = fitted_model.bandmag('sdssi', 'ab', t_range)
magi_idx = np.argmin(imag_range)
magi_max = imag_range[magi_idx]
timei_max = t_range[magi_idx]

magB15_time = timeB_max + 15
Bmag_later = fitted_model.bandmag('bessellb', 'ab', magB15_time)
mB_15 = magB_max - Bmag_later

maxpoint = tot_tbl[np.abs(tot_tbl['obsdate']-int(t_max))<1]
e_magerr_max = np.median(maxpoint['e_mag'])

nu = magB_max - (-19.31) + 0.13 * x1 - 1.77 * c # Guy et al. 2007, SALT2 paper +- 0.03, +- 0.013, +- 0.16
e_nu = np.sqrt(e_magerr_max**2 + 0.03**2 + (np.sqrt((0.013/0.13)**2+(e_x1/x1)**2))**2 + (np.sqrt((0.16/1.77)**2+(np.abs(e_c/c))**2))**2)
#nu = magB_max - (-19.31 + 5*np.log10(const_hubble/70)) + 1.52*(param_stretch-1) - 1.57 * c # P. Astier, et al. 2005
#e_nu = np.sqrt(e_magB_max**2 + 0.03**2 + (np.sqrt((0.14/1.52)**2+(e_param_stretch/param_stretch)**2))**2 + (np.sqrt((0.15/1.57)**2+(np.abs(e_c/c))**2))**2)

const_hubble = 70 

#nu = magB_max - (-19.218) + 1.295*(param_stretch-1) - 3.181*c # Guy et al. 2010
distance = 10**((nu +5)/5) # unit = pc
e_distance = 10**((nu+e_nu +5)/5) - distance
#%%
print(f'distance = {distance}+-{e_distance}')
print(f't_max_all = t_max_all = {t_max}')
print(f't_max_filt = B={timeB_max}, V = {timeV_max}, R = {timeR_max}, g={timeg_max}, r = {timer_max}, i = {timei_max}')
print(f'mag_max = B={magB_max}, V = {magV_max}, R = {magR_max}, g={magg_max}, r = {magr_max}, i = {magi_max}')
print(f'ABSmag_max = B={magB_max-nu}, V = {magV_max-nu}, R = {magR_max-nu}, g={magg_max-nu}, r = {magr_max-nu}, i = {magi_max-nu}')
print(f'stretch_param = {param_stretch}+-{e_param_stretch}')
print(f'delmag = {mB_15}+-{e_delmag}')
print(f'magB_max = {magB_max}')
print(f'ABSmagB_max = {magB_max-nu}+-{e_magerr_max}')
print(f'nu = {nu}+-{e_nu}')
#%%
from HHsupport_analysis import flux_to_mag
#observed_data.exclude_filter(remove_filter)
plt.figure(figsize = (10, 8))
observed_data.show_lightcurve(figsize =(7, 5), color_BV = False, color_gr = False)
colors, offsets, filter_key_sncosmo, _, labels = load_filt_keys()
template_tbl = show_tbl.copy()
template_tbl.sort('obsdate')
template_tbl['magsys'] = 'AB'
template_tbl['mag'] = fitted_model.bandmag(template_tbl['filter_sncosmo'], template_tbl['magsys'], template_tbl['obsdate'])
template_tbl['observatory'] = 'Template'
filterset = list(set(template_tbl['filter']))


template_data = ObservedPhot(template_tbl, observer = 'Template')
#template_data.exclude_filter(remove_filter)

for filter_ in template_data.get_defined_filter():
    filt_data = template_data.get_filt_data(data_tbl = template_data.data,filters = template_data.get_defined_filter())[filter_]
    filt_data.sort('obsdate')
    time_range = (filt_data['obsdate'][0], filt_data['obsdate'][-1])
    x_range = np.linspace(59529.3315, phase_max, 500)
    mag = fitted_model.bandmag(filt_data['filter_sncosmo'][0], 'AB', x_range)
    plt.plot(x_range, mag+offsets[filter_], color = colors[filter_], linewidth = 2)
plt.xlim(59530 - 2, 59530+30)
plt.ylim(25, 5)
#%%
from wiserepSpectrum import WiseRepSpectrum
spec = WiseRepSpectrum(specfilekey = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/spectrum/WISeREP/ascii/*FTS*')
#%%
from matplotlib import cm
cmap = cm.get_cmap('viridis')
time_range = np.arange(59531, 59535, 1)
for i, time in enumerate(time_range):
    wl = np.arange(3000, 8000, 10)
    flux = fitted_model.flux(time, wl)
    #flux_avg = np.mean(flux[245:255])
    #flux /= flux_avg
    plt.plot(wl, flux, c=cmap(i/len(time_range)), label = time)
spec.show_spec_date(59535, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
#plt.ylim(0, 10)
plt.xlim(3000, 8000)
plt.legend()
#%%
spec.show_spec_date(59529, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
spec.show_spec_date(59530, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
spec.show_spec_date(59531, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
spec.show_spec_date(59532, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
spec.show_spec_date(59533, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')

spec.show_spec_date(59535, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
spec.show_spec_date(59538, normalize= False, show_flux_unit = 'flamb', normalize_cenwl = 5500 , color='r')
#%%
#%%
t0 = result.parameters[1]
def normalize(v):
    norm=np.linalg.norm(v)
    if norm==0:
        norm=np.finfo(v.dtype).eps
    return v/norm
sdssg = sncosmo.get_bandpass('sdssg')
KCTg = sncosmo.get_bandpass('STX16803_sdss2_g')
sdssr = sncosmo.get_bandpass('sdssr')
KCTr = sncosmo.get_bandpass('STX16803_sdss2_r')
RASAr = sncosmo.get_bandpass("KL4040_sdss_r")
sdssi = sncosmo.get_bandpass('sdssi')
KCTi = sncosmo.get_bandpass('STX16803_sdss2_i')
nspec = 50
time_range = np.linspace(t0, t0+50, 100)
wl_range = np.linspace(3500,9000,6000)
fluxlist = fitted_model.flux(time = time_range, wave = wl_range)
#[i,1,1] for i in range(len(time_range))
colors = sns.color_palette("viridis", len(time_range))
#%%
f = plt.figure(dpi = 300)
ax1 = f.add_subplot()

ax1.set_xlabel('wavelength[AA]')
ax1.set_title('SALT3 spectrum template')
for i, flux in enumerate(fluxlist):
    flux_binned = gaussian_filter(flux, 30)
    flux_binned = flux_binned/(1.1*np.max(flux_binned))
    spectrum = ax1.scatter(wl_range, flux_binned, s = 1, c = colors[i], alpha = 0.01)
#ax1.set_xlim(np.min(sdssg.wave)-200, np.max(sdssg.wave)+200)
ax1.set_ylabel('relative flux')
plt.show()
#%%

f = plt.figure(dpi = 200)
ax1 = f.add_subplot()
ax1.set_ylabel('transmission')
ax1.set_xlabel('wavelength[AA]')
ax1.set_title('g band')
ax2 = ax1.twinx()
for i, flux in enumerate(fluxlist):
    flux_binned = gaussian_filter(flux, 30)
    spectrum = ax2.scatter(wl_range, flux_binned*1e13, s = 1, c = colors[i], alpha = 0.01)
ax2.set_xlim(np.min(sdssg.wave)-200, np.max(sdssg.wave)+200)
ax2.set_ylabel('flux [ergs/s/cm^2/AA]')
ax1.plot(sdssg.wave, sdssg.trans, label = 'SDSS')
ax1.plot(KCTg.wave, KCTg.trans, label = 'KCT')
ax1.legend()
plt.show()
#%%

f = plt.figure(dpi = 200)
ax1 = f.add_subplot()
ax1.set_ylabel('transmission')
ax1.set_xlabel('wavelength[AA]')
ax1.set_title('r band')
ax2 = ax1.twinx()
for i, flux in enumerate(fluxlist):
    flux_binned = gaussian_filter(flux, 30)
    spectrum = ax2.scatter(wl_range, flux_binned*1e13, s = 1, c = colors[i], alpha = 0.01)
ax2.set_xlim(np.min(sdssr.wave)-200, np.max(sdssr.wave)+200)
ax2.set_ylabel('flux [ergs/s/cm^2/AA]')
ax2.set_ylim(0,0.9)
ax1.plot(sdssr.wave, sdssr.trans, label = 'SDSS')
ax1.plot(KCTr.wave, KCTr.trans, label = 'KCT')
ax1.plot(RASAr.wave, RASAr.trans, label = 'RASA36')
ax1.legend()
plt.show()
#%%

f = plt.figure(dpi = 200)
ax1 = f.add_subplot()
ax1.set_ylabel('transmission')
ax1.set_xlabel('wavelength[AA]')
ax1.set_title('i band')
ax2 = ax1.twinx()
for i, flux in enumerate(fluxlist):
    flux_binned = gaussian_filter(flux, 30)
    spectrum = ax2.scatter(wl_range, flux_binned*1e13, s = 1, c = colors[i], alpha = 0.01)
ax2.set_xlim(np.min(sdssi.wave)-200, np.max(sdssi.wave)+200)
ax2.set_ylabel('flux [ergs/s/cm^2/AA]')
ax2.set_ylim(0,0.6)
ax1.plot(sdssi.wave, sdssi.trans, label = 'SDSS')
ax1.plot(KCTi.wave, KCTi.trans, label = 'KCT')
ax1.legend()
plt.show()

#%%
t0 = result.parameters[1]
timegrid = np.linspace(t0-16, t0+80, 100)
#%%

import matplotlib.pyplot as plt
mag_sdssg = fitted_model.bandmag(time = timegrid, band = 'sdssg', magsys='ab')
mag_kctg = fitted_model.bandmag(time = timegrid, band = 'STX16803_sdss2_g', magsys='ab')
plt.figure(dpi = 300)
plt.xlabel('Phase')
plt.ylabel('S-correction')
plt.title('SALT3, g band')
plt.plot(timegrid-t0, mag_sdssg-mag_kctg, c = 'b', label = 'SDSS-KCT')
plt.axhline(0, alpha = 0.4, linestyle = '--', c = 'k')
plt.legend()
plt.show()

mag_sdssr = fitted_model.bandmag(time = timegrid, band = 'sdssr', magsys='ab')
mag_kctr = fitted_model.bandmag(time = timegrid, band = 'STX16803_sdss2_r', magsys='ab')
mag_rasar = fitted_model.bandmag(time = timegrid, band = "KL4040_sdss_r", magsys='ab')
import matplotlib.pyplot as plt
plt.figure(dpi = 300)
plt.title('SALT3, r band')
plt.xlabel('Phase')
plt.ylabel('S-correction')
plt.plot(timegrid-t0, mag_sdssr-mag_kctr, c = 'b', label = 'SDSS-KCT')
plt.plot(timegrid-t0, mag_kctr-mag_rasar, c = 'r', label = 'KCT-RASA36')
plt.plot(timegrid-t0, mag_sdssr-mag_rasar, c = 'k', label = 'SDSS-RASA36')
plt.axhline(0, alpha = 0.4, linestyle = '--', c = 'k')
plt.legend()
plt.show()

mag_sdssi = fitted_model.bandmag(time = timegrid, band = 'sdssi', magsys='ab')
mag_kcti = fitted_model.bandmag(time = timegrid, band = 'STX16803_sdss2_i', magsys='ab')
plt.figure(dpi = 300)
plt.xlabel('Phase')
plt.ylabel('S-correction')
plt.title('SALT3, i band')
plt.axhline(0, alpha = 0.4, linestyle = '--', c = 'k')
plt.plot(timegrid-t0, mag_sdssi-mag_kcti, c = 'b', label = 'SDSS-KCT')
plt.legend()
plt.show()





#plt.plot(timegrid, mag_kctg, c = 'k')




#%%stretch parameter & delmag
import numpy as np
from HHsupport_analysis import flux_to_mag
result.parameters
z = result.parameters[0]
t0 = result.parameters[1]
x0 = result.parameters[2]
x1 = result.parameters[3]
c = result.parameters[4]
e_x1 = result.errors['x1']
e_c = result.errors['c']
param_stretch = 0.98+ 0.091*x1+ 0.003*x1**2+ 0.00075*x1**3
e_param_stretch = np.sqrt((0.091*e_x1)**2+(0.003*2*e_x1)**2+(0.00075*3*e_x1)**2)
delmag = 1.09- 0.161*x1+ 0.013*x1**2- 0.00130*x1**3
e_delmag = np.sqrt((0.161*e_x1)**2+(0.013*2*e_x1)**2+(0.00130*3*e_x1)**2)
t_max = result.parameters[1]
magB_max = fitted_model.bandmag('bessellb', 'ab', t_max)
maxpoint = tot_tbl[np.abs(tot_tbl['obsdate']-int(t_max))<1]
e_magB_max = np.median(maxpoint['e_mag'])
const_hubble = 70 
#nu = magB_max - (-19.31 + 5*np.log10(const_hubble/70)) + 1.52*(param_stretch-1) - 1.57 * c # P. Astier, et al. 2005
e_nu = np.sqrt(e_magB_max**2 + 0.03**2 + (np.sqrt((0.14/1.52)**2+(e_param_stretch/param_stretch)**2))**2 + (np.sqrt((0.15/1.57)**2+(np.abs(e_c/c))**2))**2)
nu = magB_max - (-19.31) + 0.13 * x1 - 1.77 * c # Guy et al. 2007, SALT2 paper +- 0.03, +- 0.013, +- 0.16
e_nu = np.sqrt(e_magB_max**2 + 0.03**2 + (np.sqrt((0.013/0.13)**2+(e_x1/x1)**2))**2 + (np.sqrt((0.16/1.77)**2+(np.abs(e_c/c))**2))**2)

#nu = magB_max - (-19.218) + 1.295*(param_stretch-1) - 3.181*c # Guy et al. 2010
distance = 10**((nu +5)/5) # unit = pc
e_distance = 10**((nu+e_nu +5)/5) - distance
#%%
print(f'distance = {distance}+-{e_distance}')
print(f't_max = {t0}')
print(f'stretch_param = {param_stretch}+-{e_param_stretch}')
print(f'delmag = {delmag}+-{e_delmag}')
print(f'magB_max = {magB_max}+-{e_magB_max}')
print(f'nu = {nu}+-{e_nu}')
#%%
m_ej = 1.253 + 0.172 * x1
e_m_ej = np.sqrt(0.022**2 + (0.021*e_x1)**2)
m_ni = 0.478 + 0.1 * x1
e_m_ni = np.sqrt(0.023**2 + (0.02*e_x1)**2)
magB_max-nu
e_nu
e_magB_max
magB_max
t_max
e_nu
delmag
e_delmag
e_nu
magB_max - nu
distance
e_distance
#%% Early light curve fitting for estimation of first light time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
marker_key = dict(KCT='o',
                  RASA36='^',
                  LSGT='x',
                  LasCumbres1m = 'd')
from scipy.optimize import curve_fit
from matplotlib import gridspec
def model_fireball(t,a,t0):
    flux = a * (t - t0)**2
    return flux
def model_fireball_alpha(t, a, t0, alpha):
    flux = a * (t - t0 )**alpha
    return flux

def fit_fireball(obs_tbl, model, filter_, init_guess, timerange = [0,10], show = True):
    
    obs_tbl = obs_tbl[(obs_tbl['filter'] == filter_) &  (obs_tbl['status'] == 'detected')]
    timecut = [np.min(obs_tbl['obsdate'])+timerange[0],np.min(obs_tbl['obsdate'])+timerange[1]]
    early_tbl = obs_tbl[(obs_tbl['obsdate'] >= timecut[0]) & (obs_tbl['obsdate'] < timecut[1])]
    early_tbl.sort('obsdate')
    early_tbl['flux'] = 10**((early_tbl['mag'] -25)/-2.5)
    early_tbl['fluxerr'] = early_tbl['e_mag']*early_tbl['flux']*2.303/2.5
    early_tbl['observatory'][10]
    early_tbl['obsdate','mag','observatory']
    popt, pcov = curve_fit(model, early_tbl['obsdate'], early_tbl['flux'], sigma = early_tbl['fluxerr'], p0 = init_guess )
    if show:
        color_key, offset_key, _, _, label_key = load_filt_keys(['g','r','i','B','V'])
        for obsgroup in obs_tbl.group_by('observatory').groups:
            xrange = np.arange(np.min(obs_tbl['obsdate']),timecut[1], 0.1)
            observatory = obsgroup['observatory'][0]  
            plt.xlim(xrange[0]-1,xrange[-1])
            plt.scatter(obsgroup['obsdate'], obsgroup['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_], label = f'[{observatory}]{label_key[filter_]}', marker = marker_key[observatory])
            plt.errorbar(obsgroup['obsdate'], obsgroup['mag']+offset_key[filter_], obsgroup['e_mag'], c = color_key[filter_], fmt= 'none')
            plt.plot(xrange, -2.5*np.log10(model(xrange, *popt))+25+offset_key[filter_], c= color_key[filter_], linestyle = '--', label = rf'fireball[$\alpha = {popt[2]}$]')


fit_fireball(extern_B, model_fireball_alpha, 'B', [1000, 59529, 2], timerange = [0,])    
plt.legend()    
#%%


p0_fireball = [1000, 59529]
p0_fireball_alpha = [1000, 59529, 2]

xrange = np.arange(timecut[0],timecut[1], 0.1)
popt_fireball_2, pcov_fireball_2 = curve_fit(model_fireball, early_tbl_r['obsdate'], early_tbl_r['flux'], sigma = early_tbl_r['fluxerr'], p0 = p0_fireball )
popt_fireball_alpha, pcov_fireball_alpha = curve_fit(model_fireball_alpha, early_tbl_r['obsdate'], early_tbl_r['flux'], sigma = early_tbl_r['fluxerr'], p0 = p0_fireball_alpha )

plt.figure(dpi = 300)
for obsgroup in early_tbl_r.group_by('observatory').groups:
    observatory = obsgroup['observatory'][0]    
    plt.scatter(obsgroup['obsdate'], obsgroup['mag'], facecolor = 'none', edgecolor = 'k', label = observatory, marker = marker_key[observatory])
    plt.errorbar(obsgroup['obsdate'], obsgroup['mag'], obsgroup['e_mag'], fmt= 'none', c= 'k', label = observatory, marker = marker_key[observatory])
plt.plot(xrange, -2.5*np.log10(model_fireball_alpha(xrange, popt_fireball_alpha[0], popt_fireball_alpha[1], popt_fireball_alpha[2]))+25, c= 'r', linestyle = '--', label = rf'fireball[$\alpha = {round(popt_fireball_alpha[2],1)}$]')
plt.plot(xrange, -2.5*np.log10(model_fireball(xrange, popt_fireball_2[0], popt_fireball_2[1]))+25, c= 'r', linestyle = '--', label = rf'fireball[$\alpha$ = 2]')
plt.text(xrange[0]+0.3, 12.5, s = rf'$t_f$ = {round(float(popt_fireball_alpha[1]),4)}')
plt.ylim(18,12)
#plt.plot(xrange, model_fireball(xrange, popt_fireball_2[0], popt_fireball_2[1]), c= 'k', linestyle = '--', label = r'fireball[$\alpha = 2$]')
#plt.plot(xrange, model_fireball_alpha(xrange, popt_fireball_alpha[0], popt_fireball_alpha[1], popt_fireball_alpha[2]), c= 'r', linestyle = '--', label = rf'fireball[$\alpha = {round(popt_fireball_alpha[2],1)}$]')
plt.legend()
#%%
fitted_model.bandmag('sdssg', 'ab', 59530)
tot_tbl = vstack([tbl1,tbl2, tbl3])

tot_tbl = tot_tbl[(tot_tbl['obsdate'] > 59529.33) & (tot_tbl['obsdate'] < 59529.33+30)]
tot_tbl.sort('obsdate')
plt.figure(dpi = 300, figsize = (8,6))
for filt_tbl in tot_tbl.group_by('filter').groups:
    filter_ = filt_tbl['filter'][0]
    if filter_ in ['g','r','i','B','V']:
        plt.scatter(filt_tbl['obsdate'], filt_tbl['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_], label = label_key[filter_])
        plt.errorbar(filt_tbl['obsdate'], filt_tbl['mag']+offset_key[filter_], filt_tbl['e_mag'], fmt = 'none', c = color_key[filter_])
        mag_salt2 = fitted_model.bandmag(name_key[filter_], 'ab', filt_tbl['obsdate'])
        #plt.plot(filt_tbl['obsdate'], mag_salt2+offset_key[filter_], c = color_key[filter_])
plt.legend()




#%%


names = np.array(['Ia SNe', 'Tully-Fisher', 'IRAS', 'TRGB'])
dis = (np.array([16400000, 18000000, 21300000, 17900000]))
e_mag = np.array([0.22, 0.23, 0.80, 0.49])
mu = np.array([31.07, 31.28, 31.64, 31.27 ])
e_dis = (10**((mu+e_mag+5)/5)-dis)
dis
e_dis
#%%
plt.figure(dpi = 300)
plt.scatter(names, dis/1e6, c = 'k')
plt.errorbar(names, dis/1e6, e_dis/1e6, fmt = 'none', c ='k', capsize = 3)
plt.axhline(dis[0]/1e6, c= 'k', linestyle = '--', linewidth = 0.5)
plt.scatter(names[0], dis[0]/1e6, c = 'r', label = 'This work')
plt.errorbar(names[0], dis[0]/1e6, e_dis[0]/1e6, fmt = 'none', c ='r', capsize = 3)
plt.legend()
plt.ylim(10, 35)
plt.ylabel('Disance [Mpc]')
#%%