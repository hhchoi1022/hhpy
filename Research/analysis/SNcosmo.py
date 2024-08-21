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
import astropy.units as u
import sncosmo
import numpy as np
import os, glob
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import seaborn as sns
from Research.helper import Helper
#%%
helper = Helper()
color_key, offset_key, filter_key_sncosmo, _, name_key = helper.load_filt_keys()
#%% Filter registration
def registerfilter(responsefile,
                   name,
                   force = True):
    tbl = ascii.read(responsefile, format = 'fixed_width')
    band = sncosmo.Bandpass(tbl['wavelength'], tbl['response'], wave_unit = u.AA, name = name)
    sncosmo.register(band, force = force)

# KCT 
list_responsefile = glob.glob('./transmission/KCT_STX16803/STX16803_sdss2_?')
KCT_filterkey = dict()
for responsefile in list_responsefile:
    filter_ = responsefile[-1]
    filename = os.path.basename(responsefile)
    registerfilter(responsefile, filename, force = True)
    KCT_filterkey[filter_] = filename
    
# RASA36
RASA36_responsefile = glob.glob("./transmission/RASA36_KL4040/KL4040*")[0]
RASA_filterkey = dict()
filter_ = RASA36_responsefile[-1]
filename = os.path.basename(RASA36_responsefile)
registerfilter(RASA36_responsefile, filename, force = True)
RASA_filterkey[filter_] = filename

# Lascumbres1m
list_responsefile = glob.glob('./transmission/LasCumbres/Las*')
Las_filterkey = dict()
for responsefile in list_responsefile:
    filter_ = responsefile[-1]
    filename = os.path.basename(responsefile)
    registerfilter(responsefile, filename, force = True)
    Las_filterkey[filter_] = filename
#%% Data 
from observedphot import ObservedPhot
obs_tbl = ascii.read('/data1/supernova_rawdata/SN2021aefx/photometry/all_phot_MW_dereddening_Host_dereddening.dat', format = 'fixed_width')
# obs_tbl = ascii.read('/data1/supernova_rawdata/SN2021aefx/photometry/all_phot.dat', format = 'fixed_width')
observed_data = ObservedPhot(obs_tbl, target = 'SN2021aefx')
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
remove_filter = ['B','V']
fit_tbl = fit_tbl[(fit_tbl['obsdate'] > phase_min) & (fit_tbl['obsdate'] < phase_max)]
fit_tbl = helper.remove_rows_table(fit_tbl, column_key='filter', remove_keys= remove_filter )
formatted_fit_tbl = helper.SNcosmo_format(fit_tbl['obsdate'], fit_tbl['mag'], fit_tbl['e_mag'], fit_tbl['filter_sncosmo'], magsys = fit_tbl['magsys'], zp = 25)

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

#source = sncosmo.get_source('salt3')
model = sncosmo.Model(source='salt3')
dust = sncosmo.CCM89Dust()
import sfdmap
dustmap = sfdmap.SFDMap("./sfddata-master")
# ebv = dustmap.ebv(64.9708333, -54.9480556)
# model.add_effect(dust, 'mw', 'obs')
# model.set(mwebv = ebv)
# model.set(mwr_v = 3.1)
# model.add_effect(dust, 'host', 'rest')
# model.set(hostebv = 0.097)
# model.set(hostr_v = 2.3)
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
sncosmo.plot_lc(formatted_fit_tbl, model=fitted_model, errors=result.errors, figtext = figtext, ncol = 3,  xfigsize = 10, tighten_ylim=False, color = 'black')

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

maxpoint = fit_tbl[np.abs(fit_tbl['obsdate']-int(t_max))<1]
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
names = np.array(['Ia SNe', 'Tully-Fisher', 'IRAS', 'TRGB'])
dis = (np.array([distance, 18000000, 21300000, 17900000]))
e_mag = np.array([e_nu, 0.23, 0.80, 0.49])
nu = np.array([nu, 31.28, 31.64, 31.27 ])
e_dis = (10**((nu+e_mag+5)/5)-dis)
#%%
plt.figure(dpi = 300)
plt.title('Distance to NGC1566')
plt.scatter(names, dis/1e6, c = 'k')
plt.errorbar(names, dis/1e6, e_dis/1e6, fmt = 'none', c ='k', capsize = 3)
plt.axhline(dis[0]/1e6, c= 'k', linestyle = '--', linewidth = 0.5)
plt.scatter(names[0], dis[0]/1e6, c = 'r', label = 'This work')
plt.errorbar(names[0], dis[0]/1e6, e_dis[0]/1e6, fmt = 'none', c ='r', capsize = 3)
plt.legend()
plt.ylim(10, 35)
plt.ylabel('Disance [Mpc]')
#%%