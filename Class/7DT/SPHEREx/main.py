#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:01:30 2022

@author: hhchoi1022
"""
#%%
from sensitivity_7DT import synphot_7DT
from sensitivity_SPHEREx import synphot_SPHEREx
from sensitivity_KMTNet import synphot_KMTNet
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
c = 3*1e10
from astropy.table import Table
#%% Synthetic photometry for the spectrum 
tbl = Table.read('../data/ELCOSMOS_v1.part.fits')
#ID = [460674, 434152, 358577, 902734, 266256, 887708, 587750, 651347, 565897, 596510]
#ID = [793588,838914,805512,737498,897815,897774,801223,879295,840592,358577]
ID = tbl['ID']
#target_id = 1


def calc_mag_all(target_id):
    target = tbl[tbl['ID'] == ID[target_id]]
    spectrum = Table.read(f'../data/ELCOSMOS/{target["sedfile"][0]}')
    wl_AA_obs = spectrum['wavelength']
    f_nu_obs = spectrum['flux']  * wl_AA_obs * wl_AA_obs*1e-8 / c # [erg/s/cm^2/Hz]
    
    # Setting the exposure time
    texp_3yr = 180 * 0.5 * 365 / 7 * 3
    T_sens_7DT = synphot_7DT(wl_AA_obs, f_nu_obs, 180)
    T_sens_7DT_3yr = synphot_7DT(wl_AA_obs, f_nu_obs, texp_3yr)
    T_sens_SPHEREx = synphot_SPHEREx(wl_AA_obs, f_nu_obs, 2, 2, 180)
    T_sens_SPHEREx_deep = synphot_SPHEREx(wl_AA_obs, f_nu_obs, 2, 3, 5400)
    T_sens_KMTNet = synphot_KMTNet(wl_AA_obs, f_nu_obs, 180)
    return wl_AA_obs, f_nu_obs, T_sens_7DT, T_sens_7DT_3yr, T_sens_SPHEREx, T_sens_SPHEREx_deep, T_sens_KMTNet
wl_AA_obs, f_nu_obs, T_sens_7DT, T_sens_7DT_3yr, T_sens_SPHEREx, T_sens_SPHEREx_deep, T_sens_KMTNet = calc_mag_all(2)

#%% Synthetic photometry result 
plt.figure(dpi = 300)
plt.plot(wl_AA_obs*1e-4, -2.5*np.log10(f_nu_obs) -48.60, linewidth = 0.5, c= 'k', alpha = 0.3)
#plt.scatter(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['mag'], s= 10, marker = '*', c= 'r')
plt.scatter(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['mag'], s= 10, marker = '*', c= 'r')
#plt.scatter(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['mag'], s = 10, marker = '*', c='y')
plt.scatter(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['mag'], s = 10, marker = '*', c='y')
#plt.scatter(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], s = 10, marker = '*', c='b')

#plt.errorbar(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['mag'], yerr = T_sens_7DT['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c = 'r')
plt.errorbar(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['mag'], yerr = T_sens_7DT_3yr['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c = 'k')
#plt.errorbar(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['mag'], yerr = T_sens_SPHEREx['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='y')
plt.errorbar(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['mag'], yerr = T_sens_SPHEREx_deep['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='y')
#plt.errorbar(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], yerr = T_sens_KMTNet['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='b')
plt.xlim(0.4,1)
plt.show()
#%%
def mag_to_flux_ujy(mag):
    flux_ujy = 10**-((mag + 48.60)/2.5) /1e-23 / 1e-6
    return flux_ujy

def extract_SPHEREx(T_sens_SPHEREx):
    wl_cen_SPHEREx = (T_sens_SPHEREx['wavelength'] * 1e4).astype(int)
    wl_flux_SPHEREx = mag_to_flux_ujy(T_sens_SPHEREx['mag'])
    wl_fluxerr_SPHEREx = wl_flux_SPHEREx / T_sens_SPHEREx['SN']
    T_sens_SPHEREx['name'] = 'SPHEREx'
    wl_name_SPHEREx = T_sens_SPHEREx['name']
    return wl_cen_SPHEREx, wl_flux_SPHEREx, wl_fluxerr_SPHEREx, wl_name_SPHEREx
def extract_SPHEREx_deep(T_sens_SPHEREx_deep):
    wl_cen_SPHEREx_deep = (T_sens_SPHEREx_deep['wavelength'] * 1e4).astype(int)
    wl_flux_SPHEREx_deep = mag_to_flux_ujy(T_sens_SPHEREx_deep['mag'])
    wl_fluxerr_SPHEREx_deep = wl_flux_SPHEREx_deep / T_sens_SPHEREx_deep['SN']
    T_sens_SPHEREx_deep['name'] = 'SPHEREx'
    wl_name_SPHEREx_deep = T_sens_SPHEREx_deep['name']
    return wl_cen_SPHEREx_deep, wl_flux_SPHEREx_deep, wl_fluxerr_SPHEREx_deep, wl_name_SPHEREx_deep
def extract_7DT_3yr(T_sens_7DT_3yr):
    wl_cen_7DT_3yr =  T_sens_7DT_3yr['wavelength'].astype(int)
    wl_flux_7DT_3yr = mag_to_flux_ujy(T_sens_7DT_3yr['mag'])
    wl_fluxerr_7DT_3yr = wl_flux_7DT_3yr / T_sens_7DT_3yr['SN']
    T_sens_7DT_3yr['name'] = '7DT'
    wl_name_7DT_3yr = T_sens_7DT_3yr['name']
    return wl_cen_7DT_3yr, wl_flux_7DT_3yr, wl_fluxerr_7DT_3yr, wl_name_7DT_3yr
def extract_7DT(T_sens_7DT):
    wl_cen_7DT =  T_sens_7DT['wavelength'].astype(int)
    wl_flux_7DT = mag_to_flux_ujy(T_sens_7DT['mag'])
    wl_fluxerr_7DT = wl_flux_7DT / T_sens_7DT['SN']
    T_sens_7DT['name'] = '7DT'
    wl_name_7DT = T_sens_7DT['name']
    return wl_cen_7DT, wl_flux_7DT, wl_fluxerr_7DT, wl_name_7DT
def extract_KMTNet(T_sens_7DT):    
    wl_cen_KMTNet = (T_sens_KMTNet['wavelength']).astype(int)
    wl_flux_KMTNet = mag_to_flux_ujy(T_sens_KMTNet['mag'])
    wl_fluxerr_KMTNet = wl_flux_KMTNet / T_sens_KMTNet['SN']
    T_sens_KMTNet['name'] = 'KMTNet'
    wl_name_KMTNet = T_sens_KMTNet['name']
    return wl_cen_KMTNet, wl_flux_KMTNet, wl_fluxerr_KMTNet, wl_name_KMTNet
wl_cen_7DT_3yr, wl_flux_7DT_3yr, wl_fluxerr_7DT_3yr, wl_name_7DT_3yr = extract_7DT_3yr(T_sens_7DT_3yr)
wl_cen_SPHEREx_deep, wl_flux_SPHEREx_deep, wl_fluxerr_SPHEREx_deep, wl_name_SPHEREx_deep = extract_SPHEREx_deep(T_sens_SPHEREx_deep)
#%%
plt.figure(dpi = 300)
plt.scatter(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['UL5_pts'], s= 10, marker = '*', c= 'r')
plt.scatter(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['UL5_pts'], s= 10, marker = '*', c= 'r')
plt.scatter(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['UL5_pts'], s = 10, marker = '*', c='y')
plt.scatter(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['UL5_pts'], s = 10, marker = '*', c='y')
plt.ylim(24, 18)
plt.xscale('log')
T_sens_SPHEREx
plt.show()


#%% Synthetic photometry result in FLUX
plt.figure(dpi = 300)
#plt.plot(wl_AA_obs*1e-4, -2.5*np.log10(f_nu_obs) -48.60, linewidth = 0.5, c= 'k', alpha = 0.3)
#plt.scatter(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['mag'], s= 10, marker = '*', c= 'r')
plt.scatter(T_sens_7DT_3yr['wavelength']*1e-4,wl_flux_7DT_3yr, s= 10, marker = '*', c= 'r')
#plt.scatter(T_sens_SPHEREx['wavelength'],wl_flux_SPHEREx, s = 10, marker = '*', c='y')
plt.scatter(T_sens_SPHEREx_deep['wavelength'],wl_flux_SPHEREx_deep, s = 10, marker = '*', c='y')
#plt.scatter(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], s = 10, marker = '*', c='b')

#plt.errorbar(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['mag'], yerr = T_sens_7DT['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c = 'r')
plt.errorbar(T_sens_7DT_3yr['wavelength']*1e-4,wl_flux_7DT_3yr, yerr = wl_fluxerr_7DT_3yr, fmt = 'none', elinewidth = 2, capsize = 3, c = 'k')
#plt.errorbar(T_sens_SPHEREx['wavelength'],wl_flux_SPHEREx, yerr = wl_fluxerr_SPHEREx, fmt = 'none', elinewidth = 2, capsize = 3, c='y')
plt.errorbar(T_sens_SPHEREx_deep['wavelength'],wl_flux_SPHEREx_deep, yerr = wl_fluxerr_SPHEREx_deep, fmt = 'none', elinewidth = 2, capsize = 3, c='y')
#plt.errorbar(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], yerr = T_sens_KMTNet['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='b')

plt.xlim(0.4,1)
plt.ylim(4e0,1e2)
plt.yscale('log')
plt.show()

#%% Setting the spectrm 
#wl_cen = wl_cen_SPHEREx_deep.tolist()+ wl_cen_7DT_3yr.tolist()
#wl_cen = wl_cen_7DT_3yr
wl_cen = wl_cen_SPHEREx_deep
#wl_flux = wl_flux_SPHEREx_deep.tolist()+ wl_flux_7DT_3yr.tolist() 
#wl_flux = wl_flux_7DT_3yr
wl_flux = wl_flux_SPHEREx_deep
#wl_fluxerr = wl_fluxerr_SPHEREx_deep.tolist()+ wl_fluxerr_7DT_3yr.tolist()
#wl_fluxerr = wl_fluxerr_7DT_3yr
wl_fluxerr = wl_fluxerr_SPHEREx_deep
#wl_name = wl_name_SPHEREx_deep.tolist() + wl_name_7DT_3yr.tolist()
#wl_name = wl_name_7DT_3yr
wl_name = wl_name_SPHEREx_deep


path = '/Users/hhchoi1022/Program/eazy-photoz/inputs/'
#%%
def make_translate(wl_cen, wl_name, path):
# Make header and translate file
    if len(set(wl_name)) == 2:
        obskey = '7DT_SPHEREx'
    else:
        obskey = wl_name[0]
    filtlist = []
    errfiltlist = []
    filt_key_list = []
    errfiltkeylist = []
    for i, value in enumerate(zip(wl_cen, wl_name)):
        cen, name = value
        filt_name = f'{name}_{int(cen)}'
        err_filt_name = f'e_{name}_{int(cen)}'
        filt_key = f'F{i+1}'
        err_key = f'E{i+1}'
        filtlist.append(filt_name)
        errfiltlist.append(err_filt_name)
        filt_key_list.append(filt_key)
        errfiltkeylist.append(err_key)
    col1 = filtlist + errfiltlist
    col2 = filt_key_list + errfiltkeylist
    result = Table()
    result['col1'] = col1
    result['col2'] = col2
    result.write(f'{path}{obskey}.translate', format = 'ascii.fast_no_header')  
#make_translate(wl_cen, wl_name, path)
#%%
def main():
    namerow = '# id'
    keyrow = '# id'
    targetrow = ''
    for target_id in range(len(ID)):
        row = f'     {target_id+1}'
        wl_AA_obs, f_nu_obs, T_sens_7DT, T_sens_7DT_3yr, T_sens_SPHEREx, T_sens_SPHEREx_deep, T_sens_KMTNet = calc_mag_all(target_id)
        wl_cen_7DT, wl_flux_7DT, wl_fluxerr_7DT, wl_name_7DT = extract_7DT_3yr(T_sens_7DT_3yr)
        wl_cen_SPHEREx, wl_flux_SPHEREx, wl_fluxerr_SPHEREx, wl_name_SPHEREx = extract_SPHEREx_deep(T_sens_SPHEREx_deep)
        wl_cen = wl_cen_SPHEREx
        wl_flux = wl_flux_SPHEREx
        wl_fluxerr = wl_fluxerr_SPHEREx 
        wl_name = wl_name_SPHEREx
        
        #wl_cen = wl_cen_SPHEREx.tolist() + wl_cen_7DT.tolist() 
        #wl_flux = wl_flux_SPHEREx.tolist() + wl_flux_7DT.tolist() 
        #wl_fluxerr = wl_fluxerr_SPHEREx.tolist() + wl_fluxerr_7DT.tolist() 
        #wl_name = wl_name_SPHEREx.tolist() + wl_name_7DT.tolist() 

        if len(set(wl_name)) == 2:
            obskey = '7DT_SPHEREx'
        else:
            obskey = wl_name[0]
        for i, value in enumerate(zip(wl_cen, wl_flux, wl_fluxerr, wl_name)):
            cen, flux, fluxerr, name = value
            filt_name = f'{name}_{int(cen)}'
            err_filt_name = f'e_{name}_{int(cen)}'
            filt_key = f'F{i+1}'
            err_key = f'E{i+1}'
            flux_value = format(flux, '1.5e')
            fluxerr_value = format(fluxerr, '1.5e')
        
            row += f'  {flux_value}  {fluxerr_value}'
            if target_id == 0:
                namerow += f' {filt_name} {err_filt_name}'
                keyrow += f' {filt_key} {err_key}'
        row += '\n'
        targetrow += row
    namerow += '\n'
    keyrow += '\n'
    with open(f'{path}{obskey}_eazy.cat','w') as f:
        f.write(namerow)
        f.write(keyrow)
        f.write(targetrow)

main()
#%%

for id_ in ID:
    print(tbl[tbl['ID'] == id_])
ID = [460674, 434152, 358577, 902734, 266256, 887708, 587750, 651347, 565897, 596510]
speczlist = [0.02134, 0.36132, 0.59710, 0.22081, 0.3625, 0.5187, 0.6720, 0.8661, 0.1177, 1.0040]

'''
COSMOS 2783910 / 460674 0.02134 [0.00002] 16.265 [-]
[BMM2014] 545891 / 434152 0.36132 [0.00005] 20.503 [0.005]
[BMM2014] 661332 / 358577 0.59710 [-] 21.294 [0.007]
COSMOS 2783910 / 902734 0.22081 [0.00003] 21.195 [0.007]
[BMM2014] 292143 / 266256 0.3625 [-] 21.398 [0.007]
UVISTADR1 J100200.04+023838.4 / 887708 0.5187 [-] 22.22 [0.01]
[BMM2014] 983420 / 587750 0.6720 [-] 22.59 [0.01]
UVISTADR1 J095815.18+021714.2 / 651347 0.8661 [-] 22.62 [0.01]	
[BMM2014] 988127 / 565897 0.1177 [-] 22.15 [0.01]
COSMOS 1781797 / 596510 1.0040 [-] 23.26 [0.02]
tbl[tbl['ID'] == 596510]
'''
# 7DT
# id z_spec z_a z_m1 chi_a
1   -1.0000  0.025  0.023 1.079875e+01 
2   -1.0000  0.356  0.356 2.579508e+01 
3   -1.0000  0.598  0.596 2.219759e+01 
4   -1.0000  0.221  0.221 7.319649e+01 
5   -1.0000  0.362  0.361 2.443973e+01 
6   -1.0000  0.468  0.480 3.175172e+01 
7   -1.0000  0.671  0.660 3.965259e+01 
8   -1.0000  1.323  1.183 2.064394e+01 
9   -1.0000  0.122  0.123 3.223963e+01 
10  -1.0000  0.961  0.967 1.501593e+01 

#SPHEREx 
 id z_spec z_a z_m1 chi_a
1   -1.0000  0.025  0.025 2.939512e+01 
2   -1.0000  0.369  0.358 7.519076e+01 
3   -1.0000  0.590  0.589 5.328077e+01 
4   -1.0000  0.239  0.234 9.690895e+01 
5   -1.0000  0.369  0.362 8.007980e+01 
6   -1.0000  0.468  0.464 5.828341e+01 
7   -1.0000  0.655  0.656 5.452107e+01 
8   -1.0000  0.847  0.854 7.616027e+01 
9   -1.0000  0.161  0.150 7.309705e+01 
10  -1.0000  0.951  0.948 7.189164e+01 

# Both
# id z_spec z_a z_m1 chi_a
1   -1.0000  0.025  0.023 7.608544e+01 
2   -1.0000  0.349  0.347 1.294584e+02 
3   -1.0000  0.590  0.592 8.128987e+01 
4   -1.0000  0.221  0.221 2.424149e+02 
5   -1.0000  0.369  0.368 1.266841e+02 
6   -1.0000  0.468  0.469 1.103728e+02 
7   -1.0000  0.655  0.654 1.243504e+02 
8   -1.0000  0.847  0.848 9.139554e+01 
9   -1.0000  0.122  0.116 1.257949e+02 
10  -1.0000  0.951  0.949 8.638552e+01 

#%% Graph
plt.figure(dpi = 1000, figsize = (12,18))
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=1)
i = 4
for i in range(10):
    wl_AA_obs, f_nu_obs, T_sens_7DT, T_sens_7DT_3yr, T_sens_SPHEREx, T_sens_SPHEREx_deep, T_sens_KMTNet = calc_mag_all(i)
    ELCOSMOSID = ID[i]
    
    plt.subplot(10,1,i+1)
    #plt.figure(dpi = 300)
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
    plt.title(rf'$ID_{{ELCOSMOS}}$={ELCOSMOSID}')
    plt.xlabel('wavelength [$\mu$m]')
    plt.ylabel('AB mag')
    plt.plot(wl_AA_obs*1e-4, -2.5*np.log10(f_nu_obs) -48.60, linewidth = 0.3, c= 'k', alpha = 1)
    T_sens_7DT
    plt.scatter(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['mag'], s= 5, marker = '*', c= 'r')
    plt.scatter(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['mag'], s = 5, marker = '*', c='y')
    #plt.scatter(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], s = 10, marker = '*', c='b')
    err_7DT = 2.5*np.log10(1+1/T_sens_7DT['SN'])
    err_SPHEREx = 2.5*np.log10(1+1/T_sens_SPHEREx['SN'])
    plt.errorbar(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['mag'], yerr = T_sens_7DT_3yr['magerr'], fmt = 'none', elinewidth = 1, capsize = 2, c = 'r', alpha = 0.3)
    plt.errorbar(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['mag'], yerr = T_sens_SPHEREx_deep['magerr'], fmt = 'none', elinewidth = 1, capsize = 2, c='y',  alpha = 0.3)
    #plt.errorbar(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], yerr = T_sens_KMTNet['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='b')
    #plt.ylim(10,18)
    plt.xlim(0.3 ,5.2)
    #plt.xscale('log')
    #plt.xticks([0.3,0.6,0.9,1.5,2,3,4,5],labels = ['0.3','0.6','0.9','1.5','2','3','4','5'])
    plt.ylim(np.min(T_sens_SPHEREx_deep['mag'])-1,np.max(T_sens_7DT_3yr['mag'])+1)
    #plt.show()


#%% Graph
from matplotlib.ticker import StrMethodFormatter
plt.figure(dpi = 500, figsize = (12,6))
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=1)
for i in range(10):
    wl_AA_obs, f_nu_obs, T_sens_7DT, T_sens_7DT_3yr, T_sens_SPHEREx, T_sens_SPHEREx_deep, T_sens_KMTNet = calc_mag_all(i)
    ELCOSMOSID = ID[i]
    plt.subplot(5,2,i+1)
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places
    plt.title(rf'$ID_{{ELCOSMOS}}$={ELCOSMOSID}')
    if i >=8:
        plt.xlabel('wavelength [$\mu$m]')
    plt.ylabel('AB mag')
    plt.plot(wl_AA_obs*1e-4, -2.5*np.log10(f_nu_obs) -48.60, linewidth = 0.6, c= 'k', alpha = 1)
    T_sens_7DT
    #plt.scatter(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['mag'], s= 10, marker = '*', c= 'r')
    #plt.scatter(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['mag'], s = 10, marker = '*', c='y')
    #plt.scatter(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], s = 10, marker = '*', c='b')
    err_7DT = 2.5*np.log10(1+1/T_sens_7DT['SN'])
    err_SPHEREx = 2.5*np.log10(1+1/T_sens_SPHEREx['SN'])
    #plt.errorbar(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['mag'], yerr = T_sens_7DT_3yr['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c = 'r')
    #plt.errorbar(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['mag'], yerr = T_sens_SPHEREx_deep['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='y')
    #plt.errorbar(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], yerr = T_sens_KMTNet['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='b')
    #plt.ylim(10,18)
    plt.xlim(0.3,6)
    plt.xscale('log')
    plt.xticks([0.3,0.6,0.9,1.5,2,3,4,5],labels = ['0.3','0.6','0.9','1.5','2','3','4','5'])
    plt.ylim(np.min(T_sens_SPHEREx_deep['mag'])-1,np.max(T_sens_7DT_3yr['mag'])+1)
    #plt.show()


#%%
import numpy as np
def calc_rms(val1, val2):
    return np.sqrt(np.sum((np.abs(val2-val1)**2)))
  #%%
photz_result_both = ascii.read('/Users/hhchoi1022/Program/eazy-photoz/inputs/OUTPUT/7DT_SPHEREx.output')
photz_result_7DT = ascii.read('/Users/hhchoi1022/Program/eazy-photoz/inputs/OUTPUT/7DT.output')
photz_result_SPHEREx = ascii.read('/Users/hhchoi1022/Program/eazy-photoz/inputs/OUTPUT/SPHEREx.output')
photz_result_both['z_spec'] = speczlist
photz_result_7DT['z_spec'] = speczlist
photz_result_SPHEREx['z_spec'] = speczlist
rms_both = calc_rms(speczlist, photz_result_both['z_a'])
rms_7DT = calc_rms(speczlist, photz_result_7DT['z_a'])
rms_SPHEREx = calc_rms(speczlist, photz_result_SPHEREx['z_a'])
import matplotlib.pyplot as plt
plt.figure(dpi = 300)
plt.scatter(photz_result_both['z_a'], photz_result_both['z_spec'],marker = 'x',alpha = 0.3, label = f'7DT+SPHEREx [RMS = %.3f]'%rms_both)
plt.scatter(photz_result_7DT['z_a'], photz_result_7DT['z_spec'], marker = '.', alpha = 0.3, label = f'7DT [RMS = %.3f]'%rms_7DT)
plt.scatter(photz_result_SPHEREx['z_a'], photz_result_SPHEREx['z_spec'], marker = '^', alpha = 0.3, label =f'SPHEREx [RMS = %.3f]'%rms_SPHEREx)
plt.xlabel('photo-z')
plt.ylabel('spec-z')
plt.plot([0,1],[0,1], c= 'k', linewidth = 0.3)
plt.legend()
plt.grid()

  #%%

allzlist = tbl['ZPHOT']
photz_result_both = ascii.read('/Users/hhchoi1022/Program/eazy-photoz/inputs/OUTPUT/7DT_SPHEREx_all.output')
photz_result_7DT = ascii.read('/Users/hhchoi1022/Program/eazy-photoz/inputs/OUTPUT/7DT_all.output')
photz_result_SPHEREx = ascii.read('/Users/hhchoi1022/Program/eazy-photoz/inputs/OUTPUT/SPHEREx_all.output')

photz_result_both['z_spec'] = allzlist
photz_result_7DT['z_spec'] = allzlist
photz_result_SPHEREx['z_spec'] = allzlist
photz_result_7DT = photz_result_7DT[photz_result_7DT['z_a'] > 0]
rms_both = calc_rms(photz_result_both['z_spec'], photz_result_both['z_a'])
rms_7DT = calc_rms(photz_result_7DT['z_spec'], photz_result_7DT['z_a'])
rms_SPHEREx = calc_rms(photz_result_SPHEREx['z_spec'], photz_result_SPHEREx['z_a'])
import matplotlib.pyplot as plt
plt.figure(dpi = 300)
plt.scatter(photz_result_both['z_a'], photz_result_both['z_spec'],marker = 'x',alpha = 0.3, label = f'7DT+SPHEREx [RMS = %.3f]'%rms_both)
plt.scatter(photz_result_7DT['z_a'], photz_result_7DT['z_spec'], marker = '.', alpha = 0.3, label = f'7DT [RMS = %.3f]'%rms_7DT)
plt.scatter(photz_result_SPHEREx['z_a'], photz_result_SPHEREx['z_spec'], marker = '^', alpha = 0.3, label =f'SPHEREx [RMS = %.3f]'%rms_SPHEREx)
plt.xlabel('photo-z')
plt.ylabel('ELCOSMOS-z')
plt.plot([0,1],[0,1], c= 'k', linewidth = 0.3)
plt.legend()
plt.grid()

