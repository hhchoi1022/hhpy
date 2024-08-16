#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 00:22:07 2022

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

#%%
spectrum = Table.read('../data/ELCOSMOS/sed_426043.fits')
wl_AA_obs = spectrum['wavelength']
f_nu_obs = spectrum['flux']  * wl_AA_obs * wl_AA_obs*1e-8 / c # [erg/s/cm^2/Hz]



texp_3yr = 180 * 0.5 * 365 / 7 * 3
T_sens_7DT = synphot_7DT(spec_AA = wl_AA_obs, spec_fnu = f_nu_obs, texp = 100)
T_sens_7DT_3yr = synphot_7DT(spec_AA = wl_AA_obs, spec_fnu = f_nu_obs, texp = texp_3yr)
T_sens_SPHEREx = synphot_SPHEREx(wl_AA_obs, f_nu_obs, 2, 2, 180)
T_sens_KMTNet = synphot_KMTNet(wl_AA_obs, f_nu_obs, 180)

plt.figure(dpi = 300)
plt.plot(wl_AA_obs*1e-4, -2.5*np.log10(f_nu_obs) -48.60, linewidth = 0.5, c= 'k', alpha = 0.3)
T_sens_7DT
plt.scatter(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['mag'], s= 10, marker = '*', c= 'r')
plt.scatter(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['mag'], s = 10, marker = '*', c='y')
plt.scatter(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], s = 10, marker = '*', c='b')
err_7DT = 2.5*np.log10(1+1/T_sens_7DT['SN'])
err_SPHEREx = 2.5*np.log10(1+1/T_sens_SPHEREx['SN'])
plt.errorbar(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['mag'], yerr = T_sens_7DT['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c = 'r')
plt.errorbar(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['mag'], yerr = T_sens_SPHEREx['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='y')
plt.errorbar(T_sens_KMTNet['wavelength']*1e-4,T_sens_KMTNet['mag'], yerr = T_sens_KMTNet['magerr'], fmt = 'none', elinewidth = 2, capsize = 3, c='b')
#plt.ylim(10,18)
plt.xlim(0.35,5.3)
plt.ylim(17,25)
plt.xscale('log')
plt.show()
#%%
T_sens_KMTNet
T_sens_SPHEREx
plt.figure(dpi = 300, figsize = (6,4))
plt.title('Sensitivity_sky2')
plt.ylim(24,18)
plt.xlim(0.38, 0.92)
plt.scatter(T_sens_7DT['wavelength']*1e-4,T_sens_7DT['UL5_pts'], s= 15, marker = '*', c= 'r', alpha = 0.3, label = '[7DT] Single')
plt.scatter(T_sens_7DT_3yr['wavelength']*1e-4,T_sens_7DT_3yr['UL5_pts'], s= 15, marker = '*', c= 'r', label = '[7DT] Deep ')
#plt.scatter(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['UL5_pts'], s= 15, marker = '.', alpha = 0.3, c= 'b', label = '[SPHEREx] Single')
#plt.scatter(T_sens_SPHEREx_deep['wavelength'],T_sens_SPHEREx_deep['UL5_pts'], s= 15, marker = '.', c= 'b', label = '[SPHEREx] Deep')
#plt.xscale('log')
#plt.xticks([0.3,0.6,0.9,1.5,2,3,4,5],labels = ['0.3','0.6','0.9','1.5','2','3','4','5'])

color = {'B':'b',
     'V':'g',
     'R':'r',
     'I':'y'}
filter_name = ['B','V','R','I']
#for i, filt_ in enumerate(filter_name):
#    plt.scatter(T_sens_KMTNet['wavelength'][i]*1e-4,T_sens_KMTNet['UL5_pts'][i], s = 70, marker = '+', c=f'{color[filt_]}', label = f'[KMTNet] {filt_}')


#plt.scatter(T_sens_SPHEREx['wavelength'],T_sens_SPHEREx['UL5_pts'], s = 10, marker = '*', c='b')
plt.plot(np.array([0.4250, 0.6750, 0.8250]), [20.2, 19.8, 18.9], 'x', c='blue')
plt.legend(loc=4)
#for i, value in enumerate(zip(T_sens_7DT['wavelength'],T_sens_7DT['UL5_pts'])):
    #wl, ul = value
    #if i%3 == 0:
        #plt.text(wl*1e-4-0.01, ul+0.5, int(wl), fontsize = 5)
        #plt.text(wl*1e-4-0.01, ul+3.0, int(wl), fontsize = 5)

plt.xlabel('wavelength [$\mu m$]')  
plt.ylabel('AB mag')

