#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 14:57:12 2022

@author: hhchoi1022
"""
#%%

import os, glob
import mesa_reader as mesa
import matplotlib.pyplot as plt
from astropy.table import Table
import time
import numpy as np
from matplotlib.colors import ListedColormap
from astropy.table import Table
from astropy.io import ascii
from astropy.table import Table
from scipy.optimize import curve_fit
import numpy as np

#%%
######################################################################
#################            HOMEWORK 1            ###################
######################################################################
def run_MESA(initial_mass, metalicity, run = True):
    home_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/'
    origin_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams/'
    inlist_project = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams/inlist_project'
    save_model_filename_key = f"'{initial_mass}_TAMS.mod'\n"
    filename_for_profile_when_terminate_key = f"'zams_m{initial_mass}_z{metalicity}_profile.data'\n"
    Zbase_key = f'{metalicity}d0'
    initial_mass_key = f'{initial_mass}'
    initial_z_key = f'{metalicity}d0'
    folder_key = f'zams{initial_mass}_{metalicity}/'
    
    origin_inlist = open(inlist_project, 'r')
    lines = origin_inlist.readlines()
    origin_inlist.close()
    os.system(f'cp -r {origin_dir} {home_dir}{folder_key}')
    with open(home_dir+folder_key+'inlist_project', 'w') as f:
        for line in lines:
            if 'save_model_filename' in line:
                line = line.split('=')[0] + ' = ' + save_model_filename_key
            if 'filename_for_profile_when_terminate' in line:
                line = line.split('=')[0] + ' = ' + filename_for_profile_when_terminate_key
            if 'Zbase' in line:
                line = line.split('=')[0] + ' = ' + Zbase_key
            if 'initial_mass' in line:
                line = line.split('=')[0] + ' = ' + initial_mass_key
            if 'initial_z' in line:
                line = line.split('=')[0] + ' = ' + initial_z_key
            f.write(line)
    os.chdir(home_dir+folder_key)
    if run:
        os.system('./rn')
    output_last_profile = home_dir+folder_key+f'zams_m{initial_mass}_z{metalicity}_profile.data'
    output_history = home_dir+folder_key+'LOGS/history.data'
    return output_last_profile, output_history
#%%
def read_properties_MESA(output):
    h_last = mesa.MesaData(output)
    #h_last.header_names # header data indicate the resultant physical property at the end of the model 
    #h_last.bulk_names # bulk data indicate the radial profile of each physial properties 
    
    Teff = h_last.Teff # Effective temperature
    Luminosity = h_last.photosphere_L # Luminosity
    Radius = h_last.photosphere_r # Radius
    Central_temp = 10**h_last.logT[-1] # Central temperature in log scale
    Central_density = 10**h_last.logRho[-1] # Central density in log scale
    Central_pressure = 10**h_last.logP[-1] # Central pressure in log scale
    return Teff, Luminosity, Radius, Central_temp, Central_density, Central_pressure
#%%

masslist = [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 60.0, 70.0, 100.0]
metallist = [0.02]
masslist = [2.0]
#%%
Teff_list = []
Luminosity_list = []
Radius_list = []
Central_temp_list = []
Central_density_list = []
Central_pressure_list = []
result_masslist = []
result_metallist = []
for initial_mass in masslist:
    for metalicity in metallist:
        last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
        Teff, Luminosity, Radius, Central_temp, Central_density, Central_pressure = read_properties_MESA(last_model)
        Teff_list.append(Teff)
        Luminosity_list.append(Luminosity)
        Radius_list.append(Radius)
        Central_temp_list.append(Central_temp)
        Central_density_list.append(Central_density)
        Central_pressure_list.append(Central_pressure)
        result_masslist.append(initial_mass)
        result_metallist.append(metalicity)
names = ['Teff', 'Luminosity', 'Radius', 'Central_temp', 'Central_density', 'Central_pressure', 'initial_mass', 'metalicity']
result = Table(names = names, data = [Teff_list, Luminosity_list, Radius_list, Central_temp_list, Central_density_list, Central_pressure_list, result_masslist, result_metallist])

#%% low_metalicity

masslist = [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0]
metallist = [0.001]

Teff_list = []
Luminosity_list = []
Radius_list = []
Central_temp_list = []
Central_density_list = []
Central_pressure_list = []
result_masslist = []
result_metallist = []
for initial_mass in masslist:
    for metalicity in metallist:
        last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
        Teff, Luminosity, Radius, Central_temp, Central_density, Central_pressure = read_properties_MESA(last_model)
        Teff_list.append(Teff)
        Luminosity_list.append(Luminosity)
        Radius_list.append(Radius)
        Central_temp_list.append(Central_temp)
        Central_density_list.append(Central_density)
        Central_pressure_list.append(Central_pressure)
        result_masslist.append(initial_mass)
        result_metallist.append(metalicity)
names = ['Teff', 'Luminosity', 'Radius', 'Central_temp', 'Central_density', 'Central_pressure', 'initial_mass', 'metalicity']
result_low = Table(names = names, data = [Teff_list, Luminosity_list, Radius_list, Central_temp_list, Central_density_list, Central_pressure_list, result_masslist, result_metallist])




#%%


def power_law(x, amplitude, alpha):
    y = amplitude * x ** alpha 
    return y
#%%

xrange_all = np.linspace(0, 200, 1000)
xrange_small = np.linspace(0, 2, 100)
xrange_medium = np.linspace(2, 55, 100)
xrange_large = np.linspace(55, 200, 100)

#%% ML relation
small = result[(result['initial_mass'] > 0.9) & (result['initial_mass'] < 2.1)]
medium = result[(result['initial_mass'] > 2.1) & (result['initial_mass'] < 55)]
large = result[(result['initial_mass'] < 40)]
popt_small, _ = curve_fit(power_law, small['initial_mass'], small['Luminosity'], sigma = small['Luminosity']*0.1)
popt_medium, _ = curve_fit(power_law, medium['initial_mass'], medium['Luminosity'], sigma = medium['Luminosity']*0.1)
popt_large, _ = curve_fit(power_law, large['initial_mass'], large['Luminosity'], sigma = large['Luminosity']*0.1)
popt_all, pcov_all = curve_fit(power_law, result['initial_mass'], result['Luminosity'], sigma = result['Luminosity']*0.1)
popt_low, pcov_low = curve_fit(power_law, result_low['initial_mass'], result_low['Luminosity'], sigma = result_low['Luminosity']*0.1)
perr_all = np.sqrt(np.diag(pcov_all))
perr_low = np.sqrt(np.diag(pcov_low))

plt.figure(figsize = (4,3), dpi  =600)
plt.title('Mass-Luminosity')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel(r'Luminosity [$L_\odot$]')
#plt.scatter(small['initial_mass'], small['Luminosity'])
#plt.scatter(medium['initial_mass'], medium['Luminosity'])
#plt.scatter(large['initial_mass'], large['Luminosity'])
plt.scatter(result['initial_mass'], result['Luminosity'], edgecolors = 'k', facecolors = 'none')
plt.scatter(result_low['initial_mass'], result_low['Luminosity'], edgecolors = 'r', facecolors = 'none')

#plt.plot(xrange_all, power_law(xrange_all, popt_small[0], popt_small[1]),linestyle = '--', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_small[1])
#plt.plot(xrange_medium, power_law(xrange_medium, popt_medium[0], popt_medium[1]))
#plt.plot(xrange_all, power_law(xrange_all, popt_large[0], popt_large[1]),linestyle = '--', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_large[1])
plt.plot(xrange_all, power_law(xrange_all, popt_all[0], popt_all[1]), linestyle = '--', c = 'k', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_all[1])
plt.plot(xrange_all, power_law(xrange_all, popt_low[0], popt_low[1]), linestyle = '--', c = 'r', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.001]'%popt_low[1])

plt.legend()
plt.xlim(3e-1,2e2)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0,100)
#plt.ylim(0,16)

#%% MR relation
small = result[(result['initial_mass'] > 0.9) & (result['initial_mass'] < 2.1)]
medium = result[(result['initial_mass'] > 2.1) & (result['initial_mass'] < 55)]
large = result[(result['initial_mass'] < 40)]
popt_small, _ = curve_fit(power_law, small['initial_mass'], small['Radius'], sigma = small['Radius']*0.1)
popt_medium, _ = curve_fit(power_law, medium['initial_mass'], medium['Radius'], sigma = medium['Radius']*0.1)
popt_large, _ = curve_fit(power_law, large['initial_mass'], large['Radius'], sigma = large['Radius']*0.1)
popt_all, pcov_all = curve_fit(power_law, result['initial_mass'], result['Radius'], sigma = result['Radius']*0.1)
popt_low, pcov_low = curve_fit(power_law, result_low['initial_mass'], result_low['Radius'], sigma = result_low['Radius']*0.1)
perr_all = np.sqrt(np.diag(pcov_all))
perr_low = np.sqrt(np.diag(pcov_low))
plt.figure(figsize = (4,3), dpi  =600)
plt.title('Mass-Radius')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel(r'Radius [$R_\odot$]')
#plt.scatter(small['initial_mass'], small['Radius'])
#plt.scatter(medium['initial_mass'], medium['Radius'])
#plt.scatter(large['initial_mass'], large['Radius'])
plt.scatter(result['initial_mass'], result['Radius'], edgecolors = 'k', facecolors = 'none')
plt.scatter(result_low['initial_mass'], result_low['Radius'], edgecolors = 'r', facecolors = 'none')

#plt.plot(xrange_all, power_law(xrange_all, popt_small[0], popt_small[1]),linestyle = '--', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_small[1])
#plt.plot(xrange_medium, power_law(xrange_medium, popt_medium[0], popt_medium[1]))
#plt.plot(xrange_all, power_law(xrange_all, popt_large[0], popt_large[1]),linestyle = '--', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_large[1])
plt.plot(xrange_all, power_law(xrange_all, popt_all[0], popt_all[1]), linestyle = '--', c = 'k', label = r'$R$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_all[1])
plt.plot(xrange_all, power_law(xrange_all, popt_low[0], popt_low[1]), linestyle = '--', c = 'r', label = r'$R_\odot$ $\propto$ M$^{%.2f}$ [Z = 0.001]'%popt_low[1])


plt.legend()
plt.xlim(3e-1,2e2)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0,100)
#plt.ylim(0,16)

#%% MT relation
small = result[result['initial_mass'] < 2.1]
medium = result[(result['initial_mass'] > 2.1) & (result['initial_mass'] < 55)]
large = result[(result['initial_mass'] < 40)]
popt_small, _ = curve_fit(power_law, small['initial_mass'], small['Teff'], sigma = small['Teff']*0.1)
popt_medium, _ = curve_fit(power_law, medium['initial_mass'], medium['Teff'], sigma = medium['Teff']*0.1)
popt_large, _ = curve_fit(power_law, large['initial_mass'], large['Teff'], sigma = large['Teff']*0.1)
popt_all, pcov_all = curve_fit(power_law, result['initial_mass'], result['Teff'], sigma = result['Teff']*0.1)
popt_low, pcov_low = curve_fit(power_law, result_low['initial_mass'], result_low['Teff'], sigma = result_low['Radius']*0.1)
perr_all = np.sqrt(np.diag(pcov_all))
perr_low = np.sqrt(np.diag(pcov_low))
plt.figure(figsize = (4,3), dpi  =600)
plt.title('Mass-Teff')
plt.xlabel('initial_mass [$M_\odot$]')
plt.ylabel('Teff [K]')
#plt.scatter(small['initial_mass'], small['Teff'])
#plt.scatter(medium['initial_mass'], medium['Teff'])
#plt.scatter(large['initial_mass'], large['Teff'])
plt.scatter(result['initial_mass'], result['Teff'], edgecolors = 'k', facecolors = 'none')
plt.scatter(result_low['initial_mass'], result_low['Teff'], edgecolors = 'r', facecolors = 'none')

#plt.plot(xrange_small, power_law(xrange_small, popt_small[0], popt_small[1]))
#plt.plot(xrange_medium, power_law(xrange_medium, popt_medium[0], popt_medium[1]))
#plt.plot(xrange_all, power_law(xrange_all, popt_large[0], popt_large[1]),linestyle = '--', label = r'$L$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_large[1])
plt.plot(xrange_all, power_law(xrange_all, popt_all[0], popt_all[1]), linestyle = '--', c = 'k', label = r'$T_{eff}$ $\propto$ M$^{%.2f}$ [Z = 0.02]'%popt_all[1])
plt.plot(xrange_all, power_law(xrange_all, popt_low[0], popt_low[1]), linestyle = '--', c = 'r', label = r'$T_{eff}$ $\propto$ M$^{%.2f}$ [Z = 0.001]'%popt_low[1])

plt.legend()
plt.xlim(3e-1,2e2)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0,100)
#plt.ylim(0,16)

#%% MP relation
plt.figure(figsize = (4,3), dpi  =600)
small = result[result['initial_mass'] < 2.1]
medium = result[(result['initial_mass'] > 2.1) & (result['initial_mass'] < 55)]
large = result[(result['initial_mass'] > 55)]
plt.title('Mass-Pressure')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel(r'Central pressure')
#plt.scatter(small['initial_mass'], small['Central_pressure'])
#plt.scatter(medium['initial_mass'], medium['Central_pressure'])
#plt.scatter(large['initial_mass'], large['Central_pressure'])
plt.scatter(result['initial_mass'], result['Central_pressure'], c = 'k', label = 'Z = 0.02')
plt.scatter(result_low['initial_mass'], result_low['Central_pressure'], c = 'r', label = 'Z = 0.001')
plt.legend()

plt.xscale('log')
plt.yscale('log')

#%% Mrhocen relation
plt.figure(figsize = (4,3), dpi  =300)
small = result[result['initial_mass'] < 1.6]
medium = result[(result['initial_mass'] > 2.1) & (result['initial_mass'] < 55)]
large = result[(result['initial_mass'] > 1.6)]
plt.title('Mass-Central density')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel(r'$\rho_{cen} \ [g/cm^3]$')
#plt.scatter(small['initial_mass'], small['Central_density'])
#plt.scatter(medium['initial_mass'], medium['Central_density'])
#plt.scatter(large['initial_mass'], large['Central_density'])
plt.scatter(result['initial_mass'], result['Central_density'], edgecolors = 'k', facecolors = 'none', label = 'Z = 0.02')
plt.scatter(result_low['initial_mass'], result_low['Central_density'], edgecolors = 'r', facecolors = 'none', label = 'Z = 0.001')

plt.legend()
plt.ylim(0.7, 200)
plt.xscale('log')
plt.yscale('log')


#%% MTcen relation
plt.figure(figsize = (4,3), dpi  =300)
small = result[result['initial_mass'] < 1.6]
medium = result[(result['initial_mass'] > 2.1) & (result['initial_mass'] < 55)]
large = result[(result['initial_mass'] > 1.6)]
plt.title('Mass-Central Temperature')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel(r'$T_{cen}$ [K]')
#plt.scatter(small['initial_mass'], small['Central_temp'])
#plt.scatter(medium['initial_mass'], medium['Central_temp'])
#plt.scatter(large['initial_mass'], large['Central_temp'])
plt.scatter(result['initial_mass'], result['Central_temp'], label = 'Z = 0.02', edgecolors = 'k', facecolors = 'none')
plt.scatter(result_low['initial_mass'], result_low['Central_temp'], edgecolors = 'r', facecolors = 'none', label = 'Z = 0.001')
plt.legend(loc = 4)
plt.ylim(6e6, 6e7)
plt.xscale('log')
plt.yscale('log')



#%% HR diagram
result = result[result['initial_mass']<5]
result_low = result_low[result_low['initial_mass']<5]


HRdata = Table.read('/Users/hhchoi1022/Gitrepo/config/hygdata_v3-abridged.csv')
plt.figure(dpi = 600)
plt.title('HR diagram')
plt.xlim(20000, 2000)
plt.ylim(-6, 6)
plt.ylabel(r'Luminosity [$\log_{10}(L_\odot)$]')
plt.xlabel(r'$T_{eff}$ [K]')
#plt.xscale('log')
plt.scatter(HRdata['T'], np.log10(10**(-(HRdata['absmag'] - 4.83)/2.5)), c = HRdata['T'], s = 1)
plt.scatter(result['Teff'], np.log10(result['Luminosity']), marker = 'x', label = 'Z = 0.02', c = 'r')
plt.scatter(result_low['Teff'], np.log10(result_low['Luminosity']), marker = 's', label = 'Z = 0.001', c = 'r')
plt.legend()
plt.xscale('log')
#plt.yscale('log')

#%%











######################################################################
#################            HOMEWORK 2            ###################
######################################################################

#%% constant (cgs unit)
const_k = 1.3807e-16 # [cm2 g s-2 K-1]
const_m_u = 1.6605e-24 #[g]
const_a = 7.5646e-15 # [erg cm-3 K-4]
const_msun = 1.9884e+33 # [g]
const_rsun = 6.9570e10 # [cm]
const_lsun = 3.8280e33 # [erg s-1]
const_g = 6.6743e-8 #[cm-3 g-1 s-2]
const_c = 2.9979e10 #[cm s-1]
const_Kn = 1.0035e13
const_Ke = 1.2435e15
# No reaction -> No heavy elements 
M_ch = 5.83 / (1)**2
#%%

def read_properties_MESA2(output):
    h_last = mesa.MesaData(output)
    #h_last.header_names # header data indicate the resultant physical property at the end of the model 
    #h_last.bulk_names # bulk data indicate the radial profile of each physial properties 
    
    result_tbl = Table()
    '''
    m = h_last.mass*const_msun 
    r = 10**h_last.logR * const_rsun
    t = 10**h_last.logT
    rho = 10**h_last.logRho
    p = 10**h_last.logP
    h = h_last.x_mass_fraction_H
    he = h_last.y_mass_fraction_He
    metals = h_last.z_mass_fraction_metals
    pressure_tot = h_last.pressure
    pressure_rad = h_last.prad
    pressure_gas = h_last.pgas
    mu = h_last.mu
    opacity = h_last.opacity
    luminosity = h_last.luminosity
    energy_total = h_last.total_energy
    energy_intern = h_last.energy
    energy_cellIE = h_last.cell_specific_IE
    energy_cellKE = h_last.cell_specific_KE
    frac_KEIE = h_last.cell_internal_energy_fraction
    '''
    colnames = h_last.bulk_names
    for colname in colnames:
        data = h_last.bulk_data[colname]
        result_tbl[colname] = data
    return result_tbl

#%%
masslist = [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 60.0, 70.0, 100.0]
metalicity = 0.02
plt.figure(dpi = 300)
plt.xlabel(r'M/$M_{\odot}[cumulative]$')
plt.ylabel(r'$\beta$')
betalist = []
mulist = []
for initial_mass in masslist:
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)

    #problem 1
    result_tbl['beta'] = result_tbl['pgas']/(result_tbl['prad'] + result_tbl['pgas'])
    plt.scatter(result_tbl['mass'], result_tbl['beta'], s = 1)
    beta = np.mean(result_tbl['beta'])
    mu = np.mean(result_tbl['mu'])
    print(beta)
    betalist.append(beta)
    mulist.append(mu)


x = np.linspace(betalist[0],betalist[-1],100)
y = np.log10(18.26*((1-x)/x**4)**0.5)

plt.figure(dpi = 300)
plt.scatter(betalist, np.log10(np.array(masslist)*(np.array(mulist)**2)), c = 'k', label = 'model')
plt.plot(x, y, c = 'k', linestyle = '--', label = '                 ')
plt.xlabel(r'$\beta$')
plt.ylim(-0.6, 2)
plt.ylabel(r'$\log {(\mu^2  M/M_{\odot})}$')
plt.legend()
#%% problem 2,3
'''
masslist = result_tbl['mass']*const_msun 
radiuslist = 10**result_tbl['logR']*const_rsun
shellmasslist = []
for i in range(len(masslist)):
    if i != len(masslist)-1:
        shellmasslist.append((masslist[i] - masslist[i+1])/(radiuslist[i] - radiuslist[i+1]))
    else:
        shellmasslist.append(shellmasslist[-1])


result_tbl['mass_shell'] = shellmasslist
result_tbl['potential_E'] = -const_g*masslist*result_tbl['mass_shell']/radiuslist
grav_potential = np.sum(result_tbl['potential_E'])
index = grav_potential/const_g/(  (masslist[0])**2  )*(  radiuslist[0]  ) 
'''
indexlist = []
E_gravlist = []
radiuslist = []
n_tildalist = []
for initial_mass in masslist:
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)
    
    result_tbl['potential_E'] = -const_g*result_tbl['mass']*const_msun*result_tbl['dm']/(10**result_tbl['logR']*const_rsun)
    result_tbl['potential_E']
    grav_potential = np.sum(result_tbl['potential_E'])
    index = -grav_potential/const_g/(  (result_tbl['mass'][0]*const_msun)**2  )*(  10**result_tbl['logR'][0]*const_rsun  ) 
    radius = result_tbl['logR'][0]
    indexlist.append(index)
    radiuslist.append(radius)
    E_gravlist.append(grav_potential)
    # problem3
    internal_potential = np.sum(result_tbl['cell_specific_IE']*result_tbl['dm'])
    n_tilda = -internal_potential/grav_potential * 3
    n_tildalist.append(n_tilda)
plt.figure(dpi = 300)
plt.scatter(masslist, indexlist, c = 'k')
plt.plot(masslist, indexlist, c = 'k', linestyle = '--')
plt.xlabel(r'M/$M_{\odot}[cumulative]$')
plt.ylabel(r'$\alpha$')
#plt.scatter(masslist, E_gravlist)
#%%
plt.figure(dpi = 300)
plt.scatter(masslist, n_tildalist, c = 'k')
plt.plot(masslist, n_tildalist, c = 'k', linestyle = '--')
plt.xlabel(r'M/$M_{\odot}[cumulative]$')
plt.ylabel(r'$\tilde{n}$', fontsize =13)
#%% problem 5
masslist = [5.0, 10.0, 60.0, 100.0]
for initial_mass in masslist:
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)
    result_tbl['beta'] = result_tbl['pgas']/(result_tbl['prad'] + result_tbl['pgas'])
    beta = result_tbl['beta']
    result_tbl['delta_adiabatic'] = ((1)+( ((1-beta)*(4+beta)) / (beta**2) ))/((5/2)+( (4*(1-beta)*(4+beta)) / (beta**2) ))
    result_tbl['delta_radiation'] = (3)/(16*np.pi*const_a*const_c*const_g)  *  result_tbl['pressure']/((10**result_tbl['logT'])**4)  *  result_tbl['opacity']*(result_tbl['luminosity']*const_lsun)/(result_tbl['mass']*const_msun)

    inner_tbl = result_tbl[(10**result_tbl['logR'] < np.max(10**result_tbl['logR'])/2)]
    conv_idx = np.argmin(np.abs(inner_tbl['pp'] - inner_tbl['cno']))
    plt.figure(dpi = 500, figsize = (4,4))
    plt.ylim(0,2.0)
    plt.xlim(0, np.max(10**result_tbl['logR'])*1.1)
    plt.xlabel(r'r/$R_{\odot}$', fontsize = 15)
    plt.ylabel(r'$\triangledown$', fontsize = 25)#, rotation = 0)
    plt.plot(10**result_tbl['logR'], result_tbl['delta_adiabatic'], c='k', label = r'$\triangledown_{ad}$')
    plt.plot(10**result_tbl['logR'], result_tbl['delta_radiation'], c ='k', linestyle = '--', label = r'$\triangledown_{rad}$')
    plt.scatter(10**result_tbl['logR'], result_tbl['delta_radiation'], linestyle = '--', c=result_tbl['cno']/result_tbl['pp'], s = 10)
    plt.scatter(10**inner_tbl[conv_idx]['logR'], inner_tbl[conv_idx]['delta_radiation'], linestyle = '--', c='r', s = 30, marker = '*')

    convection_radius = 10**result_tbl[(10**result_tbl['logR'] < np.max(10**result_tbl['logR'])/2) & (np.abs(result_tbl['delta_adiabatic'] - result_tbl['delta_radiation'])<0.01)][0]['logR']
    plt.fill_betweenx([0,2], convection_radius, color= 'r', alpha = 0.2)
    print(convection_radius/np.max(10**result_tbl['logR']))
    plt.legend(fontsize = 15)
    plt.colorbar(label = r'$\frac{\epsilon_{CNO}}{\epsilon_{pp}}$')

#%% problem 6a
parallel_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams1.0_0.02_parallel/LOGS/'
parallel_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams6.0_0.02_parallel/LOGS/'
h = mesa.MesaData(parallel_dir+'history.data')
10**h.log_Teff
h.star_mass
h.star_age
10**h.log_L
h.bulk_names
h.center_h1
10**h.log_R

HRdata = Table.read('/Users/hhchoi1022/Gitrepo/config/hygdata_v3-abridged.csv')

terminate_idx = -18
contractionend_idx = 48
contractionstart_idx = 28
my_cmap = plt.cm.RdBu(np.arange(plt.cm.RdBu.N))
my_cmap = ListedColormap(my_cmap)
plt.figure(dpi = 600)
plt.title(r'HR diagram [3.0$M_\odot$0.02Z]')
plt.xscale('log')
plt.xlim(1000000, 2000)
plt.ylim(-6, 6)
plt.xticks([500000, 100000,10000,3000], labels=['5e5','1e5','1e4','3e3'])
plt.ylabel(r'Luminosity [$\log_{10}(L_\odot)$]')
plt.xlabel(r'$T_{eff}$ [K]')
plt.scatter(HRdata['T'], np.log10(10**(-(HRdata['absmag'] - 4.83)/2.5)), c = 'k', s = 0.5, alpha =0.3)
plt.scatter(h.Teff[::-1], h.log_L[::-1], c = np.log10(h.star_age)[::-1], s = 10, cmap = my_cmap)
plt.colorbar(label = r'$\log{(age)}$')
#plt.scatter(h.Teff[terminate_idx], h.log_L[terminate_idx], c = 'r', s = 55, marker = '*')
#plt.scatter(h.Teff[contractionend_idx], h.log_L[contractionend_idx], c = 'r', s = 55, marker = '*')
#plt.scatter(h.Teff[contractionstart_idx], h.log_L[contractionstart_idx], c = 'r', s = 55, marker = '*')

plt.plot(h.Teff, h.log_L, linewidth  =0.5, c = 'b')

h.Teff[-18]
h.bulk_names

#%% problem 6b
contract_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams1.0_0.02_parallel/contract_LOGS/'

contractlist = glob.glob(contract_dir+'profile*.data')
contractlist.sort()
profile_data = contractlist[2]
lumlist = []
masslist = []
mulist = []
opacitylist = []
agelist = []
convectionlist = []
for profile_data in contractlist:
    h_profile = read_properties_MESA2(profile_data)
    header = mesa.MesaData(profile_data).header_data
    
    
    h_profile['pgas']/(h_profile['pgas']+h_profile['prad'])
    h_profile['beta'] = h_profile['pgas']/(h_profile['prad'] + h_profile['pgas'])
    beta = h_profile['beta']
    h_profile['delta_adiabatic'] = ((1)+( ((1-beta)*(4+beta)) / (beta**2) ))/((5/2)+( (4*(1-beta)*(4+beta)) / (beta**2) ))
    h_profile['delta_radiation'] = (3)/(16*np.pi*const_a*const_c*const_g)  *  h_profile['pressure']/((10**h_profile['logT'])**4)  *  h_profile['opacity']*(h_profile['luminosity']*const_lsun)/(h_profile['mass']*const_msun)
    
    
    lum = header['photosphere_L']
    mass = header['star_mass']
    mu = np.mean(h_profile['mu'])
    opacity = np.mean(h_profile['opacity'])
    beta = np.mean(h_profile['beta'])
    del_ad = np.mean(h_profile['delta_adiabatic'])
    del_rad = np.mean(h_profile['delta_radiation'])
    agelist.append(header['star_age'])
    lumlist.append(lum)
    masslist.append(mass)
    mulist.append(mu)
    opacitylist.append(opacity)
    convectionlist.append(del_rad/del_ad)
    plt.figure(dpi = 300)
    plt.ylim(0,2.0)
    plt.xlabel(r'r/$R_{\odot}$', fontsize = 15)
    plt.ylabel(r'$\triangledown$', fontsize = 25)#, rotation = 0)
    plt.plot(10**h_profile['logR'], h_profile['delta_adiabatic'], c='k', label = r'$\triangledown_{ad}$')
    plt.plot(10**h_profile['logR'], h_profile['delta_radiation'], c ='k', linestyle = '--', label = r'$\triangledown_{rad}$')
    plt.legend(fontsize = 15)
    np.mean(h_profile['opacity'])
    
np.array(lumlist)/(np.array(masslist)**3 * np.array(mulist)**4 * np.array(opacitylist)**1)
#%%
plt.figure(dpi = 300)
plt.scatter(np.array(agelist)[:-1]/1e6,  np.array(lumlist)[:-1]  / (np.array(masslist)[:-1]**3 * np.array(mulist)[:-1]**4 * np.array(opacitylist)[:-1]**1), c = convectionlist[:-1])
plt.yscale('log')
plt.ylabel(r'$\alpha$', fontsize = 15)
plt.xlabel('Age[Myr]')
plt.colorbar()#label = r'$\frac{\triangledown _{rad}}{\triangledown _{adia}}$')
#%% problem 6c
gplist = []
iplist = []
lmlist = []
phaselist = []
for profile_data in contractlist:
    h_profile = read_properties_MESA2(profile_data)
    header = mesa.MesaData(profile_data).header_data
    h_profile['potential_E'] = -const_g*h_profile['mass']*const_msun*h_profile['dm']/(10**h_profile['logR']*const_rsun)
    inter_potential = np.sum(h_profile['cell_specific_IE']*h_profile['dm'])
    grav_potential = np.sum(h_profile['potential_E'])
    iplist.append(inter_potential)
    gplist.append(grav_potential)
    lmlist.append(header['photosphere_L']* const_lsun)
    phaselist.append(header['star_age'])
    
dellumlist = []
delgravlist = []
deliplist = []
for i in range(len(gplist)-1):
    delgrav = gplist[i+1]-gplist[i]
    dellum = lmlist[i]
    delip = iplist[i+1]-iplist[i]
    deltime = (phaselist[i+1]- phaselist[i])*86400*365
    delgravlist.append(delgrav)
    dellumlist.append(dellum* deltime)
    deliplist.append(delip)
phaselist = phaselist[:-1]


plt.figure(dpi = 300)
plt.xlabel('Age[Myr]')
plt.ylabel(r'$\frac{\Delta E_{rad}}{\Delta E_{grav}}$', fontsize= 15)
plt.ylim(-0.35,-0.65)
plt.scatter(np.array(phaselist)/1e6, np.array(dellumlist)/np.array(delgravlist),c = convectionlist[:-1])
plt.axhline(y = -0.5, c = 'k', linestyle ='--', label = 'Virial theorem')
plt.colorbar()
plt.legend()
#plt.scatter(phaselist, np.array(deliplist)/np.array(delgravlist))

#%% problem 6d, e 
parallel_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams1.0_0.02_parallel/LOGS/'
#parallel_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/zams6.0_0.02_parallel/LOGS/'

parallellist = glob.glob(parallel_dir+'profile*.data')
parallellist


ratiolist= []
modelnumlist = []
rlist = []
agelist = []
limitlist = []
rholist = []
templist = []
pressurelist = []
for profile_data in parallellist:
    h_profile = read_properties_MESA2(profile_data)
    header = mesa.MesaData(profile_data).header_data
    pressure_electron_degenerate = const_Kn*((10**h_profile['logRho'])/(h_profile['mu']))**(5/3)
    pressure_electron_degenerate_limit = const_Ke*((10**h_profile['logRho']))**(4/3)
    pressure_central = 1.1e16*(header['star_mass']**2)*((header['photosphere_r'])**-4)
    ratio = h_profile['pressure'][-1]/pressure_electron_degenerate[-1]
    ratio_limit = h_profile['pressure'][-1]/pressure_electron_degenerate_limit[-1]
    ratio_limit = pressure_central/pressure_electron_degenerate_limit[-1]
    modelnumber = header['model_number']
    rho = h_profile['logRho'][-1]
    r = header['photosphere_r']
    age = header['star_age']
    temp = h_profile['logT'][-1]
    pressure=  h_profile['pressure'][-1]
    modelnumlist.append(modelnumber)
    agelist.append(age)
    rholist.append(rho)
    templist.append(temp)
    rlist.append(r)
    ratiolist.append(ratio)
    limitlist.append(ratio_limit)
    pressurelist.append(pressure)
#%%
# If no reaction, stars with larger than ~5.83 solar mass will become 
#plt.scatter(modelnumlist, ratiolist)
start_idx = 50
end_idx = 35
plt.figure(dpi = 300)
plt.title(r'1$M_{\odot}[Z=0.02]$ star w/o nuclear reaction')
#plt.xlim(1e-3,1e10)
plt.scatter(np.array(agelist)/1e6, limitlist, c = templist)
plt.colorbar(label = r'$\log{T}$ [K]')
plt.scatter(agelist[start_idx]/1e6, limitlist[start_idx], c = 'r', marker=  '*', s = 50)
plt.scatter(agelist[end_idx]/1e6, limitlist[end_idx], c = 'r', marker=  '*', s = 50)
#plt.xlim(0, 140)
plt.xlabel('Age[yr]')
plt.ylabel(r'$\frac{P_{gas}}{P_{e,relativistic}}$', fontsize = 15)
plt.xscale('log')


#%%
rholist_3 = rholist
templist_3 = templist
agelist_3 = agelist
#%%
plt.figure(dpi = 300)
agelist
plt.ylim(3,10)
plt.xlim(-10,10)
plt.scatter(rholist, templist, c = np.log10(agelist), marker = 'o', label = r'$1M_\odot$')
plt.scatter(rholist_3, templist_3, c = np.log10(agelist_3), marker = 'x', label = r'$3M_\odot$')
plt.colorbar(label = 'log(Age)')
plt.xlabel(r'log$\rho (g/cm^{-3})$')
plt.ylabel(r'log T (K)')
plt.legend()
agelist






######################################################################
#################            HOMEWORK 3            ###################
######################################################################
#%% Problem 1 ######################################################################
masslist = [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 60.0, 70.0, 100.0]
metalicity = 0.02
Tcenlist = []
pplist = []
cnolist = []
for initial_mass in masslist:
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)
    Tcenlist.append(10**(result_tbl['logT'][-1]))
    pplist.append(result_tbl['pp'][-1])
    cnolist.append(result_tbl['cno'][-1])

# %% Elements mass & charge
m_h = 1.008
m_h2 = 2.014
m_he3 = 3.016
m_he = 4.003
m_li = 6.015
m_li7 = 7.016
m_be7 = 7.017
m_be = 8.005
m_o = 15.999
m_c = 12.000
m_c13 = 13.003
m_n13 = 13.006
m_n = 14.003
m_n15 = 15.00
m_o15 = 15.003
m_o16 = 15.995

c_h = 1
c_he = 2
c_li = 3
c_be = 4
c_c = 6
c_n = 7
c_o = 8
# %%
def calc_reduced_mass(m1, m2):
    return m1*m2/(m1+m2)

def calc_sensitivity_temp(c1, c2, Tcen, A):
    tau = 19.72* ((c1**2 * c2 **2 * A / (Tcen/1e7))**(1/3))
    nu = (tau-2)/3
    return nu
# %% PP chain (Example for T = 1.5e7)
reduced_mass = calc_reduced_mass(m_h, m_h)
nu = calc_sensitivity_temp(c_h, c_h, 1.5e7, reduced_mass)
# %% CNO cycle (Example for T = 1.5e7)
reduced_mass = calc_reduced_mass(m_n, m_h)
nu = calc_sensitivity_temp(c_n, c_h, 1.5e7, reduced_mass)
#%% Mass - specific energy generation rate 
plt.figure(dpi  =600)

plt.plot(masslist, pplist, c= 'r', linestyle = '--', label = 'pp')
plt.plot(masslist, cnolist, c ='r', linestyle = ':', label = 'CNO')
plt.plot(masslist, np.array(pplist) + np.array(cnolist), c  = 'k', label = 'total')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel('Specific energy generation rate (erg/g/sec)')
#%% Tcen - specific energy generation rate 
plt.figure(dpi  =600)

plt.plot(Tcenlist, pplist, c= 'r', linestyle = '--', label = 'pp')
plt.plot(Tcenlist, cnolist, c ='r', linestyle = ':', label = 'CNO')
plt.plot(Tcenlist, np.array(pplist) + np.array(cnolist), c  = 'k', label = 'total')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$T_{cen}$ (K)')
plt.ylabel('Specific energy generation rate (erg/g/sec)')
#%% calculate nu(sensitivity to Temp) considering the ratio of each effect
pp_reduced_mass = calc_reduced_mass(m_h, m_h)
CNO_reduced_mass = calc_reduced_mass(m_n, m_h)
pp_nu = calc_sensitivity_temp(c_h, c_h, np.array(Tcenlist), pp_reduced_mass)
CNO_nu = calc_sensitivity_temp(c_n, c_h, np.array(Tcenlist), CNO_reduced_mass)
weight_pp = np.array(cnolist[:4])/np.array(pplist[:4])
weight_CNO = np.array(pplist[4:])/np.array(cnolist[4:])
weight_tot = np.append(weight_pp, weight_CNO)
tot_pp_nu = (1-weight_pp)*pp_nu[:4] + weight_pp * CNO_nu[:4]
tot_CNO_nu = (1-weight_CNO)*CNO_nu[4:] + weight_CNO * pp_nu[4:]
tot_nu = np.append(tot_pp_nu, tot_CNO_nu)
#tot_nu = np.append(pp_nu[:4],CNO_nu[4:])
#%% tot_nu vs stellar mass
plt.figure(dpi  =600)
plt.scatter(masslist, tot_nu, c= 'k', linestyle = '--', label = 'pp')
plt.plot(masslist, tot_nu, c= 'k', linestyle = '--', label = 'pp')
plt.legend()
#plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'initial_mass [$M_\odot$]')
plt.ylabel(r'$\nu$', size =15)

#%% tot_nu vs Tcen
plt.figure(dpi  =600)
plt.scatter(Tcenlist, tot_nu, c= 'k', linestyle = '--', label = 'pp')
plt.plot(Tcenlist, tot_nu, c= 'k', linestyle = '--', label = 'pp')
plt.legend()
#plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$T_{cen}$ (K)')
plt.ylabel(r'$\nu$', size =15)

#%% Problem 2 ######################################################################
masslist = [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 60.0, 70.0, 100.0]
metalicity = 0.02
Rlist = []
mulist = []
mucenlist = []
for initial_mass in masslist:
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)
    Rlist.append(10**(result_tbl['logR'][0]))
    mulist.append(np.mean(result_tbl['mu'])) # mean of mu
    mucenlist.append(result_tbl['mu'][-1])
#%%
factor = np.array(mulist)**((tot_nu-4)/(tot_nu+3)) * np.array(masslist)**((tot_nu-1)/(tot_nu+3))
scaled_factor = factor/(factor[0])
scaled_radius = Rlist[0] * scaled_factor
#%%
def power_law(x, amplitude, alpha):
    y = amplitude * x ** alpha 
    return y
#%% Comparison
i = 0
plt.figure(dpi = 600)
Rerr = np.array(Rlist)*0.1
popt, pcov=  curve_fit(power_law, masslist[i:], Rlist[i:],sigma = Rerr[i:])
plt.scatter(Rlist,scaled_radius, c = 'k')
plt.plot([0,np.max(Rlist)+1],[0,np.max(Rlist)+1], linestyle = '--', c = 'k')
#plt.scatter(masslist, scaled_factor_cen)
#plt.scatter(masslist, Rlist, c= 'k', alpha = 0.3)
xs = np.arange(0.5, 100, 0.01)
#plt.plot(xs, power_law(xs, popt[0], popt[1]), c=  'k' , linestyle = '--', linewidth = 0.1)
plt.xlabel(r'$R_{model}$', size = 15)
plt.ylabel(r'$R_{Homology}$', size = 15)
plt.xscale('log')
plt.yscale('log')

#%% Comparison

i = 0
plt.figure(dpi = 600)
#plt.plot([0,np.max(Rlist)+1],[0,np.max(Rlist)+1], linestyle = '--', c = 'k')
plt.scatter(masslist, scaled_radius, label  ='Homology', c= 'b', alpha = 0.6)
plt.scatter(masslist, Rlist, c= 'r', alpha = 0.6, label = 'Model')
xs = np.arange(0.5, 100, 0.01)
#plt.plot(xs, power_law(xs, popt[0], popt[1]), c=  'k' , linestyle = '--', linewidth = 0.1)
plt.xlabel(r'Initial mass$ (M_\odot)$', size = 15)
plt.ylabel(r'Radius$ (R_\odot)$', size = 15)
plt.xscale('log')
plt.yscale('log')
plt.legend()
#%%
import matplotlib.pylab as pl
plt.figure(dpi = 600)
colors = pl.cm.jet(np.linspace(0,1,len(masslist)))
for i, initial_mass in enumerate(masslist):
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)
    plt.plot(10**result_tbl['logR']/np.max(10**result_tbl['logR']), result_tbl['mass']/np.max(result_tbl['mass']), c =colors[i], label = initial_mass)
plt.legend(ncol = 2)
plt.xlabel(r'$R_r$', size = 15)
plt.ylabel(r'$M_r$', size = 15)
#%% Problem 3 ######################################################################
#%% run TAMS model
def run_MESA(initial_mass, metalicity, run = True):
    home_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/'
    origin_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams/'
    inlist_project = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams/inlist_project'
    save_model_filename_key = f"'{initial_mass}_TAMS.mod'\n"
    filename_for_profile_when_terminate_key = f"'tams_m{initial_mass}_z{metalicity}_profile.data'\n"
    Zbase_key = f'{metalicity}d0'
    initial_mass_key = f'{initial_mass}'
    initial_z_key = f'{metalicity}d0'
    folder_key = f'tams{initial_mass}_{metalicity}/'
    
    origin_inlist = open(inlist_project, 'r')
    lines = origin_inlist.readlines()
    origin_inlist.close()
    os.system(f'cp -r {origin_dir} {home_dir}{folder_key}')
    with open(home_dir+folder_key+'inlist_project', 'w') as f:
        for line in lines:
            if 'save_model_filename' in line:
                line = line.split('=')[0] + ' = ' + save_model_filename_key
            if 'filename_for_profile_when_terminate' in line:
                line = line.split('=')[0] + ' = ' + filename_for_profile_when_terminate_key
            if 'Zbase' in line:
                line = line.split('=')[0] + ' = ' + Zbase_key
            if 'initial_mass' in line:
                line = line.split('=')[0] + ' = ' + initial_mass_key
            if 'initial_z' in line:
                line = line.split('=')[0] + ' = ' + initial_z_key
            f.write(line)
    os.chdir(home_dir+folder_key)
    if run:
        os.system('./rn')
    output_last_profile = home_dir+folder_key+f'tams_m{initial_mass}_z{metalicity}_profile.data'
    output_history = home_dir+folder_key+'LOGS/history.data'
    return output_last_profile, output_history

masslist = [0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0]

metalicity = 0.02

for initial_mass in masslist:
    last_model, all_model = run_MESA(initial_mass, metalicity, run = False)
    result_tbl = read_properties_MESA2(last_model)


# %% Load evolutionary track for each model & HR diagram data
HRdata = Table.read('/Users/hhchoi1022/Gitrepo/config/hygdata_v3-abridged.csv')

masslist = [0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0]
metalicity = 0.02


all_tmplist = []
all_lumlist = []
all_xfraclist = []
for initial_mass in masslist:
    file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
    evolutionfilelist = glob.glob(file_key)
    tot_nummodel = len(evolutionfilelist)
    lumlist = []
    tmplist = []
    xfraclist =[]
    for i in range(tot_nummodel):
        file_ = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{i+1}.data'
        result_tbl = read_properties_MESA2(file_)
        cno_key = result_tbl['cno'][-1]
        pp_key = result_tbl['pp'][-1]
        zams_key = result_tbl['x_mass_fraction_H'][-1]
        if zams_key < 0.69:
            lumlist.append(result_tbl['luminosity'][0])
            tmplist.append(10**result_tbl['logT'][0])
    all_tmplist.append(tmplist)
    all_lumlist.append(lumlist)

#%% HR diagram of TAMS
plt.figure(dpi = 600)
plt.title('HR diagram')
plt.xlim(60000, 2000)

plt.ylim(-6, 6)
plt.ylabel(r'Luminosity [$\log_{10}(L_\odot)$]')
plt.xlabel(r'$T_{eff}$ [K]')
plt.xscale('log')
plt.xticks([50000, 20000, 10000, 5000], labels=[r'$5\times 10^4$',r'$2\times 10^4$',r'$1\times 10^4$',r'$5\times 10^3$'])
plt.scatter(HRdata['T'], np.log10(10**(-(HRdata['absmag'] - 4.83)/2.5)), c = 'k', s = 1, alpha = 0.3)
for i, (tmplist, lumlist) in enumerate(zip(all_tmplist, all_lumlist)):
    plt.plot(tmplist, np.log10(lumlist) , c = 'r', linewidth = 0.5)
    plt.scatter(tmplist[0], np.log10(lumlist[0]), c = 'r', marker = '*', s = 15)
    plt.scatter(tmplist[-1], np.log10(lumlist)[-1], c = 'b', marker = '*', s = 15)
plt.scatter(0,0,marker = '*', s = 15, label = 'ZAMS', c= 'r')
plt.scatter(0,0,marker = '*', s = 15, label = 'TAMS', c= 'b')
plt.legend()


#%% problem 4 ######################################################################

initial_mass = 10.0
file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
evolutionfilelist = glob.glob(file_key)
tot_nummodel = len(evolutionfilelist)
lumlist = []
tmplist = []
xfraclist =[]
tcenlist = []
rlist = []
nuclist = []
pcenlist = []
rholist = []
mulist = []
kappalist = []
mucenlist = []
c_cenlist = []
n_cenlist = []
convlist = []
cnolist = []
pplist = []
for i in range(tot_nummodel):
    file_ = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{i+1}.data'
    result_tbl = read_properties_MESA2(file_)
    result_header = mesa.MesaData(file_).header_data
    cno_key = result_tbl['cno'][-1]
    pp_key = result_tbl['pp'][-1]
    zams_key = result_tbl['x_mass_fraction_H'][-1]
    cnolist.append(np.sum(result_tbl['cno']))
    pplist.append(np.sum(result_tbl['pp']))
    #if zams_key < 0.69:
    lumlist.append(result_tbl['luminosity'][0])
    tmplist.append(result_tbl['logT'][0])
    rlist.append(10**result_tbl['logR'][0])
    tcenlist.append(result_tbl['logT'][-1])
    nuclist.append(result_header['power_nuc_burn'])
    pcenlist.append(result_tbl['pressure'][-1])
    rholist.append(result_tbl['logRho'][-1])
    mulist.append(np.mean(result_tbl['mu']))
    mucenlist.append(result_tbl['mu'][-1])
    kappalist.append(np.mean(result_tbl['opacity']))
    c_cenlist.append(result_tbl['c12'][-1])
    n_cenlist.append(result_tbl['n14'][-1])
    convlist.append(result_tbl['lum_conv_div_lum_rad'][-1])
    #plt.figure()
    #plt.plot(result_tbl['mass']/np.max(result_tbl['mass']),result_tbl['lum_conv_div_L'])

#%%
plt.figure(dpi = 600)
ax = plt.subplot()
ax.set_ylim(-0.1, 1.1)
ax.set_xlim(35, 60)
ax.set_xlabel('Phase')
ax.set_ylabel('Relative scale')
ax.plot(np.array(cnolist)/np.max(cnolist), c='g', label = 'CNO')
ax.plot(np.array(pplist)/np.max(pplist), c='b', label = 'pp')
ax.plot([0], c='r', label = '[N14/C12]')
ax1 = ax.twinx()
ax1.plot(np.array(n_cenlist)/np.array(c_cenlist), c = 'r')
ax1.set_ylabel('[N14/C12]')
ax.legend(loc = 4)


#%% Problem 5 ######################################################################
masslist = [0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0]
tams_ncratiolist = []
tams_ncratiolist_surf = []
zams_ncratiolist = []
zams_ncratiolist_surf = []
nlist = []
clist = []
zams_ncenlist = []
zams_ccenlist = []
tams_ncenlist = []
tams_ccenlist = []
for initial_mass in masslist:
    file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
    evolutionfilelist = glob.glob(file_key)
    tot_nummodel = len(evolutionfilelist)
    zams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile1.data'
    tams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{tot_nummodel}.data'
    zams_tbl = read_properties_MESA2(zams_file)
    tams_tbl = read_properties_MESA2(tams_file)
    zams_c_cen = zams_tbl['c12'][-1]
    zams_n_cen = zams_tbl['n14'][-1]
    zams_nc_ratio = round(zams_n_cen/zams_c_cen,5)
    zams_ncratiolist.append(zams_nc_ratio)
    zams_ncenlist.append(zams_n_cen)
    zams_ccenlist.append(zams_c_cen)
    tams_c_cen = tams_tbl['c12'][-1]
    tams_n_cen = tams_tbl['n14'][-1]
    tams_nc_ratio = round(tams_n_cen/tams_c_cen,5)
    tams_ncratiolist.append(tams_nc_ratio)
    tams_ncenlist.append(tams_n_cen)
    tams_ccenlist.append(tams_c_cen)
#%%
plt.figure(figsize = (4,3), dpi  =600)
plt.title('ZAMS')
plt.xlabel(r'stellar mass [$M_\odot$]')
plt.ylabel(r'[N/C]')
plt.scatter(masslist, zams_ncratiolist, edgecolors = 'k', facecolors = 'none')

plt.figure(figsize = (4,3), dpi  =600)
plt.title('TAMS')
plt.xlabel(r'stellar mass [$M_\odot$]')
plt.ylabel(r'[N/C]')
plt.scatter(masslist, tams_ncratiolist, edgecolors = 'k', facecolors = 'none')

#%%
plt.figure(figsize = (4,3), dpi  =600)
plt.title('ZAMS')
plt.xlabel(r'stellar mass [$M_\odot$]')
plt.ylabel(r'[N/C]')
ax = plt.subplot()
ax.scatter(masslist, zams_ncenlist, edgecolors = 'g', facecolors = 'none')
ax1= ax.twinx()
ax1.scatter(masslist, zams_ccenlist, edgecolors = 'r', facecolors = 'none')

plt.figure(figsize = (4,3), dpi  =600)
plt.title('TAMS')
plt.xlabel(r'stellar mass [$M_\odot$]')
plt.ylabel(r'fraction')
ax = plt.subplot()
ax.scatter(masslist, tams_ncenlist, edgecolors = 'g', facecolors = 'none', label  ='Nitrogen')
ax1 = ax.twinx()
ax1.scatter(0, 0, edgecolors = 'g', facecolors = 'none', label  ='Nitrogen')
ax1.scatter(masslist, tams_ccenlist, edgecolors = 'r', facecolors = 'none', label = 'Carbon')
plt.legend(loc = 4)
# %%
plt.scatter(masslist, )
#%% Convection core & carbon distribution 
import matplotlib.pylab as pl
colors = pl.cm.jet(np.linspace(0,1,11))
plt.figure(dpi =600, figsize = (6,4))
plt.title('TAMS')
plt.xlabel(r'$M_r$')
plt.ylabel(r'$\frac{L_{conv}}{L_{tot}}$', size = 15)
for i, initial_mass in enumerate(masslist):
    
    file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
    evolutionfilelist = glob.glob(file_key)
    tot_nummodel = len(evolutionfilelist)
    tams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{tot_nummodel}.data'
    result_tbl = read_properties_MESA2(tams_file)
    ax = plt.subplot()
    ax.plot(result_tbl['mass']/np.max(result_tbl['mass']),result_tbl['c12']/np.max(result_tbl['c12']), c=colors[i], alpha = 0.2, linestyle = '--')
    ax.plot(result_tbl['mass']/np.max(result_tbl['mass']),result_tbl['lum_conv_div_L'], c=colors[i], label = rf'{initial_mass}$M_\odot$')
    ax.set_xlim(0,1)
    ax1 = ax.twinx()
    ax1.set_ylabel(r'$frac_{scaled}$ [C12]')
    ax1.plot(result_tbl['mass']/np.max(result_tbl['mass']),result_tbl['c12']/np.max(result_tbl['c12']), c=colors[i], alpha = 0.2, linestyle = '--')
    #ax1.plot(result_tbl['mass']/np.max(result_tbl['mass']), result_tbl['cno'], c = 'k')
ax.legend()

#%% problem 6 ######################################################################

initial_mass = 10.0
file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
evolutionfilelist = glob.glob(file_key)
tot_nummodel = len(evolutionfilelist)
lumlist = []
tmplist = []
xfraclist =[]
tcenlist = []
rlist = []
nuclist = []
pcenlist = []
rholist = []
mulist = []
kappalist = []
mucenlist = []
c_cenlist = []
n_cenlist = []
convlist = []
cnolist = []
pplist = []
for i in range(tot_nummodel):
    file_ = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{i+1}.data'
    result_tbl = read_properties_MESA2(file_)
    result_header = mesa.MesaData(file_).header_data
    cno_key = result_tbl['cno'][-1]
    pp_key = result_tbl['pp'][-1]
    zams_key = result_tbl['x_mass_fraction_H'][-1]
    cnolist.append(np.sum(result_tbl['cno']))
    pplist.append(np.sum(result_tbl['pp']))
    if zams_key < 0.69:
        lumlist.append(result_tbl['luminosity'][0])
        tmplist.append(result_tbl['logT'][0])
        rlist.append(10**result_tbl['logR'][0])
        tcenlist.append(result_tbl['logT'][-1])
        nuclist.append(result_header['power_nuc_burn'])
        pcenlist.append(result_tbl['pressure'][-1])
        rholist.append(result_tbl['logRho'][-1])
        mulist.append(np.mean(result_tbl['mu']))
        mucenlist.append(result_tbl['mu'][-1])
        kappalist.append(np.mean(result_tbl['opacity']))
        c_cenlist.append(result_tbl['c12'][-1])
        n_cenlist.append(result_tbl['n14'][-1])
        convlist.append(result_tbl['lum_conv_div_lum_rad'][-1])
    #plt.figure()
    #plt.plot(result_tbl['mass']/np.max(result_tbl['mass']),result_tbl['lum_conv_div_L'])

#%% HR diagram of 10 solar mass
plt.figure(dpi = 600)
plt.title('HR diagram')

plt.plot(10**np.array(tmplist), np.log10(lumlist) , c = 'r', linewidth = 0.5)
plt.scatter(10**tmplist[0], np.log10(lumlist[0]), c = 'r', marker = '*', s = 15)
plt.scatter(10**tmplist[-1], np.log10(lumlist)[-1], c = 'b', marker = '*', s = 15)
plt.scatter(0,0,marker = '*', s = 15, label = 'ZAMS', c= 'r')
plt.scatter(0,0,marker = '*', s = 15, label = 'TAMS', c= 'b')

plt.xlim(30000, 15000)
plt.ylim(2, 6)
plt.ylabel(r'Luminosity [$\log_{10}(L_\odot)$]')
plt.xlabel(r'$T_{eff}$ [K]')
plt.legend()
# %% mu - (Pcen, Tcen, eta_nuc, Radius)

plt.figure(dpi = 600)
ax1 = plt.subplot()
ax1.plot(mucenlist, pcenlist/np.max(pcenlist), c='b', label = r'$P_{cen}$')
ax1.plot(mucenlist, 10**np.array(tcenlist)/np.max(10**np.array(tcenlist)), c = 'g', label = r'$T_{cen}$')
ax1.plot(mucenlist, nuclist/np.max(nuclist), c = 'r', label = r'$\eta_{nuc}$')
ax1.plot(mucenlist, rlist/np.max(rlist), c = 'k', label = 'Radius')
plt.legend()
plt.xlabel(r'$\mu_{cen}$')
plt.ylabel('Relative scale')
#%% rho_cen - T_cen to check degeneracy
plt.figure(dpi = 300)
plt.ylim(6,10)
plt.xlim(0,3)
plt.scatter(rholist, tcenlist, c =(np.arange(len(tcenlist))/len(tcenlist)))
plt.xlabel(r'log$\rho_c (g/cm^{-3})$')
plt.ylabel(r'$log T_c$ (K)')
plt.colorbar()
# %% Homology prediction with const kappa
def calc_hom_const_factor(kappa, mu, mass):
    factor = 1/kappa*mu**4*mass**3
    return factor
# %% Homology prediction with Kramers kappa
def calc_hom_Kramers_factor(mu, mass, radius):
    factor = mu**7.5*mass**5.5/radius**0.5
    return factor

# %%
kap_ratio = np.array(kappalist)/kappalist[0]
mu_ratio = np.array(mulist)/mulist[0]
r_ratio = np.array(rlist)/rlist[0]

const_factor = calc_hom_const_factor(kap_ratio,mu_ratio,1)
karmers_factor = calc_hom_Kramers_factor(mu_ratio,1,r_ratio)
hom_const = lumlist[0] * const_factor
hom_kramers = lumlist[0] * karmers_factor
# %% L_model vs L_homology
plt.figure(dpi = 600)

plt.scatter(np.log10(np.array(lumlist)), np.log10(np.array(hom_const)), c= 'k', label = r'$\overline{\kappa} = const$')
plt.scatter(np.log10(np.array(lumlist)), np.log10(np.array(hom_kramers)), c = 'g', label = r'$\overline{\kappa} \propto \rho T^{-7/2}$')
plt.plot([3.5,4.2],[3.5,4.2], linewidth = 1 ,linestyle = '--', c= 'k')
plt.scatter(np.log10(np.array(lumlist)[0]), np.log10(np.array(hom_const)[0]), s = 100, marker = '*', c= 'r', label = r'ZAMS')
plt.scatter(np.log10(np.array(lumlist)[-1]), np.log10(np.array(hom_const)[-1]), s = 100, marker = '*', c= 'b')
plt.scatter(np.log10(np.array(lumlist)[-1]), np.log10(np.array(hom_kramers)[-1]), s = 100, marker = '*', c= 'b', label = r'TAMS')

plt.xlim(3.7,4.1)
plt.ylim(3.7,4.3)
plt.xlabel(r'$\log(L_{model})$')
plt.ylabel(r'$\log(L_{Homology})$')
plt.legend()

#%% problem 7 ######################################################################
#%% constant (cgs unit)
const_k = 1.3807e-16 # [cm2 g s-2 K-1]
const_m_u = 1.6605e-24 #[g]
const_a = 7.5646e-15 # [erg cm-3 K-4]
const_msun = 1.9884e+33 # [g]
const_rsun = 6.9570e10 # [cm]
const_lsun = 3.8280e33 # [erg s-1]
const_g = 6.6743e-8 #[cm-3 g-1 s-2]
const_c = 2.9979e10 #[cm s-1]
const_Kn = 1.0035e13
const_Ke = 1.2435e15

masslist = [1.0, 8.0]
metalicity = 0.02
initial_mass = 1.0
file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
evolutionfilelist = glob.glob(file_key)
tot_nummodel = len(evolutionfilelist)
zams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile1.data'
tams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{tot_nummodel}.data'
zams_tbl = read_properties_MESA2(zams_file)
tams_tbl = read_properties_MESA2(tams_file)
zams_header = mesa.MesaData(zams_file).header_data
tams_header = mesa.MesaData(tams_file).header_data
#%% Rest mass energy difference
hydrogen_difference = zams_header['star_mass_h1']-tams_header['star_mass_h1']
mass_difference_onetime = 0.02866 * const_m_u
energy_onetime = mass_difference_onetime * const_c **2
n_fusion = (hydrogen_difference * const_msun / (1.008 * const_m_u))/4
dif_rest_energy = n_fusion * energy_onetime 
# Radiative energy consumed
age = 0
rad_energy = 0
tot_nummodel = len(evolutionfilelist)
for i in range(tot_nummodel):
    file_ = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{i+1}.data'
    tbl = read_properties_MESA2(file_)
    header = mesa.MesaData(file_).header_data
    duration = (header['star_age'] - age)*86400*365
    luminosity = header['photosphere_L'] * const_lsun
    rad_energy += luminosity * duration
    age = header['star_age']
    #print(age, header['photosphere_L'])
print(dif_rest_energy)
print(rad_energy)
# %% Problem 9 ######################################################################
initial_mass = 10.0
metalicity = 0.02

file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile*.data'
evolutionfilelist = glob.glob(file_key)
tot_nummodel = len(evolutionfilelist)
zams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile1.data'
tams_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/tams{initial_mass}_{metalicity}/LOGS/profile{tot_nummodel}.data'
zams_tbl = read_properties_MESA2(zams_file)
tams_tbl = read_properties_MESA2(tams_file)
zams_header = mesa.MesaData(zams_file).header_data
tams_header = mesa.MesaData(tams_file).header_data


file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/he_exhaustion/LOGS/profile*.data'
evolutionfilelist = glob.glob(file_key)
tot_nummodel = len(evolutionfilelist)
he_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/he_exhaustion/LOGS/profile{tot_nummodel}.data'
he_exhaustion_tbl = read_properties_MESA2(he_file)
he_exhaustion_header = mesa.MesaData(he_file).header_data
templist = []
lumlist = []
for i in range(tot_nummodel):
    he_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/he_exhaustion/LOGS/profile{i+1}.data'
    result_tbl = read_properties_MESA2(he_file)
    zams_key = result_tbl['x_mass_fraction_H'][-1]
    if zams_key < 0.699:
        templist.append(10**result_tbl['logT'][0])
        lumlist.append(result_tbl['luminosity'][0])
#%%
he_burning_timescale = he_exhaustion_header['star_age'] - tams_header['star_age']
h_burning_timescale = tams_header['star_age'] - zams_header['star_age']
print(h_burning_timescale/he_burning_timescale)
# %%
#%% HR diagram of TAMS
plt.figure(dpi = 600)
plt.title('HR diagram')
plt.xlim(30000, 2000)

plt.ylim(-6, 6)
plt.ylabel(r'Luminosity [$\log_{10}(L_\odot)$]')
plt.xlabel(r'$T_{eff}$ [K]')
plt.xscale('log')
plt.xticks([50000, 20000, 10000, 5000], labels=[r'$5\times 10^4$',r'$2\times 10^4$',r'$1\times 10^4$',r'$5\times 10^3$'])
plt.scatter(HRdata['T'], np.log10(10**(-(HRdata['absmag'] - 4.83)/2.5)), c = 'k', s = 1, alpha = 0.3)
plt.plot(templist, np.log10(lumlist) , c = 'r', linewidth = 0.5)
plt.scatter(templist[0], np.log10(lumlist[0]), c = 'r', marker = '*', s = 15)
plt.scatter(templist[-1], np.log10(lumlist)[-1], c = 'b', marker = '*', s = 15)
plt.scatter(0,0,marker = '*', s = 15, label = 'TAMS', c= 'r')
plt.scatter(0,0,marker = '*', s = 15, label = 'He exhaustion', c= 'b')
plt.legend()

# %%

#%%
######################################################################
#################            HOMEWORK 4            ###################
######################################################################
#%% Problem A
def run_MESA(initial_mass, metalicity, run = True):
    home_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/'
    origin_dir = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/trgb/'
    inlist_project = '/Users/hhchoi1022/Program/mesa-r22.05.1/star/trgb/inlist_project'
    save_model_filename_key = f"'{initial_mass}_TRGB.mod'\n"
    filename_for_profile_when_terminate_key = f"'tams_m{initial_mass}_z{metalicity}_profile.data'\n"
    Zbase_key = f'{metalicity}d0'
    initial_mass_key = f'{initial_mass}'
    initial_z_key = f'{metalicity}d0'
    folder_key = f'trgb{initial_mass}_{metalicity}/'
    
    origin_inlist = open(inlist_project, 'r')
    lines = origin_inlist.readlines()
    origin_inlist.close()
    os.system(f'cp -r {origin_dir} {home_dir}{folder_key}')
    with open(home_dir+folder_key+'inlist_project', 'w') as f:
        for line in lines:
            if 'save_model_filename' in line:
                line = line.split('=')[0] + ' = ' + save_model_filename_key
            if 'filename_for_profile_when_terminate' in line:
                line = line.split('=')[0] + ' = ' + filename_for_profile_when_terminate_key
            if 'Zbase' in line:
                line = line.split('=')[0] + ' = ' + Zbase_key
            if 'initial_mass' in line:
                line = line.split('=')[0] + ' = ' + initial_mass_key + '\n'
            if 'initial_z' in line:
                line = line.split('=')[0] + ' = ' + initial_z_key
            f.write(line)
    os.chdir(home_dir+folder_key)
    if run:
        os.system('./rn')
    output_last_profile = home_dir+folder_key+f'tams_m{initial_mass}_z{metalicity}_profile.data'
    output_history = home_dir+folder_key+'LOGS/history.data'
    return output_last_profile, output_history
masslist = [0.9, 1.2, 1.7, 2.0]
filelist = []
for mass in masslist:
    filelist.append(run_MESA(mass, 0.001, run = False)[0])
#%% Problem A.1)
def read_properties_MESA2(output):
    h_last = mesa.MesaData(output)
    #h_last.header_names # header data indicate the resultant physical property at the end of the model 
    #h_last.bulk_names # bulk data indicate the radial profile of each physial properties 
    
    result_tbl = Table()
    '''
    m = h_last.mass*const_msun 
    r = 10**h_last.logR * const_rsun
    t = 10**h_last.logT
    rho = 10**h_last.logRho
    p = 10**h_last.logP
    h = h_last.x_mass_fraction_H
    he = h_last.y_mass_fraction_He
    metals = h_last.z_mass_fraction_metals
    pressure_tot = h_last.pressure
    pressure_rad = h_last.prad
    pressure_gas = h_last.pgas
    mu = h_last.mu
    opacity = h_last.opacity
    luminosity = h_last.luminosity
    energy_total = h_last.total_energy
    energy_intern = h_last.energy
    energy_cellIE = h_last.cell_specific_IE
    energy_cellKE = h_last.cell_specific_KE
    frac_KEIE = h_last.cell_internal_energy_fraction
    '''
    colnames = h_last.bulk_names
    for colname in colnames:
        data = h_last.bulk_data[colname]
        result_tbl[colname] = data
    return result_tbl
#%%Problem A.1)
coreHelist = []
lumlist = []
for file in filelist:
    tbl = read_properties_MESA2(file)
    hdr = mesa.MesaData(file).header_data
    coreHelist.append(hdr['he_core_mass'])
    lumlist.append(hdr['photosphere_L'])
plt.figure(dpi = 600)
ax1 = plt.subplot()
ax2 = ax1.twinx()
ax2.sharex(ax1)
ax1.scatter(masslist, coreHelist, c = 'b', marker = 'o')
ax1.set_ylim(np.min(coreHelist)-0.01, np.max(coreHelist)+0.01)
ax1.set_ylabel(r'M$_{core}$($M_\odot$)')
ax2.scatter(masslist, lumlist, c = 'b', marker = 'o', label = r'$M_{core}$')
ax2.scatter(masslist, lumlist, c = 'r', marker = 'o', label = r'$Luminosity$')
ax1.set_xlabel(r'$M_{init}(M_\odot$)')
ax2.set_ylim(np.min(lumlist)-300, np.max(lumlist)+300)
ax2.set_ylabel(r'Luminosity(${L/L_\odot})$')
ax2.legend()

#%% Problem A.2)
from scipy.optimize import curve_fit
def func(x, a, b):
    return a*x + b
popt, pcov = curve_fit(func, xdata = coreHelist, ydata = lumlist)
xrange_ = np.linspace(0, 1, 100)
plt.figure(dpi = 600)
ax1 = plt.subplot()

ax1.scatter(coreHelist, lumlist, c = 'k', marker = 'o')
ax1.plot(xrange_, func(xrange_, popt[0], popt[1]), linestyle = '--',  c= 'k', label = r'')
ax1.set_xlim(np.min(coreHelist)-0.01, np.max(coreHelist)+0.01)
ax1.set_ylim(np.min(lumlist)-300, np.max(lumlist)+300)
ax1.set_xlabel(r'M$_{core}$($M_\odot$)')
ax1.set_ylabel(r'Luminosity(${L/L_\odot})$')

#%% Problem B.2.a,b)
def mass_AGB(time : float or np.array, 
             M0 : float = 2,
             eta : float = 3):
    """
    parameters
        1. time : time since zero age AGB star (yr)
    returns 
        1. mass : AGB star mass with the consideration of stellar wind (M_\odot)
    """
    return (-0.34*eta/3 * (np.exp( (8.85*time*10**-7) ) -1) + M0**2)**0.5

def masscore_AGB(time : float or np.array):
    """
    parameters
        1. time : time since zero age AGB star (yr)
    returns 
        1. mass : core mass of a AGB star (M_\odot)
    """
    return ( 1 / (5.9*10) * np.exp(5.9*time*10**-7) + 0.52)

eta = 9
t_coordinate = np.arange(0, 1e7, 1000)
AGB_mass = mass_AGB(t_coordinate, eta = eta)
lifetime = np.sum(~np.isnan(AGB_mass))
t_coordinate = np.arange(0, 1000*lifetime, 1000)
AGB_mass = mass_AGB(t_coordinate, eta = eta)
core_mass = masscore_AGB(t_coordinate)
TAGB_idx = np.argmin(np.abs(AGB_mass-core_mass))
TAGB_age = t_coordinate[TAGB_idx]
#%% Problem B.2.c)
t_coordinate_9 = t_coordinate
AGB_mass_9 = AGB_mass
core_mass_9 = core_mass
TAGB_age_9 = TAGB_age
plt.figure(dpi = 600)
plt.title(r'$2M_\odot AGB star$')
plt.plot(t_coordinate_1, AGB_mass_1, c = 'k', label = r'$M_{AGB}(\eta$ = 1)', linestyle = ':')
plt.plot(t_coordinate_1, core_mass_1, c ='b', label  = r'$M_{core}(\eta$ = 1)', linestyle = ':')
plt.axvline(TAGB_age_1, c = 'r', linestyle = ':')
plt.plot(t_coordinate_3, AGB_mass_3, c = 'k', label = r'$M_{AGB}(\eta$ = 3)', linestyle = '--')
plt.plot(t_coordinate_3, core_mass_3, c ='b', label  = r'$M_{core}(\eta$ = 3)', linestyle = '--')
plt.axvline(TAGB_age_3, c = 'r', linestyle = '--')
plt.plot(t_coordinate_9, AGB_mass_9, c = 'k', label = r'$M_{AGB}(\eta$ = 9)', linestyle = '-')
plt.plot(t_coordinate_9, core_mass_9, c ='b', label  =r'$M_{core}(\eta$ = 9)', linestyle = '-')
plt.axvline(TAGB_age_9, c = 'r', linestyle = '-')
plt.xlabel('time (yr)')
plt.ylabel(r'mass $(M_\odot)$')
plt.legend()
#%% Problem B.2.b)
def lum_from_t(time : float or np.array):
    """
    parameters
        1. time : time since zero age AGB star (yr)
    returns 
        1. luminosity : solar luminosity unit (L_\odot)
    """
    return 10**3 * np.exp(5.9*time*10**-7)
lum_from_t(TAGB_age)
#%% Problem B.2.b)
mass_AGB(TAGB_age)
# %% Problem B.3)
import numpy as np
def lum_from_mass(initial_mass : float or np.array):
    """
    parameters
        1. initial mass : solar mass unit (M_\odot)
    returns 
        1. luminosity : solar luminosity unit (L_\odot)
    """
    return 10**(0.781 + 2.760 * np.log10(initial_mass))
def wind_from_lum(luminosity : float or np.array):
    """
    parameters
        1. luminosity : solar luminosity unit (L_\odot)
    returns 
        1. mass loss rate : mass loss rate (M_\odot/yr)
    """
    return 10**(-12.76 + 1.3 * np.log10(luminosity))
def timescale_from_mass(initial_mass : float or np.array):
    """
    parameters
        1. initial mass : solar mass unit (M_\odot)
    returns 
        1. timescale : duration of the main sequence phase (yr)
    """
    return 10**(7.719 - 0.655 * np.log10(initial_mass))
#%% B.3.a)
masslist = np.array([25,40,60,85,120])
timescalelist = timescale_from_mass(masslist)
lumlist = lum_from_mass(masslist)
windlosslist = wind_from_lum(lumlist)
MS_masslosslist = windlosslist * timescalelist
fractionlist = MS_masslosslist / masslist
from astropy.table import Table
resultB3 = Table(names = ['mass','timescale','luminosity','massloss','massloss_frac'], data = [masslist, timescalelist.round(1),lumlist.round(0), windlosslist, fractionlist.round(2)])
resultB3
#%% B.3.b)
# For 85 solar mass star, massloss rate
initial_mass = 85
result_85 = resultB3[resultB3['mass'] == initial_mass]
timescale_exposecore = (1- 0.83) * result_85['mass'] /result_85['massloss']
print("Timescale = %dyr "%(timescale_exposecore))
#%% Problem C.1)
q = 15/1.4
factor = (0.49*q**(2/3)) / (0.6*q**(2/3) + np.log(1+q**(1/3)))
A = 1.07
roche_R = factor * A
au = 1.496e13 #cm
R_L = roche_R * au
R_L
# %% Problem C.2)

templist = []
lumlist = []
rlist = []
file_key = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/trgb15.0_0.02/LOGS/profile*.data'
evolutionfilelist = glob.glob(file_key)
tot_nummodel = len(evolutionfilelist)
for i in range(tot_nummodel):
    file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/trgb15.0_0.02/LOGS/profile{i+1}.data'
    tbl = read_properties_MESA2(file)
    header = mesa.MesaData(file).header_data
    templist.append(10**tbl['logT'][0])
    lumlist.append(tbl['luminosity'][0])
    rlist.append(header['photosphere_r'])
#%%
zams_idx  = 21
fill_idx = 50
max_idx = 100
templist_zams = templist[zams_idx:max_idx]
lumlist_zams = lumlist[zams_idx:max_idx]
rlist_zams = lumlist[zams_idx:max_idx]
plt.figure(dpi = 300)
plt.plot(np.log10(templist_zams), np.log10(lumlist_zams), c = 'k')
plt.scatter(np.log10(templist[fill_idx]), np.log10(lumlist[fill_idx]), c = 'r', marker = '*', s = 100, label = 'fill_R_L')
plt.legend()
plt.xlim(4.8, 3.4)
plt.ylim(4.0, 6.5)
plt.xlabel(r'$\log(T_{eff})$')
plt.ylabel(r'$\log(L/L_{\odot})}$')
#%%
fill_file = f'/Users/hhchoi1022/Program/mesa-r22.05.1/star/trgb15.0_0.02/LOGS/profile51.data'
fill_tbl = read_properties_MESA2(fill_file)
fill_header = mesa.MesaData(fill_file).header_data

# %%
fill_tbl
# %%
print(fill_header)
# %%
plt.figure(dpi = 300)
plt.plot(fill_tbl['mass'],fill_tbl['cno']/np.max(fill_tbl['cno']), label = 'Hydrogen burning')
plt.plot(fill_tbl['mass'],fill_tbl['tri_alfa']/np.max(fill_tbl['cno']), label = 'He burning')
plt.text(5, 0.8, r'Core temperature = $1.98\times 10^8K$')
plt.text(5, 0.7, r'He Core mass = $3.87M_{\odot}$')
plt.text(5, 0.6, r'CO Core mass = $0.00M_{\odot}$')
#plt.plot(tbl['mass'],tbl['c_alpha']/np.max(tbl['cno']), label = 'Carbon burning')
#plt.plot(tbl['mass'],tbl['c_alpha']/np.max(tbl['cno']), label = 'Carbon burning')
plt.legend()
plt.xlabel(r'Mass coordinate $(M_\odot)$')
plt.ylabel(r'Relative reaction rate')
#%%
env_tbl = fill_tbl[fill_tbl['mass'] > fill_header['he_core_mass']]
core_tbl = fill_tbl[fill_tbl['mass'] < fill_header['he_core_mass']]
# %%
msun = header['msun']
rsun = header['rsun']
lsun = header['lsun']
G = 6.6743e-8 # cgs
#%% binding energy
R = 10**env_tbl['logR'] * rsun
M_r = env_tbl['mass'] * msun
dM_r = env_tbl['dm']
grav_potential = G*M_r/R * dM_r
internal = env_tbl['dm'] * env_tbl['energy']
binding = internal - grav_potential
tot_binding = np.sum(binding)
# %% lambda
env_mass = (env_tbl[0]['mass'] - env_tbl[-1]['mass']) * msun
donor_mass = tbl[0]['mass'] * msun
lambda_ = ((0-tot_binding) * (G * donor_mass * env_mass / R_L)**-1)**-1
#%%
M_a = 1.4 * msun
M_cd = core_tbl['mass'][0] * msun
M_d = donor_mass
a_f = 9.28e12 #cm


# %%
G*M_a*M_cd/2
# %%
G*M_a*M_d/2/(A*au)
# %%
tot_binding
# %%

