#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 09:25:31 2022

@author: hhchoi1022
"""
#%%
import numpy as np
import extinction
import pyphot
import glob, os
from astropy.table import Table
from astropy.table import vstack
import re
from astropy.io import ascii
from convert_AB_Vega import ABVegaMagnitude
#%%
def correct_host_extinction(tbl, ebv, Rv = 3.1, filename = None, dereddening = True): 
    #ebv = 0.097 # Hosseinzadeh 2022
    #Rv = 3.1
    #Av = ebv * Rv # R_v = 3.1
    corrected_tbl = tbl.copy()
    lib = pyphot.get_library() # Filter library 
    filter_key = dict(u = 'SDSS_u',
                      g = 'SDSS_g',
                      r = 'SDSS_r',
                      i = 'SDSS_i',
                      U = 'GROUND_JOHNSON_U',
                      B = 'GROUND_JOHNSON_B',
                      V = 'GROUND_JOHNSON_V',
                      R = 'GROUND_COUSINS_R')
    filterlist = sorted(list(set(tbl['filter'])))
    corrected_tbl['[mag_cor-mag_origin]'] = -99.9
    Av = Rv * ebv
    for filtname in filterlist:
        if filtname not in filter_key.keys():
            print(f'filter "{filtname}" is not in filter_key')
            pass
        else:
            filter_ = lib[filter_key[filtname]]
            lpivot_AA = np.array([filter_.lpivot.value])
            Afilt = round(extinction.fitzpatrick99(lpivot_AA, Av, Rv)[0],3)
            if dereddening == True:
                corrected_tbl['mag'][corrected_tbl['filter'] == filtname] -= Afilt
                corrected_tbl['[mag_cor-mag_origin]'][corrected_tbl['filter'] == filtname] = -Afilt
            
            elif dereddening == False:
                corrected_tbl['mag'][corrected_tbl['filter'] == filtname] += Afilt
                corrected_tbl['[mag_cor-mag_origin]'][corrected_tbl['filter'] == filtname] = Afilt
            else:
                pass
    corrected_tbl['mag'] = corrected_tbl['mag'].round(3)
    corrected_tbl['[mag_cor-mag_origin]'] = corrected_tbl['[mag_cor-mag_origin]'].round(3)
    if not filename == None:
        corrected_tbl.write(f'{filename}', format = 'ascii.fixed_width', overwrite = True)
    return corrected_tbl

#%%
def correct_mw_extinction(tbl, mwebv = None, ra = None, dec = None, mwRv = 3.1, filename = None, dereddening = True): 
    #ebv = 0.097 # Hosseinzadeh 2022
    #Rv = 3.1
    #Av = ebv * Rv # R_v = 3.1
    import sfdmap
    dustmap = sfdmap.SFDMap("/Users/hhchoi1022/Gitrepo/config/sfddata-master")
    if ra is not None:
        try:
            mwebv = dustmap.ebv(ra, dec)
        except:
            raise ValueError('cannot calculate extinction with given coordinates')
    #RA = 64.9708333, Dec = -54.9480556 for SN2021aefx
    corrected_tbl = tbl.copy()
    
    lib = pyphot.get_library() # Filter library 
    filter_key = dict(u = 'SDSS_u',
                      g = 'SDSS_g',
                      r = 'SDSS_r',
                      i = 'SDSS_i',
                      U = 'GROUND_JOHNSON_U',
                      B = 'GROUND_JOHNSON_B',
                      V = 'GROUND_JOHNSON_V',
                      R = 'GROUND_COUSINS_R')
    filterlist = sorted(list(set(tbl['filter'])))
    corrected_tbl['[mag_cor-mag_origin]'] = -99.9
    Av = mwRv * mwebv
    for filtname in filterlist:
        if filtname not in filter_key.keys():
            print(f'filter "{filtname}" is not in filter_key')
            pass
        else:
            filter_ = lib[filter_key[filtname]]
            lpivot_AA = np.array([filter_.lpivot.value])
            Afilt = round(extinction.fitzpatrick99(lpivot_AA, Av, mwRv)[0],3)
            if dereddening == True:
                corrected_tbl['mag'][corrected_tbl['filter'] == filtname] -= Afilt
                corrected_tbl['[mag_cor-mag_origin]'][corrected_tbl['filter'] == filtname] = -Afilt
            elif dereddening == False:
                corrected_tbl['mag'][corrected_tbl['filter'] == filtname] += Afilt
                corrected_tbl['[mag_cor-mag_origin]'][corrected_tbl['filter'] == filtname] = Afilt
            else:
                pass
    corrected_tbl['mag'] = corrected_tbl['mag'].round(3)
    corrected_tbl['[mag_cor-mag_origin]'] = corrected_tbl['[mag_cor-mag_origin]'].round(3)
    if not filename == None:
        corrected_tbl.write(f'{filename}', format = 'ascii.fixed_width', overwrite = True)
    return corrected_tbl
#%%4
def correct_both_extinction(tbl, hostebv, hostRv, mwebv = None, ra = None, dec = None , mwRv = 3.1, filename = None, dereddening = True): 
    #ebv = 0.097 # Hosseinzadeh 2022
    #Rv = 3.1
    #Av = ebv * Rv # R_v = 3.1
    import sfdmap
    dustmap = sfdmap.SFDMap("/Users/hhchoi1022/Gitrepo/config/sfddata-master")
    if ra is not None:
        try:
            mwebv = dustmap.ebv(ra, dec)
        except:
            raise ValueError('cannot calculate extinction with given coordinates')
    #RA = 64.9708333, Dec = -54.9480556 for SN2021aefx
    corrected_tbl = tbl.copy()
    
    lib = pyphot.get_library() # Filter library 
    filter_key = dict(u = 'SDSS_u',
                      g = 'SDSS_g',
                      r = 'SDSS_r',
                      i = 'SDSS_i',
                      U = 'GROUND_JOHNSON_U',
                      B = 'GROUND_JOHNSON_B',
                      V = 'GROUND_JOHNSON_V',
                      R = 'GROUND_COUSINS_R')
    filterlist = sorted(list(set(tbl['filter'])))
    corrected_tbl['[mag_cor-mag_origin]'] = -99.9
    mwAv = mwRv * mwebv
    hostAv = hostRv * hostebv
    for filtname in filterlist:
        if filtname not in filter_key.keys():
            print(f'filter "{filtname}" is not in filter_key')
        else:
            filter_ = lib[filter_key[filtname]]
            lpivot_AA = np.array([filter_.lpivot.value])
            mw_Afilt = round(extinction.fitzpatrick99(lpivot_AA, mwAv, mwRv)[0],3)
            host_Afilt = round(extinction.fitzpatrick99(lpivot_AA, hostAv, hostRv)[0],3)
            tot_Afilt = mw_Afilt + host_Afilt
            if dereddening == True:
                corrected_tbl['mag'][corrected_tbl['filter'] == filtname] -= tot_Afilt
                corrected_tbl['[mag_cor-mag_origin]'][corrected_tbl['filter'] == filtname] = -tot_Afilt
            elif dereddening == False:
                corrected_tbl['mag'][corrected_tbl['filter'] == filtname] += tot_Afilt
                corrected_tbl['[mag_cor-mag_origin]'][corrected_tbl['filter'] == filtname] = tot_Afilt
            else:
                pass
    corrected_tbl['mag'] = corrected_tbl['mag'].round(3)
    corrected_tbl['[mag_cor-mag_origin]'] = corrected_tbl['[mag_cor-mag_origin]'].round(3)
    if not filename == None:
        corrected_tbl.write(f'{filename}', format = 'ascii.fixed_width', overwrite = True)
    return corrected_tbl

#%%
ebv = 0.097  # Hosseinzadeh 2022
ra = 64.9708333
dec = -54.9480556
# data from Ahsall 2022 is already corrected for MW extinction >>> noextin_dat
if __name__ == '__main__':
    IMSNG_file = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/IMSNG/original/IMSNG_all'
    IMSNG_savepath = os.path.dirname(os.path.dirname(IMSNG_file))
    IMSNG_tbl = ascii.read(IMSNG_file, format ='fixed_width')
    H22_file = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/original/Hosseinzadeh2022_all'
    H22_savepath = os.path.dirname(os.path.dirname(H22_file))
    H22_tbl_raw = ascii.read(H22_file, format ='fixed_width')
    A22_file = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Ashall2022/original/Ashall2022_all'
    A22_savepath = os.path.dirname(os.path.dirname(A22_file))
    A22_tbl_MWcor = ascii.read(A22_file, format ='fixed_width')

# For A22 table, already MW corrected. Thus reddeing A22 table to make them all same 
A22_tbl = correct_mw_extinction(A22_tbl_MWcor, None, ra, dec, 3.10, dereddening = False)
maglist = []
for mag, filter_ in zip(A22_tbl['mag'], A22_tbl['filter']):
    magsys = 'AB'
    if filter_ in ['B','U','V']:
        magsys = 'Vega'
    mag_instance = ABVegaMagnitude(magnitude = mag, magsys = magsys, filter_ = filter_)
    maglist.append(mag_instance.AB)
A22_tbl['mag'] = maglist
A22_tbl['[mag_cor-mag_origin]'] = 0
A22_tbl.write(filename = A22_savepath+'/Ashall2022_noextin.dat', format = 'ascii.fixed_width', overwrite = True)
#%%
# For H22 table, UBV band is Vega magnitude. Convert before use

H22_tbl = H22_tbl_raw[H22_tbl_raw['observatory'] == 'LasCumbres1m']
maglist = []
for mag, filter_ in zip(H22_tbl['mag'], H22_tbl['filter']):
    magsys = 'AB'
    if filter_ in ['B','U','V']:
        magsys = 'Vega'
    mag_instance = ABVegaMagnitude(magnitude = mag, magsys = magsys, filter_ = filter_)
    maglist.append(mag_instance.AB)
H22_tbl['mag'] = maglist
H22_tbl.write(filename = H22_savepath+'/Hosseinzadeh2022_noextin.dat', format = 'ascii.fixed_width', overwrite = True)
#%%
correct_mw_extinction(IMSNG_tbl, ra = ra, dec = dec, filename = IMSNG_savepath+'/IMSNG_mwextin3.10.dat')
correct_host_extinction(IMSNG_tbl, ebv, Rv = 3.10, filename = IMSNG_savepath+'/IMSNG_noextin.dat', dereddening = 'noextin')
correct_host_extinction(IMSNG_tbl, ebv, Rv = 2.36, filename = IMSNG_savepath+'/IMSNG_hostextin2.36.dat', dereddening = True)
correct_host_extinction(IMSNG_tbl, ebv, Rv = 3.10, filename = IMSNG_savepath+'/IMSNG_hostextin3.10.dat', dereddening = True)
#correct_host_extinction(IMSNG_tbl, ebv, Rv = 2.36, filename = IMSNG_savepath+'/IMSNG_hostiextin2.36_inv.dat', dereddening = False)
#correct_host_extinction(IMSNG_tbl, ebv, Rv = 3.10, filename = IMSNG_savepath+'/IMSNG_hostextin3.10_inv.dat', dereddening = False)

correct_mw_extinction(H22_tbl, ra = ra, dec = dec, filename = H22_savepath+'/Hosseinzadeh2022_mwextin3.10.dat')
correct_host_extinction(H22_tbl, ebv, Rv = 3.10, filename = H22_savepath+'/Hosseinzadeh2022_noextin.dat', dereddening = 'noextin')
correct_host_extinction(H22_tbl, ebv, Rv = 2.36, filename = H22_savepath+'/Hosseinzadeh2022_hostextin2.36.dat', dereddening = True)
correct_host_extinction(H22_tbl, ebv, Rv = 3.10, filename = H22_savepath+'/Hosseinzadeh2022_hostextin3.10.dat', dereddening = True)
#correct_host_extinction(H22_tbl, ebv, Rv = 2.36, filename = H22_savepath+'/Hosseinzadeh2022_hostiextin2.36_inv.dat', dereddening = False)
#correct_host_extinction(H22_tbl, ebv, Rv = 3.10, filename = H22_savepath+'/Hosseinzadeh2022_hostextin3.10_inv.dat', dereddening = False)

correct_mw_extinction(A22_tbl, ra = ra, dec = dec, filename = A22_savepath+'/Ashall2022_mwextin3.10.dat')
correct_host_extinction(A22_tbl, ebv, Rv = 3.10, filename = A22_savepath+'/Ashall2022_noextin.dat', dereddening = 'noextin')
correct_host_extinction(A22_tbl, ebv, Rv = 2.36, filename = A22_savepath+'/Ashall2022_hostextin2.36.dat', dereddening = True)
correct_host_extinction(A22_tbl, ebv, Rv = 3.10, filename = A22_savepath+'/Ashall2022_hostextin3.10.dat', dereddening = True)
#correct_host_extinction(A22_tbl, ebv, Rv = 2.36, filename = A22_savepath+'/Ashall2022_hostiextin2.36_inv.dat', dereddening = False)
#correct_host_extinction(A22_tbl, ebv, Rv = 3.10, filename = A22_savepath+'/Ashall2022_hostextin3.10_inv.dat', dereddening = False)


correct_both_extinction(IMSNG_tbl, hostebv = ebv, hostRv = 2.36, ra = ra, dec = dec, filename = IMSNG_savepath+'/IMSNG_hostmwextin2.36.dat', dereddening = True)
correct_both_extinction(IMSNG_tbl, hostebv = ebv, hostRv = 3.10, ra = ra, dec = dec, filename = IMSNG_savepath+'/IMSNG_hostmwextin3.10.dat', dereddening = True)
#correct_both_extinction(IMSNG_tbl, hostebv = ebv, hostRv = 2.36, ra = ra, dec = dec, filename = IMSNG_savepath+'/IMSNG_hostmwextin2.36_inv.dat', dereddening = False)
#correct_both_extinction(IMSNG_tbl, hostebv = ebv, hostRv = 3.10, ra = ra, dec = dec, filename = IMSNG_savepath+'/IMSNG_hostmwextin3.10_inv.dat', dereddening = False)

correct_both_extinction(H22_tbl, hostebv = ebv, hostRv = 2.36, ra = ra, dec = dec, filename = H22_savepath+'/Hosseinzadeh2022_hostmwextin2.36.dat', dereddening = True)
correct_both_extinction(H22_tbl, hostebv = ebv, hostRv = 3.10, ra = ra, dec = dec, filename = H22_savepath+'/Hosseinzadeh2022_hostmwextin3.10.dat', dereddening = True)
#correct_both_extinction(H22_tbl, hostebv = ebv, hostRv = 2.36, ra = ra, dec = dec, filename = H22_savepath+'/Hosseinzadeh2022_hostmwextin2.36_inv.dat', dereddening = False)
#correct_both_extinction(H22_tbl, hostebv = ebv, hostRv = 3.10, ra = ra, dec = dec, filename = H22_savepath+'/Hosseinzadeh2022_hostmwextin3.10_inv.dat', dereddening = False)

correct_both_extinction(A22_tbl, hostebv = ebv, hostRv = 2.36, ra = ra, dec = dec, filename = A22_savepath+'/Ashall2022_hostmwextin2.36.dat', dereddening = True)
correct_both_extinction(A22_tbl, hostebv = ebv, hostRv = 3.10, ra = ra, dec = dec, filename = A22_savepath+'/Ashall2022_hostmwextin3.10.dat', dereddening = True)
#correct_both_extinction(A22_tbl, hostebv = ebv, hostRv = 2.36, ra = ra, dec = dec, filename = A22_savepath+'/Ashall2022_hostmwextin2.36_inv.dat', dereddening = False)
#correct_both_extinction(A22_tbl, hostebv = ebv, hostRv = 3.10, ra = ra, dec = dec, filename = A22_savepath+'/Ashall2022_hostmwextin3.10_inv.dat', dereddening = False)
# %% 
