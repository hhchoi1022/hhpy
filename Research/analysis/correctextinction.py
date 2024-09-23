#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 09:25:31 2022

@author: hhchoi1022
"""
#%%
import numpy as np
import extinction
import glob, os
from astropy.table import Table
from astropy.table import vstack
import re
from astropy.io import ascii

from Research.helper import Helper
#%%
class CorrectExtinction:
    
    def __init__(self,
                 filename : str):
        self.helper = Helper()
        self.filename = filename
        self.original_data = self._load_table()
        self.corrected_data = self.original_data.copy()
        self.correction_history = list()
    
    def _load_table(self):
        tbl = ascii.read(self.filename, format = 'fixed_width')
        return tbl  

    def correct_host_extinction(self,  
                                ebv : float,
                                Rv : float = 3.1,
                                dereddening = True): 
        corrected_data = self.corrected_data.copy()
        import pyphot
        lib = pyphot.get_library() # Filter library 
        _, _, _, filter_key, _ = self.helper.load_filt_keys()
        filterlist = sorted(list(set(corrected_data['filter'])))
        Av = Rv * ebv
        corrected_data['A_filter'] = -99.9
        for filtname in filterlist:
            if filtname not in filter_key.keys():
                print(f'filter "{filtname}" is not in filter_key')
                pass
            else:
                filter_ = lib[filter_key[filtname]]
                lpivot_AA = np.array([filter_.lpivot.value])
                Afilt = round(extinction.fitzpatrick99(lpivot_AA, Av, Rv)[0],3)
                
                if dereddening == True:
                    corrected_data['mag'][corrected_data['filter'] == filtname] -= Afilt
                    corrected_data['A_filter'][corrected_data['filter'] == filtname] = -Afilt
                    is_corrected = 'Host_dereddening'
                elif dereddening == False:
                    corrected_data['mag'][corrected_data['filter'] == filtname] += Afilt
                    corrected_data['A_filter'][corrected_data['filter'] == filtname] = Afilt
                    is_corrected = 'Host_reddening'
                else:
                    raise ValueError(f'dereddening = {dereddening}')
        self.correction_history.append(is_corrected)
        corrected_data['mag'] = corrected_data['mag'].round(3)
        corrected_data['A_filter'] = corrected_data['A_filter'].round(3)
        self.corrected_data = corrected_data

    def correct_mw_extinction(self, 
                              ra : float = 64.9708333, 
                              dec : float = -54.9480556, 
                              mwRv : float = 3.1, 
                              dereddening = True,
                              SFDmap_path : str = "./sfddata-master"): 
        #ebv = 0.097 # Hosseinzadeh 2022
        import sfdmap
        import pyphot
        
        if not os.path.exists(SFDmap_path):
            raise ValueError(f'cannot find SFD map at {SFDmap_path}')
        dustmap = sfdmap.SFDMap(SFDmap_path)
        if ra is not None:
            try:
                mwebv = dustmap.ebv(ra, dec)
            except:
                raise ValueError('cannot calculate extinction with given coordinates')
        
        corrected_data = self.corrected_data.copy()
        lib = pyphot.get_library() # Filter library 
        _, _, _, filter_key, _ = self.helper.load_filt_keys()
        filterlist = sorted(list(set(corrected_data['filter'])))
        Av = mwRv * mwebv
        corrected_data['A_filter'] = -99.9
        for filtname in filterlist:
            if filtname not in filter_key.keys():
                print(f'filter "{filtname}" is not in filter_key')
                pass
            else:
                filter_ = lib[filter_key[filtname]]
                lpivot_AA = np.array([filter_.lpivot.value])
                Afilt = round(extinction.fitzpatrick99(lpivot_AA, Av, mwRv)[0],3)
                if dereddening == True:
                    corrected_data['mag'][corrected_data['filter'] == filtname] -= Afilt
                    corrected_data['A_filter'][corrected_data['filter'] == filtname] = -Afilt
                    is_corrected = 'MW_dereddening'
                    
                elif dereddening == False:
                    corrected_data['mag'][corrected_data['filter'] == filtname] += Afilt
                    corrected_data['A_filter'][corrected_data['filter'] == filtname] = Afilt
                    is_corrected = 'MW_reddening'
                else:
                    raise ValueError(f'dereddening = {dereddening}')
        self.correction_history.append(is_corrected)
        corrected_data['mag'] = corrected_data['mag'].round(3)
        corrected_data['A_filter'] = corrected_data['A_filter'].round(3)
        self.corrected_data = corrected_data
        
    
    def save(self):
        filename = self.filename.split('.')[0]+'_'+'_'.join(self.correction_history) + '.dat'
        self.corrected_data.write(filename = filename, format = 'ascii.fixed_width', overwrite = True)
        print(f'saved as {filename}')   

#%%
if __name__ == '__main__':
    ebv = 0.1#097  # Hosseinzadeh 2022
    ra = 64.9708333
    dec = -54.9480556
#%%
# data from Ahsall 2022 is already corrected for MW extinction >>> noextin_dat
if __name__ == '__main__':
    IMSNG_file = '/data1/supernova_rawdata/SN2021aefx/photometry/all_IMSNG.dat'
    C = CorrectExtinction(IMSNG_file)
    C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = True)
    #C.correct_host_extinction(ebv = ebv, Rv = 2.3)
    C.save()
#%% # For A22 table, already MW corrected. Thus reddeing A22 table to make them all same 
if __name__ == '__main__':
    #A22_file = '/data1/supernova_rawdata/SN2021aefx/photometry/Ashall2022.dat'
    #C = CorrectExtinction(A22_file)
    #C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = False)
    #C.save()
    
    A22_file = '/data1/supernova_rawdata/SN2021aefx/photometry/Ashall2022.dat'
    C = CorrectExtinction(A22_file)
    C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = True)
    C.correct_host_extinction(ebv = ebv, Rv = 2.3)
    C.save()
#%% 
if __name__ == '__main__':
    H22_file = '/data1/supernova_rawdata/SN2021aefx/photometry/Hosseinzadeh2022.dat'
    C = CorrectExtinction(H22_file)
    C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = True)
    C.correct_host_extinction(ebv = ebv, Rv = 2.3, dereddening= True)
    C.save()
#%%
if __name__ == '__main__':
    spec_file = '/data1/supernova_rawdata/SN2021aefx/photometry/all_spec.dat'
    C = CorrectExtinction(spec_file)
    #C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = True)
    C.correct_host_extinction(ebv = ebv, Rv = 2.3, dereddening= True)
    #C.save()
#%%
if __name__ == '__main__':
    spec_file = '/data1/supernova_rawdata/SN2021aefx/photometry/all_phot.dat'
    C = CorrectExtinction(spec_file)
    C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = True)
    C.correct_host_extinction(ebv = ebv, Rv = 2.3, dereddening= True)
    C.save()
#%%
if __name__ == '__main__':
    spec_file = '/data1/supernova_rawdata/SN2021aefx/photometry/all_photspec.dat'
    C = CorrectExtinction(spec_file)
    C.correct_mw_extinction(ra = ra, dec = dec, mwRv = 3.10, dereddening = True)
    C.correct_host_extinction(ebv = ebv, Rv = 2.3, dereddening= True)
    C.save()
# %%
