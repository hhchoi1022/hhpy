


#%%
import os
import numpy as np
from astropy.table import Table, MaskedColumn
from astropy.io import ascii
import matplotlib.pyplot as plt
from Research.analysis.observedphot import ObservedPhot
from Research.helper import Helper
from observedphot import ObservedPhot
from lmfit import Parameters, minimize
#%% Observation
helper = Helper()
filepath_all = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/all_IMSNG.dat'
filepath_atlas = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/ATLAS.dat'

# Query data
obs_tbl = ascii.read(filepath_all, format = 'fixed_width')
obs_tbl_atlas = ascii.read(filepath_atlas)
mag_max = 16.5
mag_min = 14
obs_tbl_atlas_new = Table()
obs_tbl_atlas_new['obsdate'] = obs_tbl_atlas['MJD']
obs_tbl_atlas_new['mag'] = MaskedColumn([obs_tbl_atlas['m'][i] if (obs_tbl_atlas['m'][i] > mag_min) and (obs_tbl_atlas['m'][i] < mag_max) else None for i in range(len(obs_tbl_atlas))], mask=[not (mag_min < obs_tbl_atlas['m'][i] < mag_max) for i in range(len(obs_tbl_atlas))])
obs_tbl_atlas_new['e_mag'] = MaskedColumn([obs_tbl_atlas['dm'][i] if (obs_tbl_atlas['m'][i] > mag_min) and (obs_tbl_atlas['m'][i] < mag_max) else None for i in range(len(obs_tbl_atlas))], mask=[not (mag_min < obs_tbl_atlas['m'][i] < mag_max) for i in range(len(obs_tbl_atlas))])
obs_tbl_atlas_new['magsys'] = 'AB'
obs_tbl_atlas_new['filter'] = [f'ATLAS_{filt_}' for filt_ in obs_tbl_atlas['F']]
obs_tbl_atlas_new['depth_5sig'] = [depth for depth in obs_tbl_atlas['mag5sig']]
obs_tbl_atlas_new['observatory'] = 'ATLAS'
obs_tbl_atlas_new['detected'] = [True if (obs_tbl_atlas['m'][i] > mag_min) and (obs_tbl_atlas['m'][i] < mag_max) else False for i in range(len(obs_tbl_atlas))]
plt.scatter(obs_tbl_atlas_new['obsdate'], obs_tbl_atlas_new['mag'])
obs_tbl_atlas_new.write('/mnt/data1/supernova_rawdata/SN2023rve/analysis/ATLAS_new.dat', format = 'ascii.fixed_width', overwrite = True)
#%%
# Exclude some observatories with poor photometry
observed_data = ObservedPhot(obs_tbl_atlas_new)
#%%
ax1, ax2 = observed_data.show_lightcurve(phase_binsize = 20,
                              scatter_linewidth=0.5, 
                              scatter_size=50, 
                              errorbar_linewidth=0.5, 
                              errorbar_capsize=0.1, 
                              color_UB = False,
                              color_BV = False, 
                              color_gr = False, 
                              UL = True, 
                              UL_alpha = 0.8,
                              label = True, 
                              label_location=0)
# %%
