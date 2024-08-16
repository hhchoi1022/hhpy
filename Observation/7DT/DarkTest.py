#%%
from astropy.io import fits
import ccdproc
import os
import glob
from ccdproc import CCDData
from ccdproc import Combiner
import numpy as np
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astropy.wcs import FITSFixedWarning

from ccdproc import ImageFileCollection
from HHsupport_phot import combine_img
from astropy.table import Table
#%%

class Prepcoess:
    
    def __init__(self, 
                 imdir,
                 calibdir = None
                ):
        
        self.imdir = imdir
        self.fitsinfo = self.get_info_images(self.imdir)
        self.calibdir = calibdir
        self.calibinfo = Table()
        if not calibdir == None:
            self.calibinfo = self.get_info_images(self.calibdir)
        
    def get_info_images(self, directory):
        '''
        parameters
        ----------
        {the directory absolute path} containing fits images
        
        returns 
        -------
        fits file information in astropy table format
        
        notes 
        -----
        -----
        '''
        iminfo = ImageFileCollection(directory).summary
        absfiles = []
        for file in iminfo['file']:
            absfiles.append(directory + file)
        iminfo['file'] = absfiles
        return iminfo
    
    def master_bias(self, 
                    bias_collection,
                    clip = 'sigma',
                    combine = 'median',
                    
                    # Clipping
                    clip_sigma_low = 2,
                    clip_sigma_high = 5,
                    clip_minmax_min = 3,
                    clip_minmax_max = 3,
                    clip_extrema_nlow = 1,
                    clip_extrema_nhigh = 1
                    ):
        biasfiles = bias_collection['file']
        
        mbias = combine_img(biasfiles,
                            clip = clip,
                            combine = combine,
                            clip_sigma_low = clip_sigma_low,
                            clip_sigma_high = clip_sigma_high,
                            clip_minmax_min = clip_minmax_min,
                            clip_minmax_max = clip_minmax_max,
                            clip_extrema_nlow = clip_extrema_nlow,
                            clip_extrema_nhigh = clip_extrema_nhigh
                            )
        return mbias
    
    def master_dark(self, 
                    dark_collection, 
                    mbias_key = 'com_BIAS',
                    clip = 'sigma',
                    combine = 'median',
                    save_single_calib : bool = True,             
                    # Clipping
                    clip_sigma_low = 2,
                    clip_sigma_high = 5,
                    clip_minmax_min = 3,
                    clip_minmax_max = 3,
                    clip_extrema_nlow = 1,
                    clip_extrema_nhigh = 1
                    ):
        mbias = self.imdir + mbias_key + '*.fits'
        mbias_file = glob.glob(mbias)[0]
        mbias_data = CCDData.read(mbias_file, unit='adu')
        groups_exp = dark_collection.group_by('exptime')
        list_mdark_data = []
        list_mdark_hdr = []
        for group_exp in groups_exp.groups:
            exp_filelist = group_exp['file']
            exptime = group_exp[0]['exptime']
            exp_datalist = []
            for file in exp_filelist:
                dark = CCDData.read(file, unit='adu')
                dark.header['MBIAS'] = mbias
                b_subtracted = ccdproc.subtract_bias(dark, mbias_data)
                sub_filename = os.path.dirname(file) + '/B_' + os.path.basename(file)
                fits.writeto(sub_filename, data = b_subtracted.data.astype(np.int16), header = b_subtracted.header)
                exp_datalist.append(sub_filename)
            # combine_img(exp_datalist, 
            #             clip = clip,
            #             combine = combine,
            #             clip_sigma_low = clip_sigma_low,
            #             clip_sigma_high = clip_sigma_high,
            #             clip_minmax_min = clip_minmax_min,
            #             clip_minmax_max = clip_minmax_max,
            #             clip_extrema_nlow = clip_extrema_nlow,
            #             clip_extrema_nhigh = clip_extrema_nhigh
            #             )

#%%
processed_key_combias = '/com_BIAS*'
processed_key_dark = '/B_DARK*'
path = '/home/hhchoi1022/Desktop/Dark_test/'
dirlist = os.listdir(path)
for dirname in dirlist:
    imdir = path + dirname +'/'
    p = Prepcoess(imdir = imdir)
    bias_key = path + dirname + processed_key_combias
    dark_key = path + dirname + processed_key_dark
    bias_file = glob.glob(bias_key)
    dark_file = glob.glob(dark_key)
    if len(bias_file) > 0:
        pass
    else:
        print(f'='*30+'Processing {dirname}...')
        p = Prepcoess(imdir = imdir)
        bias_collection = p.fitsinfo[p.fitsinfo['imagetyp'] == 'BIAS']
        p.master_bias(bias_collection)
    if len(dark_file) > 0:
        pass
    else:
        print('='*30+'Processing {dirname}...')
        p = Prepcoess(imdir = imdir)
        dark_collection = p.fitsinfo[p.fitsinfo['imagetyp'] == 'DARK']
        p.master_dark(dark_collection)
        
#%%


#%%
bias_collection = p.fitsinfo[p.fitsinfo['imagetyp'] == 'BIAS']
p.master_bias(bias_collection)
#%%
dark_collection = p.fitsinfo[p.fitsinfo['imagetyp'] == 'DARK']
p.master_dark(dark_collection)
# %%

file_key = '/home/hhchoi1022/Desktop/Dark_test/2023-11-20/B_DARK*.00s*.fits'
dark_filelist = glob.glob(file_key)
data_array = []
hdr_array = []
import numpy as np

# Create a sample 2D array

# Define the percentiles
lower_percentile = 99.99945
upper_percentile = 99.9995
data = fits.getdata(dark_filelist[0])
# Find the values corresponding to the percentiles
lower_threshold = np.percentile(data, lower_percentile)
upper_threshold = np.percentile(data, upper_percentile)

# Create a boolean mask for values within the desired range
mask = (data >= lower_threshold) & (data <= upper_threshold)

# Get the indices of the elements within the desired range
indices = np.column_stack(np.where(mask))
indices = np.array([[6108,7968],[1,1]])
selected_values = data[indices[:, 0], indices[:, 1]]
plt.plot(selected_values)
#%%
for file in dark_filelist:
    print(file)
    data = fits.getdata(file)
    data_array.append(data[indices[:, 0], indices[:, 1]])
    
    hdr = fits.getheader(file)
    hdr_array.append(hdr)
# %%
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import datetime
from datetime import timedelta
import numpy as np
x = 6408
y = 3371
from astropy.table import Table
tbl = Table()

vallist = []
timelist = []
dark_is_morning = []
exptimelist = []
ccdtemplist = []
for data, hdr in zip(data_array, hdr_array):
    vallist = []
    for i in range(len(data)):
        val = data[i]
        vallist.append(val)
    timelist.append(Time(hdr['DATE-OBS']).to_datetime())
    dark_hour = Time(hdr['DATE-LOC']).to_datetime().hour
    exptime = hdr['EXPTIME']
    ccdtemp = hdr['CCD-TEMP']
    if (ccdtemp > -21) & (ccdtemp < -19):
        ccdtemp = 20
    if (ccdtemp > -16) & (ccdtemp < -14):
        ccdtemp = -15
    else:
        ccdtemp = -10
    if (exptime > 8) & (exptime < 12):
        exptime = 10
    if (exptime > 13) & (exptime < 17):
        exptime = 15
    if (exptime > 28) & (exptime < 32):
        exptime = 30
    if (exptime > 58) & (exptime < 62):
        exptime= 60
        
    exptimelist.append(exptime)
    ccdtemplist.append(ccdtemp)
    if dark_hour > 12:
        is_morningdark = 'b'
    else:
        is_morningdark = 'r'
    dark_is_morning.append(is_morningdark)
tbl['exptime'] = exptimelist
tbl['ccdtemp'] = ccdtemplist
tbl['obsdate'] = timelist
tbl['dark_is_morning'] = dark_is_morning

for i in range(len(data)):
    val = np.array(data_array)[:,i]
    indice_name = str(indices[i])
    tbl[indice_name] = val
#%%
tbl.sort('obsdate')
tbl_60s = tbl[tbl['exptime'] == 60]
tbl_30s = tbl[tbl['exptime'] == 30]
tbl_15s = tbl[tbl['exptime'] == 15]
tbl_10s = tbl[tbl['exptime'] == 10]

#%%
tbl_15deg =tbl[tbl['ccdtemp'] == -15]
tbl_10deg =tbl[tbl['ccdtemp'] == -20]
tbl_20deg =tbl[tbl['ccdtemp'] == -10]
#%%
#%%
plt.figure(dpi = 150, figsize = (6,4))
#plt.title(f'CCDTEMP = {round(hdr["CCD-TEMP"],2)} OBSTIME = {timelist[0]}')
plt.title(f'Pixel valie [{str(indices[-2])}]')
plt.scatter(tbl_60s['obsdate'], tbl_60s[str(indices[-2])], facecolor = 'none', edgecolors = tbl_60s['dark_is_morning'], marker = '^')
#plt.scatter(tbl_60s['obsdate'], tbl_60s[str(indices[1])], facecolor = 'none', edgecolors = tbl_60s['dark_is_morning'], marker = 'o')
#plt.scatter(tbl_60s['obsdate'], tbl_60s[str(indices[2])], facecolor = 'none', edgecolors = tbl_60s['dark_is_morning'], marker = 's')

plt.xticks(rotation = 45)    
#plt.axhline(y = np.mean(vallist))
#plt.xlim(np.max(timelist)-timedelta(minutes = 670),np.max(timelist)-timedelta(minutes= 660))
#plt.xlim(np.max(timelist)-timedelta(hours = 3),np.min(timelist)+timedelta(hours= 3))
#plt.ylim(3000,3500)
#plt.xlim(np.min(timelist)-timedelta(days = 3),np.max(timelist)+timedelta(days = 3))
plt.ylabel('Value')
plt.xlabel('obsdate')
#%%
plt.figure(dpi = 150, figsize = (6,4))
#plt.title(f'CCDTEMP = {round(hdr["CCD-TEMP"],2)} OBSTIME = {timelist[0]}')
plt.title(f'Pixel valie [{str(indices[-2])}]')
plt.scatter(tbl['exptime'], tbl[str(indices[-2])], facecolor = 'none', edgecolors = tbl_60s['dark_is_morning'], marker = '^')
#plt.scatter(tbl_60s['obsdate'], tbl_60s[str(indices[1])], facecolor = 'none', edgecolors = tbl_60s['dark_is_morning'], marker = 'o')
#plt.scatter(tbl_60s['obsdate'], tbl_60s[str(indices[2])], facecolor = 'none', edgecolors = tbl_60s['dark_is_morning'], marker = 's')

plt.xticks(rotation = 45)    
#plt.axhline(y = np.mean(vallist))
#plt.xlim(np.max(timelist)-timedelta(minutes = 670),np.max(timelist)-timedelta(minutes= 660))
#plt.xlim(np.max(timelist)-timedelta(hours = 3),np.min(timelist)+timedelta(hours= 3))
#plt.ylim(3000,3500)
#plt.xlim(np.min(timelist)-timedelta(days = 3),np.max(timelist)+timedelta(days = 3))
plt.ylabel('Value')
plt.xlabel('exptime')
#%%
plt.title(f'CCDTEMP = {round(hdr["CCD-TEMP"],2)} OBSTIME = {timelist[0]}')
plt.scatter(exptimelist, vallist, c = dark_is_morning)
#plt.scatter(timelist, vallist, c = 'b')
plt.xticks(rotation = 45)  
plt.ylabel(f'Pixel valie [{x},{y}]')
plt.xlabel('exptime')
#plt.axhline(y = np.mean(vallist))
#plt.xlim(np.min(timelist)-timedelta(days = 3),np.max(timelist)+timedelta(days = 3))
#plt.xlim(np.min(timelist)-timedelta(minutes = 3),np.max(timelist)+timedelta(minutes= 3))

#%%
#plt.title(f'CCDTEMP = {round(hdr["CCD-TEMP"],2)} OBSTIME = {timelist[0]}')
plt.scatter(ccdtemplist, vallist, c = dark_is_morning)
#plt.scatter(timelist, vallist, c = 'b')
plt.xticks(rotation = 45)  
plt.title(f'Pixel valie [{x},{y}]')
plt.xlabel('ccdtemp')
#plt.axhline(y = np.mean(vallist))
#plt.xlim(np.min(timelist)-timedelta(days = 3),np.max(timelist)+timedelta(days = 3))
#plt.xlim(np.min(timelist)-timedelta(minutes = 3),np.max(timelist)+timedelta(minutes= 3))




# %%



import os
path = '/home/h'