#%%
from astropy.io import ascii
from observedphot import ObservedPhot
from astropy.table import vstack
import glob
from astropy.table import Table
import re
#%% IMSNG data
obslist = glob.glob('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/IMSNG/original/*.dat')
all_tbl = Table()
for obs in obslist:
    tbl = Table.read(obs, format = 'ascii.fixed_width')
    all_tbl = vstack([all_tbl, tbl])
all_tbl['status'] = 'detected'
try:
    all_tbl['status'][all_tbl['mag'] == all_tbl['UL5_4']] = 'UL'
except:
    pass
obs_tbl = Table()
obs_tbl['obsdate'] = all_tbl['obsdate'].round(4)
obs_tbl['filter'] = all_tbl['filter']
obs_tbl['mag'] = all_tbl['mag'].round(3)
obs_tbl['e_mag'] = all_tbl['e_mag'].round(3)
obs_tbl['observatory'] = all_tbl['observatory']
obs_tbl['status'] = all_tbl['status']
obs_tbl.write('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/IMSNG/original/IMSNG_all', format = 'ascii.fixed_width', overwrite = True)
#%% Hosseinzadeh 2022
"""
UBV(Vega), gri(AB). Deal with caution!!
"""
hossefile = glob.glob('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/original/*.tex')[0]
hossedata = ascii.read(hossefile)
hossedata['Filter'][hossedata['Filter'] == '\\nodata'] = 'Unfiltered'
maglist = []
errlist = []
tellist = []
for i, magerr in enumerate(hossedata['Apparent Magnitude']):
    try:
        mag = float(re.findall(pattern = r'(\d{2}.\d{1,2}) \\pm (\d.\d{1,2})', string = magerr)[0][0])
        err = float(re.findall(pattern = r'(\d{2}.\d{1,2}) \\pm (\d.\d{1,2})', string = magerr)[0][1])
    except:
        mag = float(re.findall(pattern = r'(\d{2}.\d{1,2})', string = magerr)[0])
        err = 100
    maglist.append(float(mag))
    errlist.append(err)
    tellist.append(hossedata[i]['Telescope'].replace(" ",""))
hosse_tbl = Table()
hosse_tbl['obsdate'] = hossedata['MJD']
hosse_tbl['filter'] = hossedata['Filter']
hosse_tbl['mag'] = maglist
hosse_tbl['e_mag'] = errlist
hosse_tbl['observatory'] = tellist
hosse_tbl['status'] = 'detected'
hosse_tbl.write('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/original/Hosseinzadeh2022_all', format = 'ascii.fixed_width', overwrite = True)
#%% Ashall 2022
"""
MW extinction is already corrected for this data (Ashall+ 2022)
"""
ashallfile = glob.glob('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve//Ashall2022/original/Ashall2022_raw.dat')[0]
ashalldata = ascii.read(ashallfile, format = 'fixed_width')
filterset = list(ashalldata.columns)[2:]
obsdates = ashalldata['MJD']
ashall_tbl = Table()
for filt_ in filterset:
    maglist = []
    errlist = []
    filtlist = []
    tellist = []
    for magerr in ashalldata[filt_]:
        try:
            mag = float(re.findall(pattern = '(\d{1,2}.\d{1,2})', string = magerr)[0])
            err = float(re.findall(pattern = '(\d{1,2}.\d{1,2})', string = magerr)[1])
        except:
            mag = 100
            err = 100
        filtlist.append(filt_)
        maglist.append(float(mag))
        errlist.append(err)
        tellist.append('Swope')
    ashall_tbl_tmp = Table()
    ashall_tbl_tmp['obsdate'] = obsdates
    ashall_tbl_tmp['filter'] = filtlist
    ashall_tbl_tmp['mag'] = maglist
    ashall_tbl_tmp['e_mag'] = errlist
    ashall_tbl_tmp['observatory'] = tellist
    ashall_tbl_tmp['status'] = 'detected'
    ashall_tbl = vstack([ashall_tbl, ashall_tbl_tmp])
    ashall_tbl = ashall_tbl[ashall_tbl['mag'] < 30]
ashall_tbl.write('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Ashall2022/original/Ashall2022_all', format = 'ascii.fixed_width', overwrite = True)
#%%
obs_tbl = vstack([obs_tbl, hosse_tbl, ashall_tbl])
#%%
#tbl_obs = vstack([tbl_ashall, tbl_IMSNG])
obs = ObservedPhot(data_tbl = obs_tbl, target = 'SN2021aefx', observer = 'IMSNG', MW_extinction_corrected= False, Host_extinction_corrected= False)
#obs.exclude_filter('U')
#obs.exclude_filter('B')
#obs.exclude_filter('V')
#obs.exclude_filter('r')
#obs.exclude_filter('u')
#obs.exclude_filter('i')
obs.exclude_observatory('LasCumbres0.4m')
import matplotlib.pyplot as plt
plt.figure(dpi = 500)
#ashall.show_lightcurve()
obs.show_lightcurve(label = True, color_BV=False, color_gr = False, scatter_size= 10, errorbar_capsize=1)
plt.xlim(59625, 59635)
#%%
