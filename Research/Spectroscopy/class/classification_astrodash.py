#%%
import os
import numpy as np
WD = os.getcwd()
import sncosmo
#%%
import astropy.units as u
from specutils import Spectrum1D
from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)
from astropy.coordinates import SkyCoord
from HHsupport_phot import to_skycoord
import sncosmo
from HHsupport_analysis import formatting_SNcosmo
from tns_search import tns_search
from astropy.time import Time
import numpy as np
from astropy.io import ascii
#%%#----------------------------------------------------------------------------------
tns = tns_search(photdata = True, specdata = True)
data = tns.get_obj('2022xkq')
#%%

ra = data['data']['reply']['ra']
dec = data['data']['reply']['dec']
coord = to_skycoord(ra, dec)
#%%
mjdlist = []
maglist = []
magerrlist = []
magsyslist = []
filterlist = []
observerlist = []

for photometry in data['data']['reply']['photometry']:
    if isinstance(photometry['flux'], float):
        mjd = Time(photometry['jd'],format = 'jd').mjd
        mag = photometry['flux']
        magerr = photometry['fluxerr']
        if (isinstance(magerr, str)) or (magerr == 0):
            magerr = 0.05
        magsys = photometry['flux_unit']['name']
        if 'ab' in magsys.lower():
            magsys = 'ab'
        if 'vega' in magsys.lower():
            magsys = 'vega'
        filter_ = photometry['filters']['name']
        observer = photometry['instrument']['name']
        mjdlist.append(mjd)
        maglist.append(mag)
        magerrlist.append(magerr)
        magsyslist.append(magsys)
        filterlist.append(filter_)
        observerlist.append(observer)
#%% SNcosmo
filterlist = ['gaia::g', 'gaia::g', 'ps1::w', 'atlaso', 'gaia::g']
# KCT
obstime_KCT = [Time(2459922.808570138, format = 'jd').mjd, Time(2459922.824602396, format= 'jd').mjd]
filter_KCT = ['g', 'r']
mag_kct = [17.2842,16.5076]
magerr_kct = [0.1373, 0.0891]
magsyslist_kct = ['ab','ab']
observer_kct = ['KCT','KCT']
filter_kct = ['sdssg', 'sdssr']
mjdlist = mjdlist + obstime_KCT 
maglist = maglist +mag_kct
magerrlist = magerrlist +magerr_kct
filterlist = filterlist + filter_kct
magsyslist = magsyslist + magsyslist_kct
observerlist = observerlist + observer_kct

# LSGT
obstime_LSGT = [Time(2459922.100798611, format = 'jd').mjd, Time(2459922.1233449075, format= 'jd').mjd, Time(2459922.112060185, format = 'jd').mjd]
mag_LSGT = [17.6162,16.6097, 16.7988]
magerr_LSGT = [0.1201, 0.139, 0.09]
magsyslist_LSGT = ['ab','ab', 'ab']
observer_LSGT = ['LSGT','LSGT', 'LSGT']
filter_LSGT = ['sdssg', 'sdssi', 'sdssr']
mjdlist = mjdlist + obstime_LSGT 
maglist = maglist +mag_LSGT
magerrlist = magerrlist +magerr_LSGT
filterlist = filterlist + filter_LSGT
magsyslist = magsyslist + magsyslist_LSGT
'''
# SAO_spec
obstime_LSGT = [59885.71478, 59885.71478, 59885.71478]
mag_LSGT = [15.702,15.614, 14.807]
magerr_LSGT = [0.1, 0.1, 0.1]
magsyslist_LSGT = ['ab','ab', 'ab']
observer_LSGT = ['SAO','SAO', 'SAO']
filter_LSGT = ['sdssg', 'sdssi', 'sdssr']
mjdlist = mjdlist + obstime_LSGT 
maglist = maglist +mag_LSGT
magerrlist = magerrlist +magerr_LSGT
filterlist = filterlist + filter_LSGT
magsyslist = magsyslist + magsyslist_LSGT

# ZTF_spec
obstime_LSGT = [59869.31593749998, 59869.31593749998, 59869.31593749998]
mag_LSGT = [16.201, 16.100, 16.047]
magerr_LSGT = [0.1, 0.1, 0.1]
magsyslist_LSGT = ['ab','ab', 'ab']
observer_LSGT = ['ZTF','ZTF', 'ZTF']
filter_LSGT = ['sdssg', 'sdssi', 'sdssr']
mjdlist = mjdlist + obstime_LSGT 
maglist = maglist +mag_LSGT
magerrlist = magerrlist +magerr_LSGT
filterlist = filterlist + filter_LSGT
magsyslist = magsyslist + magsyslist_LSGT
'''

fit_tbl = formatting_SNcosmo(mjd = np.array(mjdlist), mag = np.array(maglist), e_mag = np.array(magerrlist), filter_ = np.array(filterlist), magsys = np.array(magsyslist))
fit_tbl.remove_rows([0,2, 4])
#%%
b_mag = [14.7, 14.5, 14.4, 14.3, 14.0, 14.6, 14.6, 15.0, 14.0, 15.0, 15.4, 15.7, 16.1]
b_filter_ = ['Clear','Clear','Green','Clear','R','V','V','Clear','Clear','V','Clear','G','Clear']
b_date = [59879.753,59880.000,59882.135,59883.083,59884.125,59884.996, 59885.849, 59886.104, 59886.583, 59889.923, 59900.371, 59908.969, 59909.083]
b_tbl = Table()
b_tbl['date'] = b_date
b_tbl['filter'] = b_filter_
b_tbl['mag'] = b_mag
#%% fitting 
ra = coord.ra.value
dec = coord.dec.value
#%%
source = sncosmo.get_source('salt2')
model = sncosmo.Model(source=source)
dust = sncosmo.CCM89Dust()
import sfdmap
dustmap = sfdmap.SFDMap("../../../config/sfddata-master")
ebv = dustmap.ebv(ra, dec)
model.add_effect(dust, 'mw', 'obs')
model.set(mwebv = ebv)
model.set(mwr_v = 3.1)
model.add_effect(dust, 'host', 'rest')
model.set(hostebv = 0)
model.set(hostr_v = 2.7)
model.set(z = 0.007735)
result , fitted_model= sncosmo.fit_lc(
    fit_tbl, model,
    #['t0', 'amplitude', 'hostebv'], # hsiao
    ['t0', 'x0', 'x1', 'c'], #salt2
    bounds = {}
    )
sncosmo.plot_lc(fit_tbl, model=fitted_model, errors=result.errors,  ncol = 3,  xfigsize = 10, tighten_ylim=True, color = 'black')
import numpy as np
t0 = result.parameters[1]
#%%
fit_tbl['mag'] = -2.5*np.log10(fit_tbl['flux']) + 25
fit_tbl['magerr'] = [0.05,  0.05, 0.1373, 0.0891, 0.1201, 0.139, 0.09]
bandset = list(set(fit_tbl['band']))
bandset = [ 'Orange(ATLAS)', 'White(PS1)', 'g(SDSS)', 'i(SDSS)', 'r(SDSS)']
bandgroup = fit_tbl.group_by('band').groups
t_range  = np.linspace(t0-20, t0+60, 300)
for i, group in enumerate(bandgroup):
    band = group[0]['band']
    maglist = fitted_model.bandmag(band, 'ab', t_range)
    plt.scatter(group['mjd'],group['mag'], label = bandset[i], edgecolors = 'k')
    #plt.errorbar(group['mjd'],group['mag'], yerr = group['magerr'], fmt= 'none')
    plt.plot(t_range, maglist, linewidth  =0.5)
plt.scatter(b_date, b_mag, edgecolor= 'k', facecolors = 'none', label = 'Unknown')
plt.xlabel('MJD')
plt.ylabel('Apparent magnitude')
plt.ylim(20, 13.5)
plt.legend()
#%%

b_filtergroups = b_tbl.group_by('filter')
plt.figure(dpi = 500)
plt.ylim(19, 13.5)
plt.xlabel('MJD')
plt.ylabel('Apparent magnitude')
plt.scatter(fit_tbl['mjd'], -2.5*np.log10(fit_tbl['flux'])+25, c = 'r')
plt.scatter(b_date, b_mag, c = 'k')
plt.axvline(59879.52, linestyle = '--', c='r')
#%%
fit_filtergroups = fit_tbl.group_by('band')
#%%
fit_filtergroups.groups
fit_filterset = set(fit_tbl['band'])
b_filterset = set(b_tbl['filter'])
#%%

#%% Spectrum (download from TNS)
spec_data = []
for spec in data['data']['reply']['spectra']:
    spec_tbl = ascii.read(spec['asciifile'])
    mjd = Time(spec['jd'], format = 'jd').mjd
    phase = mjd-t0
    group = spec['source_group']['name']
    filepath = WD+'/'+os.path.basename(spec['asciifile'])
    spec_data.append([mjd, group, spec_tbl, phase, filepath])
    spec_tbl.write(filepath, format='ascii', overwrite = True)
    
#%%
mjd1, group1, spec_tbl1, phase1, filepath1 = spec_data[0]
mjd2, group2, spec_tbl2, phase2, filepath2 = spec_data[1]
mjd3, group3, spec_tbl3, phase3, filepath3 = spec_data[2]
mjd_SAO = Time(2459886.2147869794, format= 'jd').mjd
group_SAO = 'SAO'
filepath_SAO = WD+'/'+'pSN2022xkq_w_spec.csv'
spec_tbl_SAO = ascii.read(filepath_SAO)
phase_SAO = mjd_SAO-t0
#spec_tbl_SAO['lambsq'] = spec_tbl_SAO['col1']**2*spec_tbl_SAO['col2']
#%%
def spectrum(spec_tbl):
    spec = Spectrum1D(flux = spec_tbl['col2']*u.Unit('erg cm-2 s-1 AA-1') , spectral_axis = spec_tbl['col1']*u.AA)
    return spec
#%%
import matplotlib.pyplot as plt
plt.figure(dpi = 500, figsize = (8,6))
plt.plot(spectrum(spec_tbl1).spectral_axis, spectrum(spec_tbl1).flux/np.max(spectrum(spec_tbl1).flux)+0.25, label = group1, c= 'k')
plt.plot(spectrum(spec_tbl2).spectral_axis, box_smooth(spectrum(spec_tbl2),5).flux/np.max(box_smooth(spectrum(spec_tbl2),5).flux)+0.8, label = group2, c = 'k')
plt.plot(spectrum(spec_tbl3).spectral_axis, box_smooth(spectrum(spec_tbl3),1).flux/np.max(box_smooth(spectrum(spec_tbl3),1).flux)+1.6, c='k')
plt.plot(spectrum(spec_tbl_SAO).spectral_axis, box_smooth(spectrum(spec_tbl_SAO),5).flux/np.max(box_smooth(spectrum(spec_tbl_SAO),5).flux)-0.34, label = group_SAO, c ='r')
plt.plot(spectrum(spec_tbl_SAO).spectral_axis, box_smooth(spectrum(spec_tbl_SAO),1).flux/np.max(box_smooth(spectrum(spec_tbl_SAO),1).flux)-0.27, label = group_SAO, c ='k', alpha = 0.3)
#plt.legend()
#plt.grid()

plt.xlim(4000,8000)
plt.ylabel('Flux (Arbitrary unit)')
plt.xlabel('Wavelength (AA)')
#plt.xlim(5200, 6800)
#%%
import mosfit
# Create an instance of the `Fetcher` class.
my_fetcher = mosfit.fetcher.Fetcher()
SNname = 'SN2000cx'
fetched = my_fetcher.fetch(SNname)[0]
#%%
import json
file = fetched['path']
with open(file, "r") as st_json:
    st_python = json.load(st_json)
print(len(st_python[SNname]['spectra']))
#%%
#for wave, flux, fluxerr in 
for i in range(25):
    wave = np.array(st_python[SNname]['spectra'][i]['data'])[:,0].astype(float)
    flux = np.array(st_python[SNname]['spectra'][i]['data'])[:,1].astype(float)
    
    plt.figure()
    ax = plt.subplot()
    ax.plot(wave, flux/np.max(flux))
    ax1 = ax.twinx()
    ax1.plot(spectrum(spec_tbl_SAO).spectral_axis, box_smooth(spectrum(spec_tbl_SAO),10).flux/np.max(box_smooth(spectrum(spec_tbl_SAO),10).flux)-0.34, label = group_SAO, c ='r')
    plt.xlim(4000,8000)
    plt.plot()
#%%
wave
#%%
plt.plot(spectrum(spec_tbl_SAO).spectral_axis, box_smooth(spectrum(spec_tbl_SAO),20).flux/np.max(box_smooth(spectrum(spec_tbl_SAO),20).flux)-0.34, label = group_SAO, c ='r')

#%% Synthetic photometry 
import pyphot
from HHsupport_analysis import load_filt_keys
from astropy.table import Table
from pyphot import unit
def synth_phot(wl, f_lamb, filter_key ='UBVRIugriz'):
    lib = pyphot.get_library()
    _, _, _, pyphot_key, _ = load_filt_keys(filter_key)
    result_phot = Table()
    result_phot['filter'] = list(filter_key)
    result_phot.add_index('filter')
    magset = []
    for filt_ in filter_key:
        filt_pyphot = lib[pyphot_key[filt_]]
        flux = filt_pyphot.get_flux(wl*unit['AA'],f_lamb*unit['ergs/s/cm**2/AA'], axis = 1)
        mag = -2.5*np.log10(flux.value) - filt_pyphot.AB_zero_mag
        magset.append(mag)
    result_phot['mag'] = magset
    return result_phot
synth_phot(spec_tbl_SAO['col1'], spec_tbl_SAO['col2'])
#%%
import astrodash
import os
#%%
#%%
example = [
    (filepath2, 0.007735)]
#%%
# Create filenames and knownRedshifts lists
filenames = [i[0] for i in example]
knownRedshifts = [i[1] for i in example]
#%%
# Classify all spectra
classification = astrodash.Classify(filenames, knownRedshifts, classifyHost=True, knownZ=True, smooth=15, minWave = 4000, maxWave = 7500 )
#%%
bestFits, redshifts, bestTypes, rlapFlag, matchesFlag, redshiftErrs = classification.list_best_matches(n=10, saveFilename='example_best_fits.txt')
#%%
# Plot sn2013fs from open supernova catalog (2nd spectrum)
classification.plot_with_gui(indexToPlot=0)
# %%
from astropy.io import ascii
# %%
tbl0 = ascii.read('./ZTF_221013.dat')
tbl1 = ascii.read('./ZTF_221017.dat')
tbl2 = ascii.read('./pSN2022xkq_w_spec.csv')

# %%
import numpy as np
tbl0['relative_flux'] = tbl0['col2']/np.mean(tbl0['col2'])#[(tbl0['col1'] > 6400) & (tbl0['col1'] < 6600)])
tbl1['relative_flux'] = tbl1['col2']/np.max(tbl1['col2'])#[(tbl1['col1'] > 6400) & (tbl1['col1'] < 6600)])
tbl2['relative_flux'] = tbl2['col2']/np.max(tbl2['col2'])#[(tbl2['col1'] > 6400) & (tbl2['col1'] < 6600)])
tbl0['lambsqflux'] = tbl0['col2']*tbl0['col1']**2
tbl1['lambsqflux'] = tbl1['col2']*tbl1['col1']**2
tbl2['lambsqflux'] = tbl2['col2']*tbl2['col1']**2
#%%
import matplotlib.pyplot as plt

# %%
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

# %%
# %%
ax0 = plt.subplot()
ax0.plot(tbl0['col1'], smooth(tbl0['col2'],4), c ='r')
ax0.set_xlim(4000, 8000)
#ax0.set_ylim(0,3e-16)
ax2 = ax0.twinx()
ax0.plot(tbl1['col1'], smooth(tbl1['col2'],2), c ='k')
ax2.set_xlim(4000, 8000)
#ax2.set_ylim(0, 3e-15)
ax3 = ax0.twinx()
ax0.plot(tbl2['col1'], smooth(tbl2['col2'],15), c ='b')
ax3.set_xlim(4000, 8000)
ax0.grid(axis = 'x')
#ax3.set_ylim(0, 3e-14)
