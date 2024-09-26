
#%%
from lmfit import Parameters, minimize, report_fit
from astropy.io import ascii
from observedphot import ObservedPhot
import matplotlib.pyplot as plt
import numpy as np
from HHsupport_analysis import mag_to_flux
from HHsupport_analysis import load_filt_keys
from HHsupport_analysis import flux_to_mag
#%% Observation
DM = 31.14
filepath_all = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Alldata_H_M_cor310.dat'
tbl_obs = ascii.read(filepath_all, format = 'fixed_width')
observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
observed_data.exclude_observatory(['Swope', 'KMTNet_Ni2023'])
#%% Fireball model
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
def fireball_filter(params, time, filter_):
    exptime = params['exptime']
    #exptime = params[f'exptime_{filter_}']
    amp = params[f'amp_{filter_}']
    alpha= params[f'alpha_{filter_}']
    return fireball_model(np.array(time), amp, exptime, alpha)
def calc_chisq(params, x, y, e_y, filter_key):
    tot_chisq = []
    for i, filter_ in enumerate(filter_key):
        obstime = x[i]
        obsflux = y[i]
        obserr = e_y[i]
        modelflux = fireball_filter(params, obstime, filter_)
        chisq = (((obsflux - modelflux)/obserr)**2)
        tot_chisq.append(chisq)
    return np.concatenate(tot_chisq)
# %% Fitting
fit_filterset = 'UBgVri'
phase_min = 59529
phase_max = 59538
tbl_filt_all = observed_data.get_filt_data(observed_data.get_data_detected())
tbl_fit = tbl_obs[(tbl_obs['obsdate']>phase_min) & (tbl_obs['obsdate']<phase_max)]
tbl_filt = observed_data.get_filt_data(tbl_fit)

x_fit = []
y_fit = []
e_y_fit = []
fit_params = Parameters()
fit_params.add('exptime', value = 59529, min = 59527.3002, max = 59535)
for filter_ in fit_filterset:
    if filter_ in tbl_filt.keys():
        #fit_params.add(f'exptime_{filter_}', value = 59529, min = 59525, max = 59535)
        fit_params.add(f'amp_{filter_}', value = 1000, min = 0, max = 200000)
        fit_params.add(f'alpha_{filter_}', value = 2, min = 1, max = 10)
        obsdate_filt = tbl_filt[filter_]['obsdate']
        flux_filt = mag_to_flux(tbl_filt[filter_]['mag'], zp = 25)
        e_flux_filt = flux_filt * 2.303/2.5 * tbl_filt[filter_]['e_mag']
        x_fit.append(obsdate_filt.tolist())
        y_fit.append(flux_filt.tolist())
        e_y_fit.append(e_flux_filt.tolist())
out = minimize(calc_chisq, fit_params, args = (x_fit, y_fit, e_y_fit, fit_filterset))
report_fit(out.params)
#%% Visualization
color_key, offset_key, _, _, label_key = load_filt_keys(fit_filterset)
plt.figure(dpi = 300, figsize = (5, 8))
plt.gca().invert_yaxis()
phase_range = np.arange(out.params[f'exptime'].value, 59540, 0.1)

for filter_ in fit_filterset:
    exptime = out.params[f'exptime']
    #exptime = out.params[f'exptime_{filter_}']
    amp = out.params[f'amp_{filter_}']
    alpha= out.params[f'alpha_{filter_}']
    tbl_filter = tbl_filt_all[filter_]
    tbl_filter.sort('obsdate')
    flux_model = fireball_filter(out.params, phase_range, filter_)
    mag_model = flux_to_mag(flux_model, zp = 25)
    #plt.text(59531,10,'%.4f'%exptime)
    #plt.scatter(tbl_filter['obsdate'], tbl_filter['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_])
    plt.plot(phase_range, mag_model + offset_key[filter_], c = color_key[filter_], label = rf'[{label_key[filter_]}] $\alpha = {round(alpha.value,2)}$', linestyle= '--', linewidth = 1)
    if filter_ == 'U':
        mag_U = mag_model
    if filter_ == 'B':
        mag_B = mag_model
    if filter_ == 'V':
        mag_V = mag_model
    if filter_ == 'g':
        mag_g = mag_model
    if filter_ == 'r':
        mag_r = mag_model
        
#plt.axvline(x = out.params['exptime'].value, linestyle= '--')
#plt.plot(phase_range, mag_U - mag_B, c = 'k', label = rf'[{label_key[filter_]}] $\alpha = {round(alpha.value,2)}$', linestyle= '--', linewidth = 1)
#plt.plot(phase_range, mag_B - mag_V, c = 'b', label = rf'[{label_key[filter_]}] $\alpha = {round(alpha.value,2)}$', linestyle= '--', linewidth = 1)
#plt.plot(phase_range, mag_g - mag_r, c = 'g', label = rf'[{label_key[filter_]}] $\alpha = {round(alpha.value,2)}$', linestyle= '--', linewidth = 1)
#plt.ylim(-1, 1.7)
#plt.legend(loc = 4)
observed_data.show_lightcurve( day_binsize = 5, color_BV = False, color_gr = False, color_UB = False, UL = True, label = False, label_location=2, scatter_size= 120)
plt.xlim(phase_range[0]-3, 59540)
plt.ylim(24, 8)

# %%
