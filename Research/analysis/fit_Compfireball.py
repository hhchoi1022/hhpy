#%%
#import sys
#sys.path.append('/data7/yunyi/temp_supernova/Gitrepo/')
import os
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from Research.analysis.observedphot import ObservedPhot
from Research.helper import Helper
from observedphot import ObservedPhot
from lmfit import Parameters, minimize
#%% Observation
helper = Helper()
DM = 31.18
ZP = 25
filepath_all = '/data1/supernova_rawdata/SN2021aefx/photometry/all_phot_MW_dereddening_Host_dereddening.dat'
model_directory = '/data1/supernova_model/Comp_model'
#filepath_all = '/data7/yunyi/temp_supernova/Gitrepo/Research/analysis/all_phot_MW_dereddening_Host_dereddening.dat'
#model_directory = '/data7/yunyi/temp_supernova/DOM_model'

fit_filterset = 'UBVugri'
fit_start_mjd : int = 59529
fit_end_mjd : int = 59537

# Query data
obs_tbl = ascii.read(filepath_all, format = 'fixed_width')

# Exclude some observatories with poor photometry
observed_data = ObservedPhot(obs_tbl)
observed_data.exclude_observatory(['LasCumbres0.4m', 'Swift'])
observed_data.exclude_filter(['Unfilte'])

# Get detected data
detected_tbl = observed_data.get_data_detected()

# Construct table for fitting
fit_idx = [filter_ in fit_filterset for filter_ in detected_tbl['filter']]
fit_tbl = detected_tbl[fit_idx]
fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
fit_tbl.sort('obsdate')

# Add systematic error
systematic_error_cut = 0.03
def adjust_errors(errors, systematic_error):
    return np.sqrt(errors**2 + systematic_error**2)
obs_lowerror_tbl = fit_tbl[fit_tbl['e_mag'] < systematic_error_cut]
adjusted_errors = np.round(adjust_errors(obs_lowerror_tbl['e_mag'], systematic_error= systematic_error_cut),3)
obs_lowerror_tbl['e_mag'] = adjusted_errors
fit_tbl[fit_tbl['e_mag'] < systematic_error_cut] = obs_lowerror_tbl
fit_tbl['e_mag'] = np.round(fit_tbl['e_mag'], 3)

# Add flux 
fit_tbl['flux'] = helper.mag_to_flux(fit_tbl['mag'])
fit_tbl['e_flux'] = fit_tbl['e_mag']*helper.mag_to_flux(fit_tbl['mag'])*2.303/2.5
fit_tbl['absmag'] = (fit_tbl['mag'] - DM).round(3)

#%%
# Visualize
plt.figure(dpi = 400, figsize = (8, 8))
ax1, ax2 = observed_data.show_lightcurve(day_binsize = 5,
                              scatter_linewidth=0.5, 
                              scatter_size=50, 
                              errorbar_linewidth=0.5, 
                              errorbar_capsize=0.1, 
                              color_UB = True,
                              color_BV = True, 
                              color_gr = True, 
                              UL = True, 
                              UL_alpha = 0.8,
                              label = True, 
                              label_location=0, 
                              )
ax1.fill_betweenx(y = [ 30, 0], x1 = fit_start_mjd, x2 = fit_end_mjd, color = 'gray', alpha = 0.2)
ax1.set_ylim(22, 6)
plt.xlim(59525, 59545)
#%%
from Research.model import CompanionInteractionK10

def fireball_model(time, 
                   amplitude, 
                   exptime, 
                   alpha):
    flux = amplitude * (time - exptime )**alpha
    np.nan_to_num(flux, copy=False, nan=0.001)
    return flux

def get_CEI_spline(model_CEI,
                   exptime_CEI,
                   filterset : str = 'UBVRIugri',
                   smooth : float = 0.05):
    spl_dict = dict()
    for filter_ in filterset:
        model_mag = model_CEI[filter_]
        inf_idx = np.isinf(model_mag)
        mag_CEI = model_mag[~inf_idx]
        phase_CEI = model_CEI['phase'][~inf_idx]
        spl, _ = helper.interpolate_spline(phase_CEI + exptime_CEI, mag_CEI, show = False, smooth = smooth)
        spl_dict[filter_] = spl
    return spl_dict
#%%
def fit_both(fit_tbl,
             rstar,
             m_wd,
             v9 : float = 1.0,
             
             fit_method : str = 'leastsq'
             ):

    def chisq_both(params, x_fit, y_fit, e_y_fit, filter_key, model_CEI):
        vals = params.valuesdict()
        exptime_CEI = vals['exptime_CEI']
        exptime_FB = vals['exptime_FB']
        chisq_allfilter = []
        spl_allfilt_CEI = get_CEI_spline(model_CEI = model_CEI, exptime_CEI = exptime_CEI, filterset = filter_key)
        for mjd, obs_flux, obs_fluxerr, filter_ in zip(x_fit, y_fit, e_y_fit, filter_key):
            spl_CEI = spl_allfilt_CEI[filter_]
            fireball_alpha = vals[f'alpha_{filter_}']
            fireball_amplitude = vals[f'amplitude_{filter_}']
            fireball_flux = fireball_model(time = mjd, amplitude = fireball_amplitude, alpha = fireball_alpha, exptime = exptime_FB)
            CEI_flux = helper.mag_to_flux(spl_CEI(mjd)+DM)
            both_flux = fireball_flux + CEI_flux
            chisq_singlefilter = (((obs_flux - both_flux)/obs_fluxerr)**2)
            chisq_allfilter.append(chisq_singlefilter)
        return np.concatenate(chisq_allfilter)

    # Input
    filter_tbls = fit_tbl.group_by('filter').groups
    filter_key = filter_tbls.keys['filter']
    fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
    x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
    y_fit = [np.array((fit_table[filter_]['flux'].tolist())) for filter_ in filter_key]
    e_y_fit = [np.array((fit_table[filter_]['e_flux'].tolist())) for filter_ in filter_key]

    # Parameters
    fit_params_CEI_FB = Parameters()
    fit_params_CEI_FB.add('exptime_CEI', value = 59528.5, min = 59525, max = 59535)
    fit_params_CEI_FB.add('exptime_FB', value = 59528.5, min = 59525, max = 59535)
    for filter_ in filter_key:
        fit_params_CEI_FB.add(f'alpha_{filter_}', value = 2, min = 0, max = 4)
        fit_params_CEI_FB.add(f'amplitude_{filter_}', value = 3000, min = 1, max = 500000)
    
    # Fitting
    t_range = np.arange(0.1, 10, 0.1)
    Comp_model = CompanionInteractionK10(rstar = rstar, m_wd = m_wd, v9 = v9)
    model_CEI = Comp_model.get_LC(td = t_range, filterset = ''.join(filter_key), search_directory = model_directory, save = True)
    out = minimize(chisq_both, fit_params_CEI_FB, args = (x_fit, y_fit, e_y_fit, filter_key, model_CEI), method = fit_method)
    return out

#%%
header_parameters = ['rstar','m_wd','v9']
header_fitvalues = ['exptime_CEI', 'exptime_FB']
header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
for filter_ in fit_filterset:
    header_fitvalues.append(f'alpha_{filter_}')
    header_fitvalues.append(f'amplitude_{filter_}')
tot_header = header_parameters + header_fitvalues + header_fitconfig
result_tbl = Table(names = tot_header)
result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))
i = 0
sample_parameters = dict(rstar  = 2.1,
                         m_wd = 1.2,
                         v9 = 1.0)
result = fit_both(fit_tbl = fit_tbl, **sample_parameters, fit_method = 'leastsq')
data_parameters = sample_parameters
data_fitvalues = result.params.valuesdict()
data_fitconfig = dict(success = result.success, nfev = result.nfev, ndata = result.ndata, nvar = result.nvarys, chisq = result.chisqr, redchisqr = result.redchi, aic = result.aic, bic = result.bic)
all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
all_values = []
for colname in result_tbl.columns:
    value = all_data[colname]
    all_values.append(value)
result_tbl.add_row(vals = all_values)
#%%
#%%
header_parameters =['rstar','m_wd','v9']
header_fitvalues = ['exptime_CEI', 'exptime_FB']
header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
for filter_ in fit_filterset:
    header_fitvalues.append(f'alpha_{filter_}')
    header_fitvalues.append(f'amplitude_{filter_}')
tot_header = header_parameters + header_fitvalues + header_fitconfig
result_tbl = Table(names = tot_header)
result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))

range_rstar = np.round(np.arange(0.5, 30, 0.05),2)
range_m_wd =  np.round(np.arange(1.0, 1.5, 0.1),1)
range_v9 = np.round(np.arange(0.7, 1.4, 0.1),1)
os.makedirs('./result/Comp_fit_result', exist_ok = True)

for rstar in range_rstar:
    for m_wd in range_m_wd:
        result_tbl.write(f'./result/Comp_fit_result/Timestamp_M_wd_{m_wd}', format = 'ascii.fixed_width', overwrite = True)
        for v9 in range_v9:
            try:
                result = fit_both(rstar,
                                  m_wd,
                                  v9,
                                  fit_method = 'leastsq'
                                  )
                data_parameters = dict(rstar = rstar, m_wd = m_wd, v9 = v9)
                data_fitvalues = result.params.valuesdict()
                data_fitconfig = dict(success = result.success, nfev = result.nfev, ndata = result.ndata, nvar = result.nvarys, chisq = result.chisqr, redchisqr = result.redchi, aic = result.aic, bic = result.bic)
                all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
                all_values = []
                for colname in result_tbl.columns:
                    value = all_data[colname]
                    all_values.append(value)
                result_tbl.add_row(vals = all_values)
            except:
                data_parameters = dict(rstar = rstar, m_wd = m_wd)
                data_fitvalues = {value : 99999 for value in header_fitvalues}
                data_fitconfig = dict(success = False, nfev = 99999, ndata = 99999, nvar = 99999, chisq = 99999, redchisqr = 99999, aic = 99999, bic = 99999)
                all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
                all_values = []
                for colname in result_tbl.columns:
                    value = all_data[colname]
                    all_values.append(value)
                result_tbl.add_row(vals = all_values)   
result_tbl.remove_row(index = 0)

#%%
import numpy as np
import os
import time
import multiprocessing as mp

range_rstar = np.round(np.arange(0.5, 30, 0.05),2)
range_m_wd =  np.round(np.arange(1.0, 1.5, 0.1),1)
range_v9 = np.round(np.arange(0.7, 1.4, 0.1),1)

def process_combination(args):
    rstar, m_wd, v9, fit_tbl = args
    header_parameters =['rstar','m_wd','v9']
    header_fitvalues = ['exptime_CEI', 'exptime_FB']
    header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
    for filter_ in fit_filterset:
        header_fitvalues.append(f'alpha_{filter_}')
        header_fitvalues.append(f'amplitude_{filter_}')
    tot_header = header_parameters + header_fitvalues + header_fitconfig
    result_tbl = Table(names = tot_header)
    result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))
    try:
        print('Start calculation: ', rstar, m_wd, v9)
        result = fit_both(fit_tbl=fit_tbl,
                          rstar = rstar,
                          m_wd = m_wd,
                          v9 = v9,
                          fit_method='leastsq'
                          )
        data_parameters = dict(rstar = rstar, m_wd = m_wd, v9 = v9)
        data_fitvalues = result.params.valuesdict()
        data_fitconfig = dict(success = result.success, nfev = result.nfev, ndata = result.ndata, nvar = result.nvarys, chisq = result.chisqr, redchisqr = result.redchi, aic = result.aic, bic = result.bic)
        all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
        all_values = [all_data[colname] for colname in result_tbl.columns]
        result_tbl.add_row(vals=all_values)
    except:
        data_parameters = dict(rstar = rstar, m_wd = m_wd)
        data_fitvalues = {value : 99999 for value in header_fitvalues}
        data_fitconfig = dict(success = False, nfev = 99999, ndata = 99999, nvar = 99999, chisq = 99999, redchisqr = 99999, aic = 99999, bic = 99999)
        all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
        all_values = [all_data[colname] for colname in result_tbl.columns]
        result_tbl.add_row(vals=all_values)
    #os.makedirs(f'/data7/yunyi/temp_supernova/result/Comp_fit_result/M{m_wd}', exist_ok = True)
    os.makedirs(f'/data1/supernova_model/result/Comp_fit_result/M{m_wd}', exist_ok = True)
    result_tbl.remove_row(index = 0)
    #result_tbl.write(f'/data7/yunyi/temp_supernova/result/Comp_fit_result/M{m_wd}/Rstar_{rstar}_V9_{v9}.txt', format = 'ascii.fixed_width', overwrite = True)
    result_tbl.write(f'/data1/supernova_model/result/Comp_fit_result/M%.1f/%.2f_%.1f_%.1f.fit'%(m_wd, rstar, m_wd, v9), format = 'ascii.fixed_width', overwrite = True)

def main(fit_tbl):
    os.makedirs(f'/data1/supernova_model/result/Comp_fit_result', exist_ok = True)
    #os.makedirs(f'/data7/yunyi/temp_supernova/result/Comp_fit_result', exist_ok=True)
    
    # Prepare the list of all combinations of parameters
    all_combinations = [(rstar, m_wd, v9, fit_tbl)
                        for rstar in range_rstar
                        for m_wd in range_m_wd
                        for v9 in range_v9
                        ]
    # Use multiprocessing to process the combinations in parallel
    with mp.Pool(processes=8) as pool:
        pool.map(process_combination, all_combinations)


if __name__ == '__main__':
    pass
    #main(fit_tbl = fit_tbl)
#%%
'''
import glob
from astropy.table import vstack
result_key = '/data1/supernova_model/result/Comp_fit_result/*/*.fit'
files = glob.glob(result_key)
result_tbl = Table()
for file_ in files:
    tbl = ascii.read(file_, format = 'fixed_width')
    result_tbl = vstack([result_tbl, tbl])
result_tbl.write('/data1/supernova_model/result/Comp_fit_result.fit', format = 'ascii.fixed_width', overwrite = True)
'''
#%%
result_tbl = ascii.read('/data1/supernova_model/result/Comp_fit_result.fit', format = 'fixed_width')
fit_filterset = set(fit_tbl['filter'])
#fit_filterset = 'UBgVri'
i = 2
result_values = result_tbl[i]
exptime_CEI = result_values['exptime_FB']
exptime_FB = result_values['exptime_CEI']
filter_key = fit_tbl.group_by('filter').groups.keys['filter']
color_key, offset_key, _, _, label_key = helper.load_filt_keys()
plt.figure(dpi = 300, figsize = (4.5, 6.5))
plt.gca().invert_yaxis()
phase_min_FB = np.max([59526, result_values['exptime_FB']])
phase_min_CEI = np.max([59526, result_values['exptime_CEI']])
phase_range_FB = np.arange(phase_min_FB, 59540, 0.1)
phase_range_CEI = np.arange(phase_min_CEI, 59540, 0.1)
phase_range_CEI = np.arange(np.min([phase_min_FB, phase_min_CEI]), 59540, 0.1)


CEI_model = CompanionInteractionK10(rstar = result_values['rstar'], m_wd = result_values['m_wd'], v9 = result_values['v9'])
CEI_LC = CEI_model.get_LC(td = phase_range_CEI, filterset = 'UBVRIugri', search_directory = model_directory, save = True)
spl_allfilt_CEI = get_CEI_spline(CEI_LC, exptime_CEI = result_values['exptime_CEI'], filterset = ''.join(filter_key))

tbl_UL = observed_data.get_data_ul()
tbl_obs = observed_data.get_data_detected()
ax1, ax2 = observed_data.show_lightcurve(day_binsize = 5,
                            scatter_linewidth=0.5, 
                            scatter_size=50, 
                            scatter_alpha = 0.2,
                            errorbar_linewidth=0.5, 
                            errorbar_capsize=0.1, 
                            color_UB = True,
                            color_BV = True, 
                            color_gr = True, 
                            UL = True, 
                            UL_alpha = 0.8,
                            label = True, 
                            label_location=4, 
                            )
for filter_ in fit_filterset:
    amp = result_values[f'amplitude_{filter_}']
    alpha= result_values[f'alpha_{filter_}']
    flux_FB = fireball_model(time = phase_range_CEI, amplitude = amp, alpha = alpha, exptime = result_values['exptime_FB'])
    spl_CEI = spl_allfilt_CEI[filter_]
    flux_CEI = helper.mag_to_flux(spl_CEI(phase_range_CEI)+DM)
    flux_both = flux_FB+ flux_CEI
    mag_model = helper.flux_to_mag(flux_FB, zp = ZP)
    mag_DOM = helper.flux_to_mag(flux_CEI, zp = ZP)
    mag_both = helper.flux_to_mag(flux_both, zp = ZP)
    tbl_UL_filter = tbl_UL[tbl_UL['filter'] == filter_]
    tbl_obs_filter = tbl_obs[tbl_obs['filter'] == filter_]
    ax1.plot(phase_range_CEI, mag_model + offset_key[filter_], c = color_key[filter_], label = rf'[{label_key[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= ':', linewidth = 1, alpha = 0.4)
    ax1.plot(phase_range_CEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= '--', linewidth = 1, alpha = 0.4)
    ax1.plot(phase_range_CEI, mag_both + offset_key[filter_], c = color_key[filter_], linestyle= '-', linewidth = 1, alpha = 1)
    # For color plot 
    if filter_ == 'U':
        mag_U_model = mag_model
        mag_U_CEI = mag_DOM
        mag_U_both = mag_both
    if filter_ == 'B':
        mag_B_model = mag_model
        mag_B_CEI = mag_DOM
        mag_B_both = mag_both
    if filter_ == 'V':
        mag_V_model = mag_model
        mag_V_CEI = mag_DOM
        mag_V_both = mag_both
    if filter_ == 'g':
        mag_g_model = mag_model
        mag_g_CEI = mag_DOM
        mag_g_both = mag_both
    if filter_ == 'r':
        mag_r_model = mag_model
        mag_r_CEI = mag_DOM
        mag_r_both = mag_both

#ax1.plot(0,0, '--', c='k', label = 'Power law')
#ax1.plot(0,0, ':', c='k', label = 'Compansion-ejecta')
#ax1.plot(0,0, '-', c='k', label = 'Power law + Companion-ejecta')
#ax1.legend(loc = 3)
ax2.plot(phase_range_CEI, mag_U_model - mag_B_model -0.5, c = 'cyan', label = 'U-B', linestyle= ':', linewidth = 1, alpha = 0.4)
ax2.plot(phase_range_CEI, mag_B_model - mag_V_model + 0.5, c = 'b', label = 'B-V', linestyle= ':', linewidth = 1, alpha = 0.4)
ax2.plot(phase_range_CEI, mag_g_model - mag_r_model, c = 'g', label = 'g-r', linestyle= ':', linewidth = 1, alpha = 0.4)
ax2.plot(phase_range_CEI, mag_U_CEI - mag_B_CEI -0.5, c = 'cyan', label = 'U-B', linestyle= '--', linewidth = 1, alpha = 0.4)
ax2.plot(phase_range_CEI, mag_B_CEI - mag_V_CEI +0.5, c = 'b', label = 'B-V', linestyle= '--', linewidth = 1, alpha = 0.4)
ax2.plot(phase_range_CEI, mag_g_CEI - mag_r_CEI, c = 'g', label = 'g-r', linestyle= '--', linewidth = 1, alpha = 0.4)
ax2.plot(phase_range_CEI, mag_U_both - mag_B_both -0.5, c = 'cyan', label = 'U-B', linestyle= '-', linewidth = 1, alpha = 1)
ax2.plot(phase_range_CEI, mag_B_both - mag_V_both +0.5, c = 'b', label = 'B-V', linestyle= '-', linewidth = 1, alpha = 1)
ax2.plot(phase_range_CEI, mag_g_both - mag_r_both, c = 'g', label = 'g-r', linestyle= '-', linewidth = 1, alpha = 1)
#ax2.yli(np.min(obs_tbl['obsdate']) + np.arange(-20, 200, 5), np.arange(-20, 200, 5) )
ax1.set_xlim(phase_range_FB[0]-1, 59537)
ax2.set_xlim(phase_range_FB[0]-1, 59537)

ax1.set_ylim(22.5, 8)

ax1.clear()
show_idx = [0,3,7, 9, 10, 12]
obs_spec = ascii.read('/data1/supernova_rawdata/SN2021aefx/photometry/all_spec_MW_dereddening_Host_dereddening.dat', format = 'fixed_width')
obs_spec_phot = ObservedPhot(data_tbl  = obs_spec)
obs_spec_phot.data.sort('obsdate')
filt_spec_tbl = obs_spec_phot.get_filt_data(obs_spec_phot.data)
UB_tbl = helper.match_table(filt_spec_tbl['U'], filt_spec_tbl['B'], key = 'obsdate')
BV_tbl = helper.match_table(filt_spec_tbl['B'], filt_spec_tbl['V'], key = 'obsdate')
gr_tbl = helper.match_table(filt_spec_tbl['g'], filt_spec_tbl['r'], key = 'obsdate')

import glob
import matplotlib.cm as cm  # Import the colormap
from Research.spectroscopy.spectroscopyfile import SpectroscopyFile
from Research.spectroscopy import Spectrum
from Research.spectroscopy import TimeSeriesSpectrum
from Research.model import CompanionInteractionK10

num_files = len(show_idx)  # Determine the number of files
colormap = cm.get_cmap('cool', num_files)  # Choose a colormap and set the number of colors
bb_temp = [11000, 6000, 8000, 8000, 9000, 10000]
CEI_model = CompanionInteractionK10(rstar = 2.35, m_wd = 1.3, v9 = 0.9).get_LC(td = np.arange(0.1, 10, 0.1))
exptime_CEI = 59528.97627959757
CEI_model['phase'] = CEI_model['phase'] + exptime_CEI
spl_temp,_ = helper.interpolate_spline(list(CEI_model['phase']), list(CEI_model['Temperature_eff']), show = False)
position_txt = [3, 0.3]
for i, idx in enumerate(show_idx):
    file_ = UB_tbl[idx]['filename_1']
    obsdate = UB_tbl[idx]['obsdate_1']
    specfile = SpectroscopyFile(file_)
    flux = specfile.flux
    wl = specfile.wavelength
    spec = Spectrum(wl, flux, flux_unit='flamb')
        
    # Get a color from the colormap
    color = colormap(i)
    Planck = helper.planck
    temp = spl_temp(obsdate)
    val = Planck(temperature=temp, wl_AA=wl)
    spec_bb = Spectrum(wl, val['flamb'], flux_unit='flamb')
    
    # If spec.show() supports a color parameter
    spec.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.05, normalize_cenwl=7500, color=color, label = specfile.obsdate, offset = -2*i, axis = ax1, linewidth = 1)
    if i < 2:
        spec_bb.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.05, normalize_cenwl=7500, color='black', offset = -2*i, axis = ax1, linestyle = '--', linewidth = 0.5)
        ax1.text(5000, position_txt[i], '%.0f K'%temp)
    ax2.scatter(BV_tbl['obsdate_1'][idx], BV_tbl['mag_1'][idx] - BV_tbl['mag_2'][idx] + 0.5, facecolor = 'b', edgecolor = color, marker = '*', s = 150, alpha = 1, zorder = 10)
    ax2.scatter(gr_tbl['obsdate_1'][idx], gr_tbl['mag_1'][idx] - gr_tbl['mag_2'][idx], facecolor = 'g', edgecolor = color, marker = '*', s = 150, alpha = 1, zorder = 10)

ax1.tick_params(axis='x', which='both', direction='in', top=True)
ax2.tick_params(axis='x', which='both', direction='in', top=True)
ax1.set_ylim(-10, 5.5)
ax1.set_xlim(3000, 10000)
ax1.set_xticks(np.arange(3000, 11000, 1000), np.arange(3000, 11000, 1000))
# Move x-ticks to the top for ax1
ax1.xaxis.tick_top()  # This moves the x-tick labels to the top of ax1
ax1.xaxis.set_label_position('top')  # This moves the x-axis label to the top
ax1.set_ylabel(rf'Normalized flux ($F_\lambda$) + offset')
ax2.set_xlim(59528, 59537)
# %%
plt.figure(dpi = 300, figsize = (6,4))

plt.plot(CEI_model['phase'], CEI_model['Temperature_eff'])
plt.plot([-1,0.8855], [11391, 11391], color = 'r', linestyle = '--')
plt.plot([0.8855,0.8855], [-3000, 11391], color = 'r', linestyle = '--')
#plt.plot([0.08855,0], [0.8855, 11391], color = 'r', linestyle = '--')
plt.xlim(-1,10)
plt.ylim(-2000, 35000)
#%%
plt.figure(dpi = 300, figsize = (6,4))
obsdate = UB_tbl[0]['obsdate_1'] - exptime_CEI
CEI_model = CompanionInteractionK10(rstar = result_values['rstar'], m_wd = result_values['m_wd'], v9 = result_values['v9']).get_LC(td = np.arange(0.1, 10, 0.1))
plt.plot(CEI_model['phase'], CEI_model['Temperature_eff'], color = 'k')
temp = spl_temp(obsdate + exptime_CEI)
plt.plot([-1,obsdate], [temp, temp], color = 'k', linestyle = '--')
plt.plot([obsdate,obsdate], [-3000, temp], color = 'k', linestyle = '--')
plt.xlabel('Days since the first detection [MJD - 59529.3318]')
plt.ylabel('Effective temperature [K]')
plt.xticks(np.array([0,2,4,6,8,10])+59529.3318-exptime_CEI, np.array([0,2,4,6,8,10]))
#plt.plot([0.08855,0], [0.8855, 11391], color = 'r', linestyle = '--')
plt.xlim(-1,10)
plt.ylim(-2000, 35000)

#%%

fig = plt.figure(dpi = 300)
mag_tbl = CEI_model
ax1 = plt.subplot()
ax1.plot(mag_tbl['phase'], mag_tbl['Luminosity_shock'], c='k')
ax1.set_yscale('log')
ax1.set_xticks(np.array([0,2,4,6,8,10])+59529.3318-exptime_CEI, np.array([0,2,4,6,8,10]))

ax1.plot([-1,obsdate], [2.8e41, 2.8e41], color = 'k', linestyle = '--')
ax1.set_yticks([1e41, 5e41, 1e42, 5e42, 1e43, 5e43], [1e41, 5e41, 1e42, 5e42, 1e43, 5e43])
ax1.set_ylabel(r'$L_{shock}\ [erg/s]$', fontsize = 10)
ax1.set_xlabel('Phase [day]')
#ax1.set_ylim(5e40, 5e43)

ax2 = ax1.twinx()
ax2.plot([10,obsdate], [temp, temp], color = 'r', linestyle = '--')
ax2.plot([obsdate,obsdate], [-3000, temp], color = 'k', linestyle = '-')
ax2.plot(mag_tbl['phase'], mag_tbl['Temperature_eff'], c='r')
ax2.set_ylabel(r'$T_{eff}\ [K]$', rotation = 270, c= 'r')
ax2.set_yticks([5000,10000,15000,20000,25000,30000,35000], [5000,10000,15000,20000,25000,30000,35000], c= 'r')
ax2.set_ylim(0, 35000)

plt.plot(1, 1, c='k', label =r'$L_{shock}$')
plt.plot(1, 1, c='r', label =r'$T_{eff}$')
plt.legend()
plt.xlim(-1,10)
plt.ylim(-2000, 35000)

# %%
