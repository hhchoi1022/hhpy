#%%
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from observedphot import ObservedPhot
from lmfit import Parameters, minimize
from scipy.interpolate import UnivariateSpline
from Research.model import DOMInteractionL17
from Research.helper import Helper
import matplotlib
%matplotlib inline
#%% Observation
helper = Helper()
DM = 31.18
ZP = 25
filepath_all = '/data1/supernova_rawdata/SN2021aefx/photometry/all_phot_MW_dereddening_Host_dereddening.dat'
model_directory = '/data1/supernova_model/DOM_model/'
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
#%% Define model
def fireball_model(time, 
                   amplitude, 
                   exptime, 
                   alpha):
    flux = amplitude * (time - exptime )**alpha
    np.nan_to_num(flux, copy=False, nan=0.0001)
    return flux

def get_DEI_spline(model_DEI,
                   exptime_DEI,
                   filterset : str = 'UBVRIugri',
                   smooth : float = 0.05):
    spl_dict = dict()
    for filter_ in filterset:
        model_mag = model_DEI[filter_]
        inf_idx = np.isinf(model_mag)
        mag_DEI = model_mag[~inf_idx]
        phase_DEI = model_DEI['phase'][~inf_idx]
        spl, _ = helper.interpolate_spline(phase_DEI + exptime_DEI, mag_DEI, show = False, smooth = smooth)
        spl_dict[filter_] = spl
    return spl_dict
#%%
def fit_both(fit_tbl,
             E_exp,
             M_ej,
             kappa,
             M_dom,
             V_dom,
             f_dom,
             t_delay,
             f_comp,
             
             fit_method = 'leastsq'
             ):

    def chisq_both(params, x_fit, y_fit, e_y_fit, filter_key, model_DEI):
        vals = params.valuesdict()
        exptime_DEI = vals['exptime_DEI']
        exptime_FB = vals['exptime_FB']
        chisq_allfilter = []
        spl_allfilt_DEI = get_DEI_spline(model_DEI = model_DEI, exptime_DEI = exptime_DEI, filterset = filter_key)
        for mjd, obs_flux, obs_fluxerr, filter_ in zip(x_fit, y_fit, e_y_fit, filter_key):
            spl_DEI = spl_allfilt_DEI[filter_]
            fireball_alpha = vals[f'alpha_{filter_}']
            fireball_amplitude = vals[f'amplitude_{filter_}']
            fireball_flux = fireball_model(time = mjd, amplitude = fireball_amplitude, alpha = fireball_alpha, exptime = exptime_FB)
            DEI_flux = helper.mag_to_flux(spl_DEI(mjd)+DM)
            both_flux = fireball_flux + DEI_flux
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
    fit_params_DEI_FB = Parameters()
    fit_params_DEI_FB.add('exptime_DEI', value = 59528.5, min = 59525, max = 59535)
    fit_params_DEI_FB.add('exptime_FB', value = 59528.5, min = 59525, max = 59535)
    for filter_ in filter_key:
        fit_params_DEI_FB.add(f'alpha_{filter_}', value = 2, min = 0, max = 4)
        fit_params_DEI_FB.add(f'amplitude_{filter_}', value = 3000, min = 1, max = 500000)
    
    # Fitting
    t_range = np.arange(0.1, 10, 0.1)
    DOM_model = DOMInteractionL17(E_exp = E_exp, M_ej = M_ej, kappa = kappa, M_dom = M_dom, V_dom = V_dom, f_dom = f_dom, t_delay = t_delay, f_comp = f_comp)
    model_DEI = DOM_model.get_LC(td = t_range, filterset = ''.join(filter_key), search_directory = model_directory, save = True)
    out = minimize(chisq_both, fit_params_DEI_FB, args = (x_fit, y_fit, e_y_fit, filter_key, model_DEI), method = fit_method)
    return out

#%%
header_parameters = ['E_exp','M_ej','kappa','t_delay','f_comp','M_dom','V_dom','f_dom']
header_fitvalues = ['exptime_DEI', 'exptime_FB']
header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
for filter_ in fit_filterset:
    header_fitvalues.append(f'alpha_{filter_}')
    header_fitvalues.append(f'amplitude_{filter_}')
tot_header = header_parameters + header_fitvalues + header_fitconfig
result_tbl = Table(names = tot_header)
result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))
i = 0
sample_parameters = dict(E_exp = 1.4,
                         M_ej = 1.0,
                         kappa = 0.05,
                        
                         M_dom = 0.11,
                         V_dom = 5e3,
                         f_dom = 0.14,
                        
                         t_delay = 200,
                         f_comp = 1.5)
import time
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
header_parameters = ['E_exp','M_ej','kappa','t_delay','f_comp','M_dom','V_dom','f_dom']
header_fitvalues = ['exptime_DEI', 'exptime_FB']
header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
for filter_ in fit_filterset:
    header_fitvalues.append(f'alpha_{filter_}')
    header_fitvalues.append(f'amplitude_{filter_}')
tot_header = header_parameters + header_fitvalues + header_fitconfig
result_tbl = Table(names = tot_header)
result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))
#%%

import numpy as np

# Define the ranges with rounding or conversion to integers where necessary
range_E_exp = np.round(np.arange(0.8, 1.6, 0.2), 2)  # 10^51 ergs, rounded to 2 decimal places
range_M_ej = np.round(np.arange(0.8, 1.6, 0.2), 2)   # solar mass, rounded to 2 decimal places
range_kappa = np.round([0.03, 0.05], 2)              # cm^2/g, rounded to 2 decimal places
range_t_delay = np.round(np.arange(1e2, 1e4, 200), 0)  # s, rounded to 0 decimal places (integer-like)
range_f_comp = np.round([1.5], 2)                    # compress fraction, rounded to 2 decimal places
range_M_dom = np.round(np.arange(0.08, 0.16, 0.02), 2)  # solar mass, rounded to 2 decimal places
range_v_dom = np.int_(np.round(np.arange(5e3, 10e3, 1e3), 0))  # km/s, rounded and converted to integer
range_f_dom = np.round(np.arange(0.03, 0.16, 0.02), 2)  # fraction of DOM mass, rounded to 2 decimal places

import time
start = time.time()
os.makedirs('./result/DOM_fit_result', exist_ok = True)
for E_exp in range_E_exp:
    result_tbl.write(f'./result/DOM_fit_result/Timestamp_E_exp_{E_exp}', format = 'ascii.fixed_width', overwrite = True)
    for M_ej in range_M_ej:
        result_tbl.write(f'./result/DOM_fit_result/Timestamp_M_ej_{M_ej}', format = 'ascii.fixed_width', overwrite = True)
        for kappa in range_kappa:
            for t_delay in range_t_delay:
                result_tbl.write(f'./result/DOM_fit_result/Timestamp_t_delay_{t_delay}', format = 'ascii.fixed_width', overwrite = True)
                for f_comp in range_f_comp:
                    for M_dom in range_M_dom:
                        for V_dom in range_v_dom:
                            for f_dom in range_f_dom:
                                try:
                                    result = fit_both(fit_tbl = fit_tbl,
                                                      E_exp = E_exp,
                                                      M_ej = M_ej,
                                                      kappa = kappa,
                                                      M_dom = M_dom,
                                                      V_dom = V_dom,
                                                      f_dom = f_dom,
                                                      t_delay = t_delay,
                                                      f_comp = f_comp,
                                                      fit_method = 'leastsq'
                                                      )
                                    data_parameters = dict(E_exp = E_exp, M_ej = M_ej, kappa = kappa, t_delay = t_delay, f_comp = f_comp, M_dom = M_dom, V_dom = V_dom, f_dom = f_dom)
                                    data_fitvalues = result.params.valuesdict()
                                    data_fitconfig = dict(success = result.success, nfev = result.nfev, ndata = result.ndata, nvar = result.nvarys, chisq = result.chisqr, redchisqr = result.redchi, aic = result.aic, bic = result.bic)
                                    all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
                                    all_values = []
                                    for colname in result_tbl.columns:
                                        value = all_data[colname]
                                        all_values.append(value)
                                    result_tbl.add_row(vals = all_values)
                                except:
                                    data_parameters = dict(E_exp = E_exp, M_ej = M_ej, kappa = kappa, t_delay = t_delay, f_comp = f_comp, M_dom = M_dom, V_dom = V_dom, f_dom = f_dom)
                                    data_fitvalues = {value : 99999 for value in header_fitvalues}
                                    data_fitconfig = dict(success = False, nfev = 99999, ndata = 99999, nvar = 99999, chisq = 99999, redchisqr = 99999, aic = 99999, bic = 99999)
                                    all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
                                    all_values = []
                                    for colname in result_tbl.columns:
                                        value = all_data[colname]
                                        all_values.append(value)
                                    result_tbl.add_row(vals = all_values)
result_tbl.remove_row(index = 0)
print(time.time() - start)
#%%
import numpy as np
import os
import time
import multiprocessing as mp

# Define the ranges with rounding or conversion to integers where necessary
range_E_exp = np.round(np.arange(0.8, 1.6, 0.2), 2)  # 10^51 ergs, rounded to 2 decimal places
range_M_ej = np.round(np.arange(0.8, 1.6, 0.2), 2)   # solar mass, rounded to 2 decimal places
range_kappa = np.round([0.03, 0.05], 2)              # cm^2/g, rounded to 2 decimal places
range_t_delay = np.round(np.arange(1e2, 1e4, 200), 0)  # s, rounded to 0 decimal places (integer-like)
range_f_comp = np.round([1.5], 2)                    # compress fraction, rounded to 2 decimal places
range_M_dom = np.round(np.arange(0.08, 0.16, 0.02), 2)  # solar mass, rounded to 2 decimal places
range_v_dom = np.int_(np.round(np.arange(5e3, 10e3, 1e3), 0))  # km/s, rounded and converted to integer
range_f_dom = np.round(np.arange(0.03, 0.16, 0.02), 2)  # fraction of DOM mass, rounded to 2 decimal places

def process_combination(args):
    E_exp, M_ej, kappa, t_delay, f_comp, M_dom, V_dom, f_dom, fit_tbl = args
    header_parameters = ['E_exp','M_ej','kappa','t_delay','f_comp','M_dom','V_dom','f_dom']
    header_fitvalues = ['exptime_DEI', 'exptime_FB']
    header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
    for filter_ in fit_filterset:
        header_fitvalues.append(f'alpha_{filter_}')
        header_fitvalues.append(f'amplitude_{filter_}')
    tot_header = header_parameters + header_fitvalues + header_fitconfig
    result_tbl = Table(names = tot_header)
    result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))
    try:
        print('Start calculation: ', E_exp, M_ej, kappa, t_delay, f_comp, M_dom, V_dom, f_dom)
        result = fit_both(
                          fit_tbl=fit_tbl,
                          E_exp=E_exp,
                          M_ej=M_ej,
                          kappa=kappa,
                          M_dom=M_dom,
                          V_dom=V_dom,
                          f_dom=f_dom,
                          t_delay=t_delay,
                          f_comp=f_comp,
                          fit_method='leastsq'
                          )
        data_parameters = dict(E_exp=E_exp, M_ej=M_ej, kappa=kappa, t_delay=t_delay, f_comp=f_comp, M_dom=M_dom, V_dom=V_dom, f_dom=f_dom)
        data_fitvalues = result.params.valuesdict()
        data_fitconfig = dict(success=result.success, nfev=result.nfev, ndata=result.ndata, nvar=result.nvarys, chisq=result.chisqr, redchisqr=result.redchi, aic=result.aic, bic=result.bic)
        all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
        all_values = [all_data[colname] for colname in result_tbl.columns]
        result_tbl.add_row(vals=all_values)
    except:
        data_parameters = dict(E_exp=E_exp, M_ej=M_ej, kappa=kappa, t_delay=t_delay, f_comp=f_comp, M_dom=M_dom, V_dom=V_dom, f_dom=f_dom)
        data_fitvalues = {value: 99999 for value in header_fitvalues}
        data_fitconfig = dict(success=False, nfev=99999, ndata=99999, nvar=99999, chisq=99999, redchisqr=99999, aic=99999, bic=99999)
        all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
        all_values = [all_data[colname] for colname in result_tbl.columns]
        result_tbl.add_row(vals=all_values)
    os.makedirs(f'./result/DOM_fit_result/kappa{kappa}/E{E_exp}')
    result_tbl.write(f'./result/DOM_fit_result/kappa{kappa}/E{E_exp}/{E_exp}_{M_ej}_{kappa}_{t_delay}_{f_comp}_{M_dom}_{V_dom}_{f_dom}.fit', format='ascii.fixed_width', overwrite=True)

def main(fit_tbl):
    start = time.time()
    os.makedirs('./result/DOM_fit_result', exist_ok=True)
    
    # Prepare the list of all combinations of parameters
    all_combinations = [(E_exp, M_ej, kappa, t_delay, f_comp, M_dom, V_dom, f_dom, fit_tbl)
                        for E_exp in range_E_exp
                        for M_ej in range_M_ej
                        for kappa in range_kappa
                        for t_delay in range_t_delay
                        for f_comp in range_f_comp
                        for M_dom in range_M_dom
                        for V_dom in range_v_dom
                        for f_dom in range_f_dom]

    # Use multiprocessing to process the combinations in parallel
    with mp.Pool(processes=8) as pool:
        pool.map(process_combination, all_combinations)

    result_tbl.remove_row(index=0)
    print(time.time() - start)

if __name__ == '__main__':
    main(fit_tbl = fit_tbl)



#%%
result_2 = result_tbl
#%%
result_tbl.write('./DOM_fire_Result_3.txt', format ='ascii.fixed_width')
#%%
result_tbl.sort('chisq')
#%%
result_tbl = ascii.read('./DOM_fire_Result_3.txt', format = 'fixed_width')
fit_filterset = 'UBVgri'
result_tbl.sort('chisq')
#%%
import matplotlib
matplotlib.use('TkAgg')  # Or 'Agg', 'Qt5Agg', etc. depending on your system
import matplotlib.pyplot as plt

result_tbl.sort('chisq')
i = 1
result_values = result_tbl[i]
exptime_DEI = result_values['exptime_DEI']
exptime_FB = result_values['exptime_FB']
filter_key = fit_tbl.group_by('filter').groups.keys['filter']
color_key, offset_key, _, _, label_key = helper.load_filt_keys()
plt.figure(dpi = 300, figsize = (5, 8))
plt.gca().invert_yaxis()
phase_min_FB = np.max([59526, result_values['exptime_FB']])
phase_min_DEI = np.max([59526, result_values['exptime_DEI']])
phase_range_FB = np.arange(phase_min_FB, 59540, 0.1)
phase_range_DEI = np.arange(phase_min_DEI, 59540, 0.1)

DOM_model = DOMInteractionL17(E_exp = result_values['E_exp'], M_ej = result_values['M_ej'], kappa = result_values['kappa'], t_delay = result_values['t_delay'], f_comp = result_values['f_comp'], M_dom = result_values['M_dom'], V_dom = result_values['V_dom'], f_dom = result_values['f_dom'])
DOM_LC = DOM_model.get_LC(td = np.arange(0.1, 10, 0.1), filterset  = ''.join(filter_key), search_directory = model_directory, save = True, force_calculate= True)
spl_allfilt_DEI = get_DEI_spline(DOM_LC, exptime_DEI = result_values['exptime_DEI'], filterset = ''.join(filter_key))

tbl_UL = observed_data.get_data_ul()
tbl_obs = observed_data.get_data_detected()
for filter_ in fit_filterset:
    #exptime_FB = out.params[f'exptime_{filter_}']
    amp = result_values[f'amplitude_{filter_}']
    alpha= result_values[f'alpha_{filter_}']
    flux_FB = fireball_model(time = phase_range_DEI, amplitude = amp, alpha = alpha, exptime = result_values['exptime_FB'])
    spl_DEI = spl_allfilt_DEI[filter_]
    flux_DEI = helper.mag_to_flux(spl_DEI(phase_range_DEI)+DM)
    flux_both = flux_FB+ flux_DEI
    mag_model = helper.flux_to_mag(flux_FB, zp = ZP)
    mag_DOM = helper.flux_to_mag(flux_DEI, zp = ZP)
    mag_both = helper.flux_to_mag(flux_both, zp = ZP)
    tbl_UL_filter = tbl_UL[tbl_UL['filter'] == filter_]
    tbl_obs_filter = tbl_obs[tbl_obs['filter'] == filter_]
    observed_data.show_lightcurve(day_binsize = 5, scatter_linewidth=0.5, scatter_size=40, errorbar_linewidth=0.5, errorbar_capsize=3, color_BV = False, color_gr = False, UL = True, label = True, label_location=0, color_UB = False)
    #plt.axvline(exptime_FB, linestyle = '--', c='r')
    #plt.axvline(exptime_DEI, linestyle = '--', c='b')
    #plt.text(exptime_FB+0.1,20,'FB_EXPTIME = %.4f'%exptime_FB, c='r')
    #plt.text(exptime_DEI+0.1,21,'DEI_EXPTIME = %.4f'%exptime_DEI, c='b')
    #plt.scatter(tbl_obs_filter['obsdate'], tbl_obs_filter['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_])
    #plt.scatter(tbl_UL_filter['obsdate'], tbl_UL_filter['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_], alpha = 0.2)
    plt.plot(phase_range_DEI, mag_model + offset_key[filter_], c = color_key[filter_], label = rf'[{label_key[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
    plt.plot(phase_range_DEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= ':', linewidth = 1)
    plt.plot(phase_range_DEI, mag_both + offset_key[filter_], c = color_key[filter_], linestyle= '-', linewidth = 1)
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
#plt.legend(loc = 4)
#observed_data.show_lightcurve( day_binsize = 5, color_BV = False, color_gr = False, color_UB = False, UL = True, label = False, label_location=2, scatter_size= 120)
#observed_data.show_lightcurve( day_binsize = 5, color_BV = True, color_gr = True, color_UB = True, UL = True, label = False, label_location=2, scatter_size= 120)
'''
plt.plot(phase_range_DEI, mag_U_model - mag_B_model, c = 'cyan', label = 'U-B', linestyle= '--', linewidth = 1)
plt.plot(phase_range_DEI, mag_B_model - mag_V_model, c = 'b', label = 'B-V', linestyle= '--', linewidth = 1)
plt.plot(phase_range_DEI, mag_g_model - mag_r_model, c = 'g', label = 'g-r', linestyle= '--', linewidth = 1)
plt.plot(phase_range_DEI, mag_U_CEI - mag_B_CEI, c = 'cyan', label = 'U-B', linestyle= ':', linewidth = 1)
plt.plot(phase_range_DEI, mag_B_CEI - mag_V_CEI, c = 'b', label = 'B-V', linestyle= ':', linewidth = 1)
plt.plot(phase_range_DEI, mag_g_CEI - mag_r_CEI, c = 'g', label = 'g-r', linestyle= ':', linewidth = 1)
'''
#plt.plot(phase_range_DEI, mag_U_both - mag_B_both, c = 'cyan', label = 'U-B', linestyle= '-', linewidth = 1)
#plt.plot(phase_range_DEI, mag_B_both - mag_V_both, c = 'b', label = 'B-V', linestyle= '-', linewidth = 1)
#plt.plot(phase_range_DEI, mag_g_both - mag_r_both, c = 'g', label = 'g-r', linestyle= '-', linewidth = 1)
#plt.ylim(-1, 1.7)

plt.xticks(np.min(obs_tbl['obsdate']) + np.arange(-20, 200, 5), np.arange(-20, 200, 5) )

# U band
'''
filter_ = 'U'
tbl_filter = observed_data.get_filt_data(obs_tbl)[filter_]
tbl_obs = tbl_filter[tbl_filter['status'] =='detected']
tbl_UL = tbl_filter[tbl_filter['status'] =='UL']
plt.scatter(tbl_obs['obsdate'], tbl_obs['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_])
plt.scatter(tbl_UL['obsdate'], tbl_UL['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_], alpha = 0.2)


spl_DEI = spl_allfilt_DEI[filter_]
flux_DEI = mag_to_flux(spl_DEI(phase_range_DEI)+DM)
mag_DOM = flux_to_mag(flux_DEI, zp = ZP)
plt.plot(phase_range_DEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= ':', linewidth = 1)
'''

plt.xlim(phase_range_FB[0]-3, 59540)
plt.ylim(22.5, 8)
plt.show()
#plt.legend(loc = 4)
# %%

#%%
pandas_tbl = result_tbl.to_pandas()
import seaborn as sns
drop_colnames = ['exptime_DEI', 'exptime_FB', 'alpha_B', 'amplitude_B',
       'alpha_V', 'amplitude_V', 'alpha_g', 'amplitude_g', 'alpha_r',
       'amplitude_r', 'alpha_i', 'amplitude_i', 'success', 'nfev', 'ndata',
       'nvar', 'aic', 'bic']
pd_tbl = pandas_tbl.drop(columns = drop_colnames)
sns.heatmap(pd_tbl.corr())
# %%
plt.figure(figsize=(8, 8))
# Store heatmap object in a variable to easily access it when you want to include more features (such as title).
# Set the range of values to be displayed on the colormap from -1 to 1, and set the annotation to True to display the correlation values on the heatmap.
heatmap = sns.heatmap(pd_tbl.corr(), vmin=-1, vmax=1, annot=True)
# Give a title to the heatmap. Pad defines the distance of the title from the top of the heatmap.
heatmap.set_title('Correlation Heatmap', fontdict={'fontsize':10}, pad=12);
# %%
