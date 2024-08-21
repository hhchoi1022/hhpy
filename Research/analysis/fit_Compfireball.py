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
    start = time.time()
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

    result_tbl.remove_row(index=0)
    print(time.time() - start)

if __name__ == '__main__':
    main(fit_tbl = fit_tbl)



#%%
result_tbl.write('./CEI_fire_Result_1.txt', format = 'ascii.fixed_width')
print(time.time() - start)
#%%
result_tbl = ascii.read('./CEI_fire_Result_1.txt', format = 'fixed_width')
result_tbl.sort('chisq')
#%%
fit_filterset = 'UBgVri'
i = 0
result_values = result_tbl[i]
exptime_CEI = result_values['exptime_FB']

exptime_FB = result_values['exptime_CEI']
color_key, offset_key, _, _, label_key = load_filt_keys()
plt.figure(dpi = 300, figsize = (5, 8))
#plt.gca().invert_yaxis()

phase_min_FB = np.max([59526, result_values['exptime_FB']])
phase_min_CEI = np.max([59526, result_values['exptime_CEI']])
phase_range_FB = np.arange(phase_min_FB, 59540, 0.1)
phase_range_CEI = np.arange(phase_min_CEI, 59540, 0.1)
phase_range_CEI = np.arange(np.min([phase_min_FB, phase_min_CEI]), 59540, 0.1)


CEI_spl = simulate_CEI_model(rstar = result_values['rstar'], m_wd = result_values['m_wd'], v9 = result_values['v9'])
spl_allfilt_CEI = get_CEI_spline(CEI_spl, exptime_CEI = result_values['exptime_CEI'])

for filter_ in fit_filterset:
    #exptime_FB = out.params[f'exptime_{filter_}']
    amp = result_values[f'amplitude_{filter_}']
    alpha= result_values[f'alpha_{filter_}']
    tbl_filter = observed_data.get_filt_data(observed_data.data)[filter_]
    tbl_filter.sort('obsdate')
    tbl_UL = tbl_filter[tbl_filter['status'] =='UL']
    tbl_obs = tbl_filter[tbl_filter['status'] =='detected']
    flux_FB = fireball_model(time = phase_range_CEI, amplitude = amp, alpha = alpha, exptime = exptime_FB)
    spl_CEI = spl_allfilt_CEI[filter_]
    flux_CEI = mag_to_flux(spl_CEI(phase_range_CEI)+DM)
    flux_both = flux_FB+ flux_CEI
    mag_model = flux_to_mag(flux_FB, zp = ZP)
    mag_DOM = flux_to_mag(flux_CEI, zp = ZP)
    mag_both = flux_to_mag(flux_both, zp = ZP)
    plt.plot(phase_range_CEI, mag_model + offset_key[filter_], c = color_key[filter_], label = rf'[{label_key[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
    plt.plot(phase_range_CEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= ':', linewidth = 1)
    plt.plot(phase_range_CEI, mag_both + offset_key[filter_], c = color_key[filter_], linestyle= '-', linewidth = 1)
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
observed_data.show_lightcurve( day_binsize = 5, color_BV = False, color_gr = False, color_UB = False, UL = True, label = False, label_location=2, scatter_size= 40)
#observed_data.show_lightcurve( day_binsize = 5, color_BV = True, color_gr = True, color_UB = True, UL = True, label = False, label_location=2, scatter_size= 120)

'''
plt.plot(phase_range_CEI, mag_U_model - mag_B_model, c = 'cyan', label = 'U-B', linestyle= '--', linewidth = 1)
plt.plot(phase_range_CEI, mag_B_model - mag_V_model, c = 'b', label = 'B-V', linestyle= '--', linewidth = 1)
plt.plot(phase_range_CEI, mag_g_model - mag_r_model, c = 'g', label = 'g-r', linestyle= '--', linewidth = 1)
plt.plot(phase_range_CEI, mag_U_CEI - mag_B_CEI, c = 'cyan', label = 'U-B', linestyle= ':', linewidth = 1)
plt.plot(phase_range_CEI, mag_B_CEI - mag_V_CEI, c = 'b', label = 'B-V', linestyle= ':', linewidth = 1)
plt.plot(phase_range_CEI, mag_g_CEI - mag_r_CEI, c = 'g', label = 'g-r', linestyle= ':', linewidth = 1)
'''
#plt.plot(phase_range_CEI, mag_U_both - mag_B_both, c = 'cyan', label = 'U-B', linestyle= '-', linewidth = 1)
#plt.plot(phase_range_CEI, mag_B_both - mag_V_both, c = 'b', label = 'B-V', linestyle= '-', linewidth = 1)
#plt.plot(phase_range_CEI, mag_g_both - mag_r_both, c = 'g', label = 'g-r', linestyle= '-', linewidth = 1)
#plt.xticks(np.min(obs_tbl['obsdate']) + np.arange(-20, 200, 5), np.arange(-20, 200, 5) )

# U band
'''
filter_ = 'U'
tbl_filter = observed_data.get_filt_data(obs_tbl)[filter_]
tbl_obs = tbl_filter[tbl_filter['status'] =='detected']
tbl_UL = tbl_filter[tbl_filter['status'] =='UL']
plt.scatter(tbl_obs['obsdate'], tbl_obs['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_])
plt.scatter(tbl_UL['obsdate'], tbl_UL['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_], alpha = 0.2)


spl_DEI = spl_allfilt_CEI[filter_]
flux_DEI = mag_to_flux(spl_DEI(phase_range_CEI)+DM)
mag_DOM = flux_to_mag(flux_DEI, zp = ZP)
plt.plot(phase_range_CEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= ':', linewidth = 1)
'''
plt.xlim(phase_range_FB[0]-1, 59539)
plt.ylim(22.5, 8)
#plt.ylim(-1, 1.7)
#%%
plt.figure(dpi  =300)
plt.plot(0, 0, '--', c='k', label = 'Power law')
plt.plot(0, 0, ':', c='k', label = 'DEI')
plt.plot(0, 0, '-', c='k', label = 'Power law + DEI')
plt.legend()
# %%


#%%
pandas_tbl = result_tbl.to_pandas()
import seaborn as sns
drop_colnames = ['exptime_CEI', 'exptime_FB', 'alpha_B', 'amplitude_B',
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
