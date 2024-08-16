#%%
import numpy as np
from HHsupport_analysis import interpolate_spline
from HHsupport_analysis import load_filt_keys
from astropy.table import Table
from DOM_interaction_L17 import DOMInteractionL17
from astropy.io import ascii
import matplotlib.pyplot as plt
from observedphot import ObservedPhot
from HHsupport_analysis import mag_to_flux, flux_to_mag
from lmfit import Parameters, minimize
from HHsupport_analysis import read_Polin2019, synphot_Polin2019_spec
#%% Observation
DM = 31.14
ZP = 25
filepath_all = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Alldata_H_M_cor310.dat'
tbl_obs = ascii.read(filepath_all, format = 'fixed_width')
tbl_obs['absmag'] = (tbl_obs['mag'] - DM).round(3)
observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
observed_data.exclude_observatory(['Swope'])
plt.figure(dpi = 400, figsize =(4,7))
observed_data.show_lightcurve(day_binsize = 20, scatter_linewidth=0.5, scatter_size=40, errorbar_linewidth=0.5, errorbar_capsize=3, color_BV = True, color_gr = True, UL = True, label = True, label_location=0, color_UB = True)
plt.xlim(59525, 59540)
obs_tbl = observed_data.get_data_detected()
#%%

def load_polin_synphot(wd_mass : float,
                    he_shell_mass : float):
    model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/spectrum/Polin/ddet_Polin2019/good_{wd_mass}_{shell}_doubledet.h5')
    model_spec = synphot_Polin2019_spec(model_file_key(wd_mass, he_shell_mass))
    return model_spec


def load_polin_phot(wd_mass : float,
                    he_shell_mass : float,
                    edge_lit : bool = False):
    model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/lightcurve/Polin/ddet_Polin2019/{shell}HeShell/good_{wd_mass}_{shell}_doubledet_ugrizUBVRI.mag')
    if edge_lit:
        model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/lightcurve/Polin/ddet_Polin2019/{shell}HeShell/good_{wd_mass}_{shell}_edgelit_ugrizUBVRI.mag')
    model_phot = read_Polin2019(model_file_key(wd_mass, he_shell_mass))
    return model_phot

def get_DDet_spline(exptime : float,
                    model_tbl,
                    phase_max : int = 15,
                    filterset : str = 'UBVRIugri',
                    smooth : float = 0.1
                    ):
    #time_range = np.arange(phase_range[0], phase_range[1], 0.05)
    spl_dict = dict()
    for filter_ in filterset:
        model_mag = model_tbl[filter_]
        inf_idx = np.isinf(model_mag)
        time_idx = (model_tbl['phase'][~inf_idx] < phase_max)
        mag_polin = model_mag[~inf_idx][time_idx]
        phase_polin = model_tbl['phase'][~inf_idx][time_idx]
        spl_polin, _ = interpolate_spline(phase_polin + exptime, mag_polin, show = False, smooth = smooth)
        spl_dict[filter_] = spl_polin
    return spl_dict

def fit_ddet(wdmass,
             he_shell_mass,
             fit_filterset = 'UBgVri',
             fit_start_mjd : int = 59529,
             fit_end_mjd : int = 59538,
             fit_method = 'brute',
             synphot : bool = False,
             edge_lit : bool = True):
    if synphot:
        model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
    else:
        model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass, edge_lit= edge_lit)
    phase_max = fit_end_mjd - fit_start_mjd

    fit_idx = [filter_ in fit_filterset for filter_ in obs_tbl['filter']]
    fit_tbl = obs_tbl[fit_idx]
    early_fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
    early_fit_tbl.sort('obsdate')
    #early_fit_tbl['flux'] = mag_to_flux(early_fit_tbl['mag'])
    #early_fit_tbl['e_flux'] = early_fit_tbl['e_mag']*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    #early_fit_tbl['e_flux'] = 0.05*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    filter_tbls = early_fit_tbl.group_by('filter').groups
    filter_key = filter_tbls.keys['filter']
    fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
    x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
    y_fit = [np.array((fit_table[filter_]['mag'].tolist())) for filter_ in filter_key]
    e_y_fit = [np.array((fit_table[filter_]['e_mag'].tolist())) for filter_ in filter_key]

    fit_params_DD = Parameters()
    fit_params_DD.add('magoffset', value = 0, min = 0.5, max = 3, brute_step= 0.01)
    fit_params_DD.add('exptime_polin', value = 59528.5, min = 59525, max = 59530, brute_step= 0.01)

    def chisq_ddet(params, x_fit, y_fit, e_y_fit, filter_key, 
                model_polin,
                phase_max):
        vals = params.valuesdict()
        exptime_polin = vals['exptime_polin']
        magoffset = vals['magoffset']
        chisq_allfilter = []
        #print('Exptime:%.3f,Offset:%.3f'%(exptime_polin, magoffset))
        polin_spl_allfilt = get_DDet_spline(exptime = exptime_polin , model_tbl = model_polin, filterset = filter_key, phase_max = phase_max)
        for mjd, obs_mag, obs_magerr, filter_ in zip(x_fit, y_fit, e_y_fit, filter_key):
            polin_spl = polin_spl_allfilt[filter_]
            model_mag = polin_spl(mjd) + magoffset + DM
            chisq_singlefilter = (np.abs((obs_mag - model_mag)/obs_magerr))
            chisq_allfilter.append(chisq_singlefilter)
        #print(np.sum(np.concatenate(chisq_allfilter)))
        return np.concatenate(chisq_allfilter)

    out = minimize(chisq_ddet, fit_params_DD, args = (x_fit, y_fit, e_y_fit, filter_key, model_polin, phase_max), method = fit_method)
    return out

#%% Fitting parameters
range_wdmass = np.arange(0.6, 1.5, 0.1)
range_he_shell_mass = np.arange(0.01, 0.11, 0.01)
range_edge_lit  = [True, False]
fit_filterset = 'BVgri'
fit_timerange = (59529, 59540)
synphot = False
#%% Fitting
header_parameters =  ['wdmass','he_shell_mass','edge_lit']
header_fitvalues = ['exptime_polin', 'magoffset']
header_fitconfig = ['success','nfev', 'ndata', 'nvar', 'chisq', 'redchisqr', 'aic', 'bic']
tot_header = header_parameters + header_fitvalues + header_fitconfig
result_tbl = Table(names = tot_header)
result_tbl.add_row(vals = np.zeros(len(result_tbl.colnames)))
for wdmass in range_wdmass:
    for he_shell_mass in range_he_shell_mass:
        print(he_shell_mass)
        result_tbl.write('./DDet_Result.txt', format = 'ascii.fixed_width', overwrite = True)
        for edge_lit in range_edge_lit:
            try:
                wdmass = np.round(wdmass,1)
                he_shell_mass = np.round(he_shell_mass,2)
                result = fit_ddet(wdmass = wdmass,
                                he_shell_mass = he_shell_mass,
                                fit_filterset = fit_filterset,
                                fit_start_mjd = fit_timerange[0],
                                fit_end_mjd = fit_timerange[1],
                                fit_method = 'brute',
                                synphot = synphot,
                                edge_lit = edge_lit
                                )
                data_parameters = dict(wdmass = wdmass, he_shell_mass = he_shell_mass, edge_lit = edge_lit)
                data_fitvalues = result.params.valuesdict()
                data_fitconfig = dict(success = True, nfev = result.nfev, ndata = result.ndata, nvar = result.nvarys, chisq = result.chisqr, redchisqr = result.redchi, aic = result.aic, bic = result.bic)
                all_data = {**data_parameters, **data_fitvalues, **data_fitconfig}
                all_values = []
                for colname in result_tbl.columns:
                    value = all_data[colname]
                    all_values.append(value)
                result_tbl.add_row(vals = all_values)
            except:
                pass
result_tbl.remove_row(index = 0)
result_tbl.sort('chisq')
result_tbl.write('./DDet_Result.txt', format = 'ascii.fixed_width', overwrite = True)
#%%
result_tbl = ascii.read('./result_fit/DDet_Result.txt', format = 'fixed_width')
result_tbl.sort('chisq')
#%% Visualization
i =3
result_values = result_tbl[i]
exptime_DDet = result_values['exptime_polin']
magoffset = result_values['magoffset']
wdmass = result_values['wdmass']
he_shell_mass = result_values['he_shell_mass']
edge_lit = result_values['edge_lit']
time_range = np.arange(fit_timerange[0], fit_timerange[1], 0.01)
color_key, offset_key, _, _, label_key = load_filt_keys()
#wdmass = 1.1
#magoffset = 0.5336917553977
#he_shell_mass = 0.05
#exptime_DDet = 59529.05954157


plt.figure(dpi = 300, figsize = (5, 8))
plt.gca().invert_yaxis()

if synphot:
    tbl_DDet = load_polin_synphot(wd_mass= wdmass, he_shell_mass= he_shell_mass)
else:
    tbl_DDet = load_polin_phot(wd_mass=wdmass, he_shell_mass = he_shell_mass, edge_lit = edge_lit)
spl_allfilt_DDet = get_DDet_spline(exptime = exptime_DDet, model_tbl = tbl_DDet, phase_max = fit_timerange[1]-fit_timerange[0])
for filter_ in 'UBVgri':
    tbl_filter = observed_data.get_filt_data(observed_data.data)[filter_]
    tbl_filter.sort('obsdate')
    tbl_UL = tbl_filter[tbl_filter['status'] =='UL']
    tbl_obs = tbl_filter[tbl_filter['status'] =='detected']
    spl = spl_allfilt_DDet[filter_]
    mag = spl(time_range) + magoffset + DM
    #plt.axvline(exptime_FB, linestyle = '--', c='r')
    #plt.axvline(exptime_DEI, linestyle = '--', c='b')
    #plt.text(exptime_FB+0.1,20,'FB_EXPTIME = %.4f'%exptime_FB, c='r')
    #plt.text(exptime_DEI+0.1,21,'DEI_EXPTIME = %.4f'%exptime_DEI, c='b')
    #plt.scatter(tbl_obs['obsdate'], tbl_obs['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_])
    #plt.scatter(tbl_UL['obsdate'], tbl_UL['mag']+offset_key[filter_], facecolor = 'none', edgecolor = color_key[filter_], alpha = 0.2)
    plt.plot(time_range, mag + offset_key[filter_], c = color_key[filter_], linestyle= '-', linewidth = 1)
    if filter_ == 'U':
        mag_U_model = mag
    if filter_ == 'B':
        mag_B_model = mag
    if filter_ == 'V':
        mag_V_model = mag
    if filter_ == 'g':
        mag_g_model = mag
    if filter_ == 'r':
        mag_r_model = mag
#observed_data.show_lightcurve( day_binsize = 5, color_BV = False, color_gr = False, color_UB = False, UL = True, label = False, label_location=2, scatter_size= 120)

observed_data.show_lightcurve( day_binsize = 5, color_BV = True, color_gr = True, color_UB = True, UL = True, label = False, label_location=2, scatter_size= 120)
plt.plot(time_range, mag_U_model - mag_B_model, c = 'cyan', label = 'U-B', linestyle= '-', linewidth = 1)
plt.plot(time_range, mag_B_model - mag_V_model, c = 'b', label = 'B-V', linestyle= '-', linewidth = 1)
plt.plot(time_range, mag_g_model - mag_r_model, c = 'g', label = 'g-r', linestyle= '-', linewidth = 1)

plt.xlim(fit_timerange[0]-3, 59540)
plt.ylim(22.5, 8)
plt.ylim(-1, 1.7)
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
