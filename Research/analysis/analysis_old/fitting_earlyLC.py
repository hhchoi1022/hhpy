
#%%
from HHsupport_analysis import read_Polin2019_spec
from scipy.ndimage import median_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from HHsupport_analysis import load_filt_keys
from HHsupport_analysis import read_Polin2019
from HHsupport_analysis import mag_to_flux
from HHsupport_analysis import interpolate_spline
from astropy.table import Table
from astropy.table import vstack
from lmfit import Parameters, minimize
from HHsupport_analysis import read_Polin2019
from astropy.table import vstack
from observedphot import ObservedPhot
from HHsupport_analysis import synphot_Polin2019_spec
from companion_interaction_K10 import CompanionInteractionK10
distancemodulus = 31.17
#%% Observation
# noextin

filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_noextin.dat'
filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_noextin.dat'
filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_noextin.dat'

# extin2.36

filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat'
filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostmwextin3.10.dat'
filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostmwextin3.10.dat'
filepath_all = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Alldata_No_cor.dat'

tbl1 = ascii.read(filepath_1, format = 'fixed_width')
tbl2 = ascii.read(filepath_2, format = 'fixed_width')
tbl3 = ascii.read(filepath_3, format = 'fixed_width')

tbl_obs = vstack([tbl1, tbl2])
tbl_obs = ascii.read(filepath_all, format ='fixed_width')
observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
observed_data.data.sort('obsdate')
#observed_data.data = observed_data.get_data_detected()
err_tbl_key = observed_data.data['e_mag'] < 0.01
observed_data.data['e_mag'][err_tbl_key] = 0.03
observed_data.exclude_observatory(['LasCumbres0.4m', 'Swift', 'Swope'])
#observed_data.exclude_filter('u')
plt.figure(dpi = 500, figsize = (7.2,4))
plt.gca().invert_yaxis()
plt.grid()
observed_data.show_lightcurve(label = True, color_BV = True, color_gr = True, color_UB = True)
obs_tbl = observed_data.data
plt.xlim(59525, 59540)

#%%
'''from HHsupport_analysis import match_table
from HHsupport_analysis import binning_table
#observed_data.exclude_observatory('LasCumbres1m')
obs_tbl = observed_data.get_data_detected()
filt_data = observed_data.get_filt_data(obs_tbl)
color_tbl_BV = match_table(filt_data['B'] , filt_data['V'], 'obsdate', 0.05)
color_tbl_BV = binning_table(color_tbl_BV, key = 'obsdate_1', tolerance= 0.1)
color_tbl_gr = match_table(filt_data['g'] , filt_data['r'], 'obsdate', 0.05)
plt.figure(figsize = (5.6,5))
plt.scatter(color_tbl_BV['obsdate_1'], color_tbl_BV['mag_1']-color_tbl_BV['mag_2'] + 0.1059, c = 'g', marker ='*', s = 80, label = '2021aefx')
#plt.scatter(color_tbl _gr['obsdate_1'], color_tbl_gr['mag_1']-color_tbl_gr['mag_2'])
plt.xlim(59529.3315-0.25, 59529.3315+10.15)
plt.ylim(-0.35, 0.73)
plt.savefig('myfigure.png', transparent = True)'''
# %% Load model Polin 2019

def load_polin_spec(wd_mass : float,
                    he_shell_mass : float):
    model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/spectrum/Polin/ddet_Polin2019/good_{wd_mass}_{shell}_doubledet.h5')
    Lnu, Llamb, lamb, time, mu = read_Polin2019_spec(model_file_key(wd_mass, he_shell_mass))
    return Lnu, Llamb, lamb, time, mu

def get_polin_spec_day(wd_mass : float,
                       he_shell_mass : float,
                       td : float,
                       fnu : bool = None):
    Lnu, Llamb, lamb, time, mu = load_polin_spec(wd_mass, he_shell_mass)
    spec_idx = np.argmin(np.abs(time - td))
    phase = time[spec_idx]
    flux = Llamb
    if fnu:
        flux = Lnu
    spec = flux[spec_idx]
    wl = lamb
    return spec, wl, phase

def load_polin_synphot(wd_mass : float,
                    he_shell_mass : float):
    model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/spectrum/Polin/ddet_Polin2019/good_{wd_mass}_{shell}_doubledet.h5')
    model_spec = synphot_Polin2019_spec(model_file_key(wd_mass, he_shell_mass))
    return model_spec


def load_polin_phot(wd_mass : float,
                    he_shell_mass : float,
                    edge_lit : bool = True):
    model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/lightcurve/Polin/ddet_Polin2019/{shell}HeShell/good_{wd_mass}_{shell}_doubledet_ugrizUBVRI.mag')
    if edge_lit:
        model_file_key = (lambda wd_mass, shell: f'/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/lightcurve/Polin/ddet_Polin2019/{shell}HeShell/good_{wd_mass}_{shell}_edgelit_ugrizUBVRI.mag')
    model_phot = read_Polin2019(model_file_key(wd_mass, he_shell_mass))
    return model_phot
#%%
#%%
def get_comp_interaction(rstar : float,
                         wdmass : float,
                         exptime : float,
                         filterset :str,
                         phase_max : int,
                         smooth : float = 0.001):
    trange = np.arange(0.01, phase_max, 0.2)
    Comp = CompanionInteractionK10(rstar = rstar, m_wd = wdmass, commonangle = False)
    Comp_tbl, _, _ = Comp.calc_magnitude(trange, filterset = filterset)
    spl_dict = dict()
    for filter_ in filterset:
        mag_kasen = Comp_tbl[filter_]
        spl_kasen_filt, _ = interpolate_spline(trange + exptime, mag_kasen, show = False, smooth = smooth)
        spl_dict[filter_] = spl_kasen_filt
    return spl_dict
#%%
def get_ddet(exptime : float,
             model_tbl,
             filterset : str,
             phase_max : int,
             smooth : float = 0.05
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

def fit_both(wdmass,
             he_shell_mass,
             fit_filterset = 'BVgr',
             fit_start_mjd : int = 59529,
             fit_end_mjd : int = 59539,
             fit_method = 'leastsq',
             synphot : bool = False):
    model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
    if synphot:
        model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)

    phase_max = fit_end_mjd - fit_start_mjd

    fit_idx = [filter_ in fit_filterset for filter_ in obs_tbl['filter']]
    fit_tbl = obs_tbl[fit_idx]
    early_fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
    early_fit_tbl.sort('obsdate')
    early_fit_tbl['flux'] = mag_to_flux(early_fit_tbl['mag'])
    early_fit_tbl['e_flux'] = early_fit_tbl['e_mag']*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    filter_tbls = early_fit_tbl.group_by('filter').groups
    filter_key = filter_tbls.keys['filter']
    fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
    x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
    y_fit = [np.array((fit_table[filter_]['flux'].tolist())) for filter_ in filter_key]
    e_y_fit = [np.array((fit_table[filter_]['e_flux'].tolist())) for filter_ in filter_key]

    fit_params_DD_CEI = Parameters()
    fit_params_DD_CEI.add('magoffset', value = 0.3, min = 0, max = 1, brute_step = 0.1)
    fit_params_DD_CEI.add('exptime_polin', value = 59528.5, min = 59525, max = 59535, brute_step = 0.5)
    fit_params_DD_CEI.add('exptime_kasen', value = 59528.5, min = 59525, max = 59535, brute_step = 0.5)
    fit_params_DD_CEI.add('rstar', value = 5, min = 1, max = 15, brute_step = 1)

    def chisq_both(params, x_fit, y_fit, y_err, filter_key, 
                    model_polin,
                    phase_max,
                    wdmass):
        vals = params.valuesdict()
        exptime_polin = vals['exptime_polin']
        exptime_kasen = vals['exptime_kasen']
        magoffset = vals['magoffset']
        rstar = vals['rstar']
        chisq_allfilter = []
        print('Exptime:%.3f,Rstar:%.2f,Offset:%.3f'%(exptime_polin, rstar, magoffset))
        kasen_spl_allfilt = get_comp_interaction(rstar = rstar, wdmass = wdmass, exptime = exptime_kasen, filterset = filter_key, phase_max = phase_max)
        polin_spl_allfilt = get_ddet(exptime = exptime_polin , model_tbl = model_polin, filterset = filter_key, phase_max = phase_max)
        for mjd, obs_flux, obs_fluxerr, filter_ in zip(x_fit, y_fit, y_err, filter_key):
            kasen_spl = kasen_spl_allfilt[filter_]
            polin_spl = polin_spl_allfilt[filter_]
            model_flux = mag_to_flux(polin_spl(mjd)+distancemodulus + magoffset ) + mag_to_flux(kasen_spl(mjd)+distancemodulus + magoffset)
            chisq_singlefilter = (((obs_flux - model_flux)/obs_fluxerr)**2)
            chisq_allfilter.append(chisq_singlefilter)
        print(np.sum(np.concatenate(chisq_allfilter)))
        return np.concatenate(chisq_allfilter)
    
    out = minimize(chisq_both, fit_params_DD_CEI, args = (x_fit, y_fit, e_y_fit, filter_key, model_polin, phase_max, wdmass), method = fit_method)
    return out


#%% Observation data
# Filterting fit table
wdmasslist = [0.7, 0.8, 0.85, 0.9, 0.95, 0.96, 1.0, 1.1, 1.2]
he_shell_masslist = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
#wdmasslist = [1.1]
#he_shell_masslist = [0.01]

synphot = False
fit_start_mjd = 59529
fit_end_mjd = 59539
exptime_kasenlist = []
exptime_polinlist = []
rstarlist = []
magoffsetlist = []
chisqlist = []
mass_wdlist = []
he_shelllist = []
for wdmass in wdmasslist:
    for he_shell_mass in he_shell_masslist:
        try:
            out = fit_both(wdmass, he_shell_mass, 'UBVgri', synphot = synphot, fit_start_mjd =fit_start_mjd, fit_end_mjd= fit_end_mjd)
            params = out.params
            exptime_polin = params['exptime_polin'].value
            exptime_kasen = params['exptime_kasen'].value
            rstar = params['rstar'].value
            magoffset = params['magoffset'].value
            chisq = out.redchi
            magoffsetlist.append(magoffset)
            rstarlist.append(rstar)
            exptime_polinlist.append(exptime_polin)
            exptime_kasenlist.append(exptime_kasen)
            chisqlist.append(chisq)
            he_shelllist.append(he_shell_mass)
            mass_wdlist.append(wdmass)
        except:
            pass
result_tbl = Table()
result_tbl['wdmass'] = mass_wdlist
result_tbl['heshellmass'] = he_shelllist
result_tbl['rstar'] = rstarlist
result_tbl['magoffset'] = magoffsetlist
result_tbl['exptime_polin'] = exptime_polinlist
result_tbl['exptime_kasen'] = exptime_kasenlist
result_tbl['chisq'] = chisqlist
result_tbl.sort('chisq')
#%%
plt.figure(dpi =300, figsize = (5, 5))
i = 0
show_filters = 'UBVgri'
rstar = result_tbl[i]['rstar']
wdmass = result_tbl[i]['wdmass']
magoffset = result_tbl[i]['magoffset']
he_shell_mass = result_tbl[i]['heshellmass']
exptime_polin = result_tbl[i]['exptime_polin']
exptime_kasen = result_tbl[i]['exptime_kasen']
model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
if not synphot:
    model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)

time_range_polin = np.arange(exptime_kasen + 0.1, 59570, 0.01)
time_range_kasen = np.arange(exptime_kasen + 0.1, 59570, 0.01)
color, offset_key, _, _, label_key = load_filt_keys()
scatterlinewidth = 1
mksize = 50
smooth = 0.01

observed_data.show_lightcurve(color_BV = False, color_gr = False, label = False, label_location= 4, day_binsize= 10)
from HHsupport_analysis import flux_to_mag
kasen_spl_allfilt = get_comp_interaction(rstar = rstar, wdmass = wdmass, exptime = exptime_kasen, filterset = show_filters, phase_max = 40, smooth = smooth)
polin_spl_allfilt = get_ddet(exptime = exptime_polin, model_tbl = model_polin, filterset = show_filters, phase_max = 40, smooth  = smooth)
#plt.text(exptime_polin+4, 20.3, s = r'$R_{companion}$ = %.2f $R_{\odot}$'%(round(rstar,2)))
#plt.text(exptime_polin+4, 21.0, s = r'$t_{ddet} = %.4f$' %exptime_polin)
#plt.text(exptime_polin+4, 21.7, s = r'$t_{comp}$  = %.4f' %exptime_kasen)
plt.yticks(np.arange(0, 26, 2),np.arange(0, 26, 2))
#plt.xticks(exptime_kasen + np.arange(0, 12, 2),['0', '2', '4', '6', '8', '10'])
for filter_ in show_filters:
    kasen_spl = kasen_spl_allfilt[filter_]
    polin_spl = polin_spl_allfilt[filter_]
    both = flux_to_mag(mag_to_flux(polin_spl(time_range_polin)+distancemodulus + magoffset) + mag_to_flux(kasen_spl(time_range_polin)+distancemodulus + magoffset))
    plt.plot(time_range_polin, both + offset_key[filter_], c = color[filter_], label = label_key[filter_])
    plt.plot(time_range_kasen, kasen_spl(time_range_kasen) + distancemodulus + offset_key[filter_] + magoffset, c = color[filter_], linestyle = ':')
    plt.plot(time_range_polin, polin_spl(time_range_polin) + distancemodulus+ offset_key[filter_] + magoffset, c = color[filter_], linestyle = '--')
plt.plot([0],[0], linestyle= ':', c='k', label = 'Companion interaction')
plt.plot([0],[0], linestyle= '--', c='k', label = 'Double detonation')
plt.plot([0],[0], linestyle= '-', c='k', label = 'Combined')
#plt.legend()
#plt.title(r'Double detonation [$M_{wd}$ = %.2f, $M_{He}$ = %.2f]'%(wdmass, he_shell_mass), fontsize = 10)
plt.grid(axis = 'y')
plt.xlim(exptime_kasen-1.5, exptime_kasen+40.2)
plt.ylim(25, 7)

# %% Color
filt_data = observed_data.get_data_detected()
obs_filter_data = observed_data.get_filt_data(filt_data)
obs_U = obs_filter_data['U']
obs_B = obs_filter_data['B']
obs_V = obs_filter_data['V']
obs_g = obs_filter_data['g']
obs_r = obs_filter_data['r']
time_U = obs_U['obsdate']
time_B = obs_B['obsdate']
time_V = obs_V['obsdate']
time_g = obs_g['obsdate']
time_r = obs_r['obsdate']

obs_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], show = True, k = 2, smooth = 1)
obs_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'], show = True, k = 2, smooth = 0.3)
obs_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'], show = True, k = 2, smooth = 0.3)
obs_g,_ = interpolate_spline(obs_g['obsdate'], obs_g['mag'], show = True, k = 2, smooth = 1)
obs_r,_ = interpolate_spline(obs_r['obsdate'], obs_r['mag'], show = True, k = 2, smooth = 1)
obs_UB = obs_U(time_U) - obs_B(time_U)
obs_BV = obs_B(time_B) - obs_V(time_B)
obs_gr = obs_g(time_g) - obs_r(time_g)

polin_U = polin_spl_allfilt['U'](time_range_polin)
polin_B = polin_spl_allfilt['B'](time_range_polin)
polin_V = polin_spl_allfilt['V'](time_range_polin)
polin_g = polin_spl_allfilt['g'](time_range_polin)
polin_r = polin_spl_allfilt['r'](time_range_polin)
polin_UB = polin_U - polin_B
polin_BV = polin_B - polin_V
polin_gr = polin_g - polin_r

kasen_U = kasen_spl_allfilt['U'](time_range_polin)
kasen_B = kasen_spl_allfilt['B'](time_range_polin)
kasen_V = kasen_spl_allfilt['V'](time_range_polin)
kasen_g = kasen_spl_allfilt['g'](time_range_polin)
kasen_r = kasen_spl_allfilt['r'](time_range_polin)
kasen_UB = kasen_U - kasen_B
kasen_BV = kasen_B - kasen_V
kasen_gr = kasen_g - kasen_r

both_U = flux_to_mag(mag_to_flux(polin_U) + mag_to_flux(kasen_U))
both_B = flux_to_mag(mag_to_flux(polin_B) + mag_to_flux(kasen_B))
both_V = flux_to_mag(mag_to_flux(polin_V) + mag_to_flux(kasen_V))
both_g = flux_to_mag(mag_to_flux(polin_g) + mag_to_flux(kasen_g))
both_r = flux_to_mag(mag_to_flux(polin_r) + mag_to_flux(kasen_r))
both_UB = both_U - both_B
both_BV = both_B - both_V
both_gr = both_g - both_r
# %% B-V
plt.figure(dpi = 300, figsize = (5, 2))
plt.axhline(0, c= 'k')

plt.plot(time_range_polin-exptime_kasen, both_UB, c= 'b')
plt.scatter(time_U-exptime_kasen, obs_UB, label = 'U-B', marker = 'o', edgecolors='b', facecolors =  'None')

plt.plot(time_range_polin-exptime_kasen, both_BV, c= 'g')
plt.scatter(time_B-exptime_kasen, obs_BV, label = 'B-V', marker = 'o', edgecolors='g', facecolors =  'None')

plt.plot(time_range_polin-exptime_kasen, both_gr, c='r')
plt.scatter(time_g-exptime_kasen, obs_gr, label = 'g-r', marker = 'o', edgecolors='r', facecolors =  'None')

plt.xlim(-1.5, +10.2)
plt.ylim(-0.8, 1)
plt.ylabel('color')
plt.xlabel('days since explosion')
plt.legend()


#%%TMP Kasen 2010

plt.figure(dpi =300, figsize = (5, 5))
i = 0
show_filters = 'UBVgri'
wdmass = result_tbl[i]['wdmass']
magoffset = result_tbl[i]['magoffset']
he_shell_mass = result_tbl[i]['heshellmass']
exptime_kasen = result_tbl[i]['exptime_kasen']
plt.figure(dpi = 300, figsize = (5, 4))
rstarlist = np.array([2,5,10,100])#np.array([3e11, 5e11, 2e12, 2e13])/2/7e10
from matplotlib import cm
cmap = cm.get_cmap('viridis')
for i, rstar in enumerate(rstarlist):

    time_range_kasen = np.arange(exptime_kasen + 0.01, exptime_kasen+5, 0.01)
    color, offset_key, _, _, label_key = load_filt_keys()
    scatterlinewidth = 1
    mksize = 50
    smooth = 0.01

    from HHsupport_analysis import flux_to_mag
    kasen_spl_allfilt = get_comp_interaction(rstar = rstar, wdmass = wdmass, exptime = exptime_kasen, filterset = show_filters, phase_max = 40, smooth = smooth)
    for filter_ in show_filters:
        kasen_spl = kasen_spl_allfilt[filter_]

    filt_data = observed_data.get_data_detected()
    obs_filter_data = observed_data.get_filt_data(filt_data)
    obs_U = obs_filter_data['U']
    obs_B = obs_filter_data['B']
    obs_V = obs_filter_data['V']
    obs_g = obs_filter_data['g']
    obs_r = obs_filter_data['r']
    time_U = obs_U['obsdate']
    time_B = obs_B['obsdate']
    time_V = obs_V['obsdate']
    time_g = obs_g['obsdate']
    time_r = obs_r['obsdate']

    obs_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], show = False, k = 2, smooth = 0.3)
    obs_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'], show = False, k = 2, smooth = 0.3)
    obs_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'], show = False, k = 2, smooth = 0.3)
    obs_g,_ = interpolate_spline(obs_g['obsdate'], obs_g['mag'], show = False, k = 2, smooth = 0.5)
    obs_r,_ = interpolate_spline(obs_r['obsdate'], obs_r['mag'], show = False, k = 2, smooth = 0.7)
    obs_UB = obs_U(time_U) - obs_B(time_U)
    obs_BV = obs_B(time_B) - obs_V(time_B)
    obs_gr = obs_g(time_g) - obs_r(time_g)

    kasen_U = kasen_spl_allfilt['U'](time_range_kasen)
    kasen_B = kasen_spl_allfilt['B'](time_range_kasen)
    kasen_V = kasen_spl_allfilt['V'](time_range_kasen)
    kasen_g = kasen_spl_allfilt['g'](time_range_kasen)
    kasen_r = kasen_spl_allfilt['r'](time_range_kasen)
    kasen_UB = kasen_U - kasen_B
    kasen_BV = kasen_B - kasen_V
    kasen_gr = kasen_g - kasen_r



    plt.plot(time_range_kasen-exptime_kasen, kasen_BV, c = cmap(i/len(rstarlist)), label =r'$R_{comp} = $%.1f$R_\odot$'%rstar)
    plt.xlim(0, +10.2)
    plt.ylim(-0.8, 1)
    plt.ylabel('B-V')
    plt.xlabel('days since explosion')
plt.legend()

#%% Polin 2019
plt.figure(dpi =300, figsize = (5, 5))
show_filters = 'UBVgri'
magoffset = result_tbl[i]['magoffset']
he_shell_masslist = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1])
he_shell_masslist = np.array([0.01])
exptime_kasen = result_tbl[i]['exptime_kasen']
plt.figure(dpi = 300, figsize = (5, 4))
rstarlist = np.array([2,5,10,100])#np.array([3e11, 5e11, 2e12, 2e13])/2/7e10
from matplotlib import cm
cmap = cm.get_cmap('hot')
for i, he_shell_mass in enumerate(he_shell_masslist):
    try:
        model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
        time_range_kasen = np.arange(exptime_polin + 0.5, 59570, 0.01)
        color, offset_key, _, _, label_key = load_filt_keys()
        scatterlinewidth = 1
        mksize = 50
        smooth = 0.01
        kasen_spl_allfilt = get_comp_interaction(rstar = rstar, wdmass = wdmass, exptime = exptime_kasen, filterset = show_filters, phase_max = 40, smooth = smooth)
        polin_spl_allfilt = get_ddet(exptime = exptime_polin, model_tbl = model_polin, filterset = show_filters, phase_max = 40, smooth  = smooth)

        from HHsupport_analysis import flux_to_mag
        for filter_ in show_filters:
            polin_spl = polin_spl_allfilt[filter_]

        filt_data = observed_data.get_data_detected()
        obs_filter_data = observed_data.get_filt_data(filt_data)
        obs_U = obs_filter_data['U']
        obs_B = obs_filter_data['B']
        obs_V = obs_filter_data['V']
        obs_g = obs_filter_data['g']
        obs_r = obs_filter_data['r']
        time_U = obs_U['obsdate']
        time_B = obs_B['obsdate']
        time_V = obs_V['obsdate']
        time_g = obs_g['obsdate']
        time_r = obs_r['obsdate']

        obs_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], show = False, k = 2, smooth = 0.3)
        obs_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'], show = False, k = 2, smooth = 0.3)
        obs_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'], show = False, k = 2, smooth = 0.3)
        obs_g,_ = interpolate_spline(obs_g['obsdate'], obs_g['mag'], show = False, k = 2, smooth = 0.5)
        obs_r,_ = interpolate_spline(obs_r['obsdate'], obs_r['mag'], show = False, k = 2, smooth = 0.7)
        obs_UB = obs_U(time_U) - obs_B(time_U)
        obs_BV = obs_B(time_B) - obs_V(time_B)
        obs_gr = obs_g(time_g) - obs_r(time_g)

        kasen_U = kasen_spl_allfilt['U'](time_range_kasen)
        kasen_B = kasen_spl_allfilt['B'](time_range_kasen)
        kasen_V = kasen_spl_allfilt['V'](time_range_kasen)
        kasen_g = kasen_spl_allfilt['g'](time_range_kasen)
        kasen_r = kasen_spl_allfilt['r'](time_range_kasen)
        kasen_UB = kasen_U - kasen_B
        kasen_BV = kasen_B - kasen_V
        kasen_gr = kasen_g - kasen_r
        
        polin_U = polin_spl_allfilt['U'](time_range_kasen)
        polin_B = polin_spl_allfilt['B'](time_range_kasen)
        polin_V = polin_spl_allfilt['V'](time_range_kasen)
        polin_g = polin_spl_allfilt['g'](time_range_kasen)
        polin_r = polin_spl_allfilt['r'](time_range_kasen)
        polin_UB = polin_U - polin_B
        polin_BV = polin_B - polin_V
        polin_gr = polin_g - polin_r
        
        both_B = flux_to_mag(mag_to_flux(polin_B) + mag_to_flux(kasen_B))
        both_V = flux_to_mag(mag_to_flux(polin_V) + mag_to_flux(kasen_V))
        both_BV = both_B - both_V


        plt.plot(time_range_kasen-exptime_polin, polin_BV, c = cmap(i/len(he_shell_masslist)), label =r'[Poiin]$M_{He}$ = %.2f$M_\odot$'%( he_shell_mass))
        plt.plot(time_range_kasen-exptime_polin, kasen_BV, c = cmap(i/len(he_shell_masslist)), label =r'[Kasen]$R_{*}$ = %.2f$M=R_\odot$'%( rstar))
        plt.plot(time_range_kasen-exptime_polin, both_BV, c = cmap(i/len(he_shell_masslist)), label =r'[Both]')
        #plt.scatter(time_B - exptime_polin, obs_BV)
        plt.xlim(0, +10.2)
        plt.ylim(-0.8, 3)
        plt.ylabel('B-V')
        plt.xlabel('days since explosion')
    except:
        pass
plt.legend()
#%%

def fit_ddet(wdmass,
             he_shell_mass,
             fit_filterset = 'UBgVri',
             fit_start_mjd : int = 59529,
             fit_end_mjd : int = 59538,
             fit_method = 'leastsq',
             synphot : bool = False,
             edge_lit : bool = False):
    if synphot:
        model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
    else:
        model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass, edge_lit= edge_lit)
    phase_max = fit_end_mjd - fit_start_mjd

    fit_idx = [filter_ in fit_filterset for filter_ in obs_tbl['filter']]
    fit_tbl = obs_tbl[fit_idx]
    early_fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
    early_fit_tbl.sort('obsdate')
    early_fit_tbl['flux'] = mag_to_flux(early_fit_tbl['mag'])
    early_fit_tbl['e_flux'] = early_fit_tbl['e_mag']*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    early_fit_tbl['e_flux'] = 0.05*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    filter_tbls = early_fit_tbl.group_by('filter').groups
    filter_key = filter_tbls.keys['filter']
    fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
    x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
    y_fit = [np.array((fit_table[filter_]['flux'].tolist())) for filter_ in filter_key]
    e_y_fit = [np.array((fit_table[filter_]['e_flux'].tolist())) for filter_ in filter_key]

    fit_params_DD = Parameters()
    fit_params_DD.add('magoffset', value = -3, min = 0, max = 3, brute_step = 0.01)
    fit_params_DD.add('exptime_polin', value = 59528.5, min = 59525.0, max = 59530, brute_step = 0.01)

    def chisq_ddet(params, x_fit, y_fit, y_err, filter_key, 
                model_polin,
                phase_max):
        vals = params.valuesdict()
        exptime_polin = vals['exptime_polin']
        magoffset = vals['magoffset']
        chisq_allfilter = []
        print('Exptime:%.3f,Offset:%.3f'%(exptime_polin, magoffset))
        polin_spl_allfilt = get_ddet(exptime = exptime_polin , model_tbl = model_polin, filterset = filter_key, phase_max = phase_max)
        for mjd, obs_flux, obs_fluxerr, filter_ in zip(x_fit, y_fit, y_err, filter_key):
            polin_spl = polin_spl_allfilt[filter_]
            model_flux = mag_to_flux(polin_spl(mjd)+distancemodulus + magoffset)
            chisq_singlefilter = (((obs_flux - model_flux)/obs_fluxerr)**2)
            chisq_allfilter.append(chisq_singlefilter)
        print(np.sum(np.concatenate(chisq_allfilter)))
        return np.concatenate(chisq_allfilter)

    out = minimize(chisq_ddet, fit_params_DD, args = (x_fit, y_fit, e_y_fit, filter_key, model_polin, phase_max), method = fit_method)
    return out
#%% Observation data
# Filterting fit table
wdmasslist = [0.95, 0.96, 1.0, 1.1, 1.2, 1.3]
he_shell_masslist = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
#wdmasslist = [ 1.1]
#he_shell_masslist = [0.05]
synphot = False
edge_lit = False
fit_start_mjd = 59529
fit_end_mjd = 59538
exptime_kasenlist = []
exptime_polinlist = []
rstarlist = []
magoffsetlist = []
chisqlist = []
mass_wdlist = []
he_shelllist = []
for wdmass in wdmasslist:
    for he_shell_mass in he_shell_masslist:
        try:
            out = fit_ddet(wdmass, he_shell_mass, 'BVgri', synphot = synphot, fit_start_mjd= fit_start_mjd, fit_end_mjd = fit_end_mjd, edge_lit = edge_lit)
            params = out.params
            exptime_polin = params['exptime_polin'].value
            magoffset = params['magoffset'].value
            chisq = out.redchi
            magoffsetlist.append(magoffset)
            exptime_polinlist.append(exptime_polin)
            chisqlist.append(chisq)
            he_shelllist.append(he_shell_mass)
            mass_wdlist.append(wdmass)
        except:
            pass
#%%
result_tbl_ddet = Table()
result_tbl_ddet['wdmass'] = mass_wdlist
result_tbl_ddet['heshellmass'] = he_shelllist
result_tbl_ddet['magoffset'] = magoffsetlist
result_tbl_ddet['exptime_polin'] = exptime_polinlist
result_tbl_ddet['chisq'] = chisqlist
result_tbl_ddet.sort('chisq')
#%%

plt.figure(dpi =300, figsize = (5, 5))
i = 0
show_filters = 'UBVgri'
wdmass = result_tbl_ddet[i]['wdmass']
magoffset = result_tbl_ddet[i]['magoffset']
he_shell_mass = result_tbl_ddet[i]['heshellmass']
exptime_polin = result_tbl_ddet[i]['exptime_polin']

#wdmass = 1.1
#magoffset = 0.5336917553977
#he_shell_mass = 0.05
#exptime_polin = 59529.05954157

if synphot:
    model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
else:
    model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
#model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass, edge_lit= True)

time_range_polin = np.arange(exptime_polin + 0.1, 59570, 0.01)
color, offset_key, _, _, label_key = load_filt_keys()
scatterlinewidth = 1
mksize = 50
smooth = 0.01

observed_data.show_lightcurve(color_BV = False, color_gr = False, label = False, label_location= 4)
from HHsupport_analysis import flux_to_mag
polin_spl_allfilt = get_ddet(exptime = exptime_polin, model_tbl = model_polin, filterset = show_filters, phase_max = 40, smooth  = smooth)
#plt.text(exptime_polin+4, 21.0, s = r'$t_{ddet} = %.4f$' %exptime_polin)
plt.yticks(np.arange(0, 26, 2),np.arange(0, 26, 2))
plt.xticks(exptime_polin + np.arange(0, 12, 2),['0', '2', '4', '6', '8', '10'])
for filter_ in show_filters:
    polin_spl = polin_spl_allfilt[filter_]
    plt.plot(time_range_polin, polin_spl(time_range_polin) + distancemodulus+ offset_key[filter_] + magoffset, c = color[filter_], linestyle = '--')
plt.plot([0],[0], linestyle= ':', c='k', label = 'Companion interaction')
plt.plot([0],[0], linestyle= '--', c='k', label = 'Double detonation')
plt.plot([0],[0], linestyle= '-', c='k', label = 'Combined')
#plt.legend()
#plt.title(r'Double detonation [$M_{wd}$ = %.2f, $M_{He}$ = %.2f]'%(wdmass, he_shell_mass), fontsize = 10)
plt.grid(axis = 'y')
plt.xlim(exptime_polin-1.5, exptime_polin+10.2)
plt.ylim(25, 7)
#%%
plt.figure(dpi = 300)
plt.plot([0],[0], linestyle= ':', c='k', label = r'$R_{comp} = 3.3R_\odot$')
plt.plot([0],[0], linestyle= '--', c='k', label = 'Double detonation')
plt.plot([0],[0], linestyle= '-', c='k', label = 'Combined')
plt.legend()

#%%

show_filters = 'UBVgri'
i = 4
params = out.params
exptime_polin = 59528.8#result_tbl_ddet[i]['exptime_polin']
magoffset = result_tbl_ddet[i]['magoffset']
wdmass = result_tbl_ddet[i]['wdmass']
he_shell_mass = result_tbl_ddet[i]['heshellmass']

if synphot:
    model_polin = load_polin_synphot(wd_mass = wdmass, he_shell_mass= he_shell_mass)
else:
    model_polin = load_polin_phot(wd_mass = wdmass, he_shell_mass= he_shell_mass)    
time_range_polin = np.arange(exptime_polin + 0.1, 59570, 0.01)
color, offset_key, _, _, label_key = load_filt_keys()
scatterlinewidth = 1
mksize = 50
observed_data.show_lightcurve()
from HHsupport_analysis import flux_to_mag
polin_spl_allfilt = get_ddet(exptime = exptime_polin, model_tbl = model_polin, filterset = show_filters, phase_max = 40, smooth = 0.03)
plt.title('WDmass = %.2f, Heshellmass = %.2f'%(wdmass, he_shell_mass))

for filter_ in show_filters:
    polin_spl = polin_spl_allfilt[filter_]
    plt.plot(time_range_polin, polin_spl(time_range_polin) + distancemodulus+ offset_key[filter_] + magoffset, c = color[filter_], linestyle = '--')

plt.xlim(59525, 59550)
plt.ylim(25, 5)

# %% Color

# %% Color
filt_data = observed_data.get_data_detected()
obs_filter_data = observed_data.get_filt_data(filt_data)
obs_U = obs_filter_data['U']
obs_B = obs_filter_data['B']
obs_V = obs_filter_data['V']
obs_g = obs_filter_data['g']
obs_r = obs_filter_data['r']
time_U = obs_U['obsdate']
time_B = obs_B['obsdate']
time_V = obs_V['obsdate']
time_g = obs_g['obsdate']
time_r = obs_r['obsdate']

obs_U,_ = interpolate_spline(obs_U['obsdate'], obs_U['mag'], show = False, k = 2, smooth = 0.3)
obs_B,_ = interpolate_spline(obs_B['obsdate'], obs_B['mag'], show = False, k = 2, smooth = 0.3)
obs_V,_ = interpolate_spline(obs_V['obsdate'], obs_V['mag'], show = False, k = 2, smooth = 0.3)
obs_g,_ = interpolate_spline(obs_g['obsdate'], obs_g['mag'], show = False, k = 2, smooth = 0.5)
obs_r,_ = interpolate_spline(obs_r['obsdate'], obs_r['mag'], show = False, k = 2, smooth = 0.7)
obs_UB = obs_U(time_U) - obs_B(time_U)
obs_BV = obs_B(time_B) - obs_V(time_B)
obs_gr = obs_g(time_g) - obs_r(time_g)

polin_U = polin_spl_allfilt['U'](time_range_polin)
polin_B = polin_spl_allfilt['B'](time_range_polin)
polin_V = polin_spl_allfilt['V'](time_range_polin)
polin_g = polin_spl_allfilt['g'](time_range_polin)
polin_r = polin_spl_allfilt['r'](time_range_polin)
polin_UB = polin_U - polin_B
polin_BV = polin_B - polin_V
polin_gr = polin_g - polin_r

# %% B-V
plt.figure(dpi = 300, figsize = (5, 2))
plt.axhline(0, c= 'k')

plt.plot(time_range_polin-exptime_polin, polin_UB, c= 'b')
plt.scatter(time_U-exptime_polin, obs_UB, label = 'U-B', edgecolors='b', marker = 'o', facecolors = 'None')

plt.plot(time_range_polin-exptime_polin, polin_BV, c= 'g')
plt.scatter(time_B-exptime_polin, obs_BV, label = 'B-V', edgecolors='g', marker = 'o', facecolors =  'None')

plt.plot(time_range_polin-exptime_polin, polin_gr, c='r')
plt.scatter(time_g-exptime_polin, obs_gr, label = 'g-r', marker = 'o', edgecolors='r', facecolors =  'None')
plt.xlim(-1.5, +10.2)
plt.ylim(-0.8, 1.2)
plt.ylabel('color')
plt.xlabel('days since explosion')
plt.legend()






# %%
rstar = 9.10128965629746
wdmass = 1.1
he_mass = 0.01

#rstar = 0
from HHsupport_analysis import mag_to_fnu
from HHsupport_analysis import fnu_to_mag
wl_range = (3000, 9000)
kasen_spec_spl = CompanionInteraction(rstar = rstar, m_wd = wdmass).calc_spectrum(td = 2)
spec, wl, phase = get_polin_spec_day(wd_mass = wdmass, he_shell_mass = he_mass, td = 2, fnu = True)
idx = np.where((wl > wl_range[0]) & (wl < wl_range[1]))
#%%
spec_polin = spec[idx] 
wl_polin = wl[idx]
flux_polin = spec_polin
flux_kasen = mag_to_fnu(kasen_spec_spl(wl_polin))
flux_both = flux_polin + flux_kasen

# %%



#%%
from HHsupport_analysis import nuflux_to_lflux
median_filtersize = 10
f_lamb = False
if f_lamb:
    flux_both = nuflux_to_lflux(flux_both,wl_polin)
    flux_polin = nuflux_to_lflux(flux_polin, wl_polin)
    flux_kasen = nuflux_to_lflux(flux_kasen, wl_polin)
    


plt.figure(figsize = (10,2))
plt.xlim(3500, 9500)
flux_both = median_filter(flux_both , median_filtersize)
flux_polin = median_filter(flux_polin , median_filtersize)
norm_idx = (wl_polin>6995) & (wl_polin< 7005)
norm_factor = np.mean(flux_both[norm_idx])
plt.plot(wl_polin, flux_both, c= 'b')
plt.plot(wl_polin, flux_polin )
plt.plot(wl_polin, flux_kasen )



#plt.plot(wl_polin, flux_polin, c= 'k')

# %%
from wiserepSpectrum import WiseRepSpectrum
plt.figure(figsize = (6, 3))
obs_spec = WiseRepSpectrum('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/spectrum/WISeREP/ascii/*.txt')
#%%
plt.figure(figsize = (6, 8))


for phase in obs_spec.phase[:4]:
    data = obs_spec.get_spec_date(phase)#, normalize= False, show_flux_unit = 'flamb', label = round(phase-59529.859,2))
    plt.step(data.wavelength, data.fnu.value, label =round(phase-59548.470,1))

#plt.ylim(0, 7e-15)
plt.xlim(4000, 6000)
plt.axvline(5875, c='k', linestyle= '--')
plt.axvline(6673, c='k', linestyle= '--')
plt.axvline(6563, c='k', linestyle= '--')
plt.xlabel('Wavelength[AA]')
plt.ylabel(r'$f_{\nu}$ [erg/s/cm2/Hz]')
#plt.grid()
plt.legend()

#plt.step(wl_polin, flux_both/norm_factor)
# %%
from polinspectrum import PolinSpectrum
polin_spec = PolinSpectrum('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/spectrum/Polin/ddet_Polin2019/*.h5')
# %%
obs_spec.show_spec_date(59529, normalize= True, color = 'r')
#obs_spec.show_spec_date(59533, normalize= True)
polin_model = polin_spec.timeSeriesSpectrum(wdmass = 1.2, he_shellmass = 0.01)
color = plt.cm.jet(np.linspace(0,1,10))
for i, time in enumerate(range(1, 4)):
    polin_model.show_spec_date(time, smooth_factor = 11, show_flux_unit='fnu', normalize= True, color=color[i], label = time)
plt.xlim(3000, 9000)
plt.legend()
# %%

plt.xlim(5300, 6500)

cmap = cm.get_cmap('viridis')
show_obs_spec_phase = obs_spec.phase[:6]
for i, phase in enumerate(show_obs_spec_phase):
    obs_spec.show_spec_date(phase, normalize= True, label = f'{round(phase-59548.470,1)}d', normalize_cenwl=6355, show_flux_unit= 'jy', smooth_factor= 5, color = cmap((len(show_obs_spec_phase)-i)/len(show_obs_spec_phase)))
#plt.axvline(5875, c='k', linestyle= '--')
#plt.axvline(6673, c='k', linestyle= '--')
#plt.axvline(6563, c='k', linestyle= '--')
plt.legend()
#%%
from spectrum import Spectrum
for td in np.arange(1, 5, 0.05):
    kasen_spec_spl = CompanionInteraction(rstar = 6, m_wd = 1.2, commonangle= False).calc_spectrum(td = td)
    from HHsupport_analysis import mag_to_fnu
    spec_polin = polin_model.get_spec_date(1)
    flux_polin, wl_polin = spec_polin.flux, spec_polin.wavelength.value
    wl_range = (3000, 10000)
    idx = np.where((wl_polin > wl_range[0]) & (wl_polin < wl_range[1]))
    flux_polin = flux_polin[idx] 
    wl_polin = wl_polin[idx]
    flux_kasen = mag_to_fnu(kasen_spec_spl(wl_polin))
    flux_both = flux_polin.value + flux_kasen
    spec_kasen = Spectrum(wl_polin, flux_kasen, flux_unit= 'fnu')
    spec_polin = Spectrum(wl_polin, flux_polin.value, flux_unit= 'fnu')
    spec_both = Spectrum(wl_polin, flux_both, flux_unit= 'fnu')
    spec_kasen.show(smooth_factor = 5, normalize = False, show_flux_unit='fnu')
    #spec_polin.show(smooth_factor = 11, normalize = True, show_flux_unit= 'flamb')
    #spec_both.show(smooth_factor = 11, color = None, label =td, normalize = True, show_flux_unit= 'flamb')
plt.legend()
#obs_spec.show_spec_date(59529, normalize= True, label = round(phase-59529.859,2), show_flux_unit= 'flamb')

# %%
plt.step(wl_polin, flux_polin)

#plt.step(wl_polin, flux_polin+flux_kasen)
# %%
from wiserepSpectrum import WiseRepSpectrum
# %%

cbv = WiseRepSpectrum('/Users/hhchoi1022/Gitrepo/data/SN2017cbv/observation/spectrum/WISeREP/ascii/*.ascii')
# %%
cbv.show_spec_date(57825, show_flux_unit='flamb', normalize = True, smooth_factor= 41, color ='b', normalize_cenwl= 6400)
obs_spec.show_spec_date(59529, show_flux_unit = 'flamb', normalize = True, smooth_factor= 21, color='r', normalize_cenwl= 6200)
# %%
