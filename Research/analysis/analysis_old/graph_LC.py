#%%
from astropy.io import ascii
import matplotlib.pyplot as plt
from observedphot import ObservedPhot
from HHsupport_analysis import load_filt_keys
import numpy as np
import matplotlib.patches as mpatches
from HHsupport_analysis import match_table






#%% Load raw data
DM = 31.14
ZP = 25
filepath_all = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Alldata_H_M_cor310.dat'
tbl_obs = ascii.read(filepath_all, format = 'fixed_width')
tbl_obs['absmag'] = (tbl_obs['mag'] - DM).round(3)
observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
observed_data.data  = observed_data.data[observed_data.data['obsdate'] > 59526]
observed_data.exclude_observatory(['Swope','KMTNet_Ni2023'])
data_ul = observed_data.get_data_ul()
data_detected = observed_data.get_data_detected()
# %% Configuration
filters = set(observed_data.data['filter'])
name_telescopes = dict(KCT='KCT', RASA36 = 'RASA36',LSGT = 'LSGT',LasCumbres1m = 'LCO [Hosseinzadeh+22]')
markers = dict(KCT = 's', RASA36 = 'd', LSGT = '^', LasCumbres1m = '.')
sizes = dict(KCT = 40, RASA36 = 30, LSGT = 40, LasCumbres1m = 100)
colors, offsets, _, _, labels = load_filt_keys(filters)
linewidth = 1
UL_length = 1
phase_max = 59547.2325


#%%


# %% LC [All]
plt.figure(figsize = (6.5,4.7), dpi = 400)
plt.gca().invert_yaxis()
ax1 = plt.subplot()
for filter_ in colors.keys():
    for observatory in ['KCT','LSGT','RASA36']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_])
            ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            ax1.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
        if len(show_data_ul) > 0:
            for x,y in zip(show_data_ul['obsdate'], show_data_ul['mag']+offsets[filter_]):
                ax1.arrow(x=x,y=y,dx=0,dy=0.3, linewidth = 0.5, head_width = 2, head_length = 0.15, color = colors[filter_], shape = 'full', alpha = 0.5)
                ax1.arrow(x=x-1,y=y,dx=2,dy=0,linewidth = 0.5, head_width = 0, color = colors[filter_], alpha = 0.5)            
    for observatory in ['LasCumbres1m']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_], alpha  = 0.3)
            ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            ax1.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    ax1.set_ylim(23, 8)
    ax1.set_ylabel('Apparent Magnitude[AB]')
    ax1.set_xlabel('Phase [day]')
    
    #ax2 = ax1.twinx()
    #ax2.set_ylim(23-DM, 7-DM)
    #ax2.set_ylabel('Absolute Magnitude[AB]')
    
    plt.xticks(phase_max + np.arange(-40, 200, 20), np.arange(-40, 200, 20))
plt.xlim(59529.3315 - 8, 59710)
# legends
#sorted_keys = {k: v for k, v in sorted(labels.items(), key=lambda item: item[1])}
sorted_keys = {k: v for k, v in labels.items()}
rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
name_rows = [labels[clr] for clr in sorted_keys]
columns = [plt.plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in list(name_telescopes.keys())]
name_columns = list(name_telescopes.values())
#plt.legend(rows+columns, name_rows+name_columns, loc=1, ncol = 1, fontsize = 9)
plt.legend(columns, name_columns, loc=1, ncol = 1, fontsize = 9)
#%%



# %% LC [Early]
plt.figure(figsize = (4,6), dpi = 400)
plt.gca().invert_yaxis()
for filter_ in colors.keys():
    for observatory in ['KCT','LSGT','RASA36']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_])
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            plt.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
        if len(show_data_ul) > 0:
            for x,y in zip(show_data_ul['obsdate'], show_data_ul['mag']+offsets[filter_]):
                plt.arrow(x=x,y=y,dx=0,dy=0.1, linewidth = 0.5, head_width = 0.25, head_length = 0.3, color = colors[filter_], shape = 'full', alpha = 0.5)
                plt.arrow(x=x-0.125,y=y-0.01,dx=0.25,dy=0,linewidth = 0.5, head_width = 0, color = colors[filter_], alpha = 0.5)            
    for observatory in ['LasCumbres1m']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_], alpha  = 0.3)
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            plt.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    plt.xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
plt.xlim(59529.3315 - 2, 59529.3315+10)
plt.ylim(23, 8)

plt.xlabel('Days since first detection')
plt.ylabel('Apparent Magnitude [AB]')
# legends
#sorted_keys = {k: v for k, v in sorted(labels.items(), key=lambda item: item[1])}
sorted_keys = {k: v for k, v in labels.items()}
rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
name_rows = [labels[clr] for clr in sorted_keys]
columns = [plt.plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in list(name_telescopes.keys())]
name_columns = list(name_telescopes.values())
plt.legend(columns, name_columns, loc=4, ncol = 1, fontsize = 9)

#%%



# %% Load Extinction Corrected data
DM = 31.14
ZP = 25
filepath_all = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Alldata_H_M_cor310.dat'
#filepath_all = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/Alldata_no_cor.dat'

tbl_obs = ascii.read(filepath_all, format = 'fixed_width')
tbl_obs['absmag'] = (tbl_obs['mag'] - DM).round(3)
observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
observed_data.data  = observed_data.data[observed_data.data['obsdate'] > 59526]
observed_data.exclude_observatory(['Swope','KMTNet_Ni2023'])
data_ul = observed_data.get_data_ul()
data_detected = observed_data.get_data_detected()
# %% LC [Early with fireball]
from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
plt.figure(figsize = (4,6), dpi = 400)
plt.gca().invert_yaxis()
for filter_ in colors.keys():
    for observatory in ['KCT','LSGT','RASA36']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_])
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            plt.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
        if len(show_data_ul) > 0:
            for x,y in zip(show_data_ul['obsdate'], show_data_ul['mag']+offsets[filter_]-DM):
                plt.arrow(x=x,y=y,dx=0,dy=0.1, linewidth = 0.5, head_width = 0.25, head_length = 0.3, color = colors[filter_], shape = 'full', alpha = 0.5)
                plt.arrow(x=x-0.125,y=y-0.01,dx=0.25,dy=0,linewidth = 0.5, head_width = 0, color = colors[filter_], alpha = 0.5)            
    for observatory in ['LasCumbres1m']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_], alpha  = 0.3)
            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            plt.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    plt.xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
amplitudelist = dict(U = 35.7756453, B = 101.575783, g = 189.019474, V = 333.903886, r = 319.233812, i = 174.071520)
alphalist = dict(U = 3.14959759, B = 3.01872108, g = 2.74397461, V = 2.38187287, r = 2.42068873, i = 2.4828142)
for filter_ in colors.keys():
    exptime = 59527.3002
    phase_range = np.arange(exptime, 59540, 0.1)
    amp = amplitudelist[filter_]
    alpha= alphalist[filter_]
    flux_model = fireball_model(phase_range, amp, exptime, alpha)
    mag_model = flux_to_mag(flux_model, zp = 25)-DM
    plt.plot(phase_range, mag_model + offsets[filter_], c = colors[filter_], label = rf'[{labels[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
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
        
plt.xlim(59529.3315 - 2, 59529.3315+10)
plt.ylim(23-DM, 8-DM)
plt.xlabel('Days since first detection')
plt.ylabel('Apparent Magnitude [AB]')
# legends
#sorted_keys = {k: v for k, v in sorted(labels.items(), key=lambda item: item[1])}
sorted_keys = {k: v for k, v in labels.items()}
rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
name_rows = [labels[clr] for clr in sorted_keys]
columns = [plt.plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in list(name_telescopes.keys())]
name_columns = list(name_telescopes.values())
plt.legend(columns, name_columns, loc=4, ncol = 1, fontsize = 9)


#%%



# %% Color [Early with fireball]
from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
plt.figure(figsize = (5,2), dpi = 400)
plt.gca().invert_yaxis()

for observatory in ['KCT','LSGT','RASA36']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        plt.scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        plt.scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 1)
        plt.errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        plt.scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        plt.scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        plt.errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        plt.scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        plt.scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        plt.errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
for observatory in ['LasCumbres1m']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        plt.scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        plt.scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 1)
        plt.errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        plt.scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        plt.scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        plt.errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        plt.scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        plt.scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        plt.errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    
plt.xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
amplitudelist = dict(U = 35.7756453, B = 101.575783, g = 189.019474, V = 333.903886, r = 319.233812, i = 174.071520)
alphalist = dict(U = 3.14959759, B = 3.01872108, g = 2.74397461, V = 2.38187287, r = 2.42068873, i = 2.4828142)

for filter_ in colors.keys():
    exptime = 59527.3002
    phase_range = np.arange(exptime, 59540, 0.1)
    amp = amplitudelist[filter_]
    alpha= alphalist[filter_]
    flux_model = fireball_model(phase_range, amp, exptime, alpha)
    mag_model = flux_to_mag(flux_model, zp = 25)
    #plt.plot(phase_range, mag_model + offsets[filter_], c = colors[filter_], label = rf'[{labels[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
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
plt.plot(phase_range, mag_U-mag_B, c = 'cyan', linestyle = '--', linewidth = 1)
plt.plot(phase_range, mag_B-mag_V, c = 'b', linestyle = '--', linewidth = 1)
plt.plot(phase_range, mag_g-mag_r, c = 'k', linestyle = '--', linewidth = 1, label = 'Power law fit')
plt.plot(phase_range, mag_g-mag_r, c = 'r', linestyle = '--', linewidth = 1)

plt.xlim(59529.3315 - 2.5, 59529.3315+10.1)
plt.ylim(-1, 1.5)
plt.xlabel('Days since first detection')
plt.ylabel('Color')
plt.legend()
#%%



# %% LC +Color [Early with fireball]
from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
#plt.figure(figsize = (5,5), dpi = 400)
#plt.gca().invert_yaxis()
fig, axes = plt.subplots(
    nrows=2, ncols=1, 
    sharex=True, # sharing properties among x axes
    figsize=(5, 7),
    dpi = 500,
    gridspec_kw={'height_ratios': [3, 1]})
for filter_ in colors.keys():
    for observatory in ['KCT','LSGT','RASA36']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_])
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            axes[0].errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
        if len(show_data_ul) > 0:
            for x,y in zip(show_data_ul['obsdate'], show_data_ul['mag']+offsets[filter_]-DM):
                axes[0].arrow(x=x,y=y,dx=0,dy=0.1, linewidth = 0.5, head_width = 0.25, head_length = 0.3, color = colors[filter_], shape = 'full', alpha = 0.5)
                axes[0].arrow(x=x-0.125,y=y-0.01,dx=0.25,dy=0,linewidth = 0.5, head_width = 0, color = colors[filter_], alpha = 0.5)            
    for observatory in ['LasCumbres1m']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, label = labels[filter_], alpha  = 0.3)
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            axes[0].errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    axes[0].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
amplitudelist = dict(U = 35.7756453, B = 101.575783, g = 189.019474, V = 333.903886, r = 319.233812, i = 174.071520)
alphalist = dict(U = 3.14959759, B = 3.01872108, g = 2.74397461, V = 2.38187287, r = 2.42068873, i = 2.4828142)
for filter_ in colors.keys():
    exptime = 59527.3002
    phase_range = np.arange(exptime, 59540, 0.1)
    amp = amplitudelist[filter_]
    alpha= alphalist[filter_]
    flux_model = fireball_model(phase_range, amp, exptime, alpha)
    mag_model = flux_to_mag(flux_model, zp = 25-DM)
    axes[0].plot(phase_range, mag_model + offsets[filter_], c = colors[filter_], label = rf'[{labels[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
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
axes[0].set_ylabel('Absolute Magnitude [AB]')
axes[0].set_xlim(59529.3315 - 2, 59529.3315+10.1)
axes[0].set_ylim(23-DM, 8-DM)
axes[0].set_yticks(np.arange(-24, -8, 2), np.arange(-24, -8, 2) )

# legend
legend1_1 = axes[0].plot(0, 0, c = 'k', linestyle = '--', linewidth = 1, label = 'Power-law rise')
legend2_1 = axes[0].scatter(0, 0, marker = markers['LasCumbres1m'], facecolors = 'none',  edgecolors = 'k', s = sizes['LasCumbres1m'], linewidth = linewidth*0.5, label = 'H22')
legend2_2 = axes[0].scatter(0, 0, marker = markers['KCT'], facecolors = 'none',  edgecolors = 'k', s = sizes['KCT'], linewidth = linewidth*0.5, label = 'KCT')
legend2_3 = axes[0].scatter(0, 0, marker = markers['RASA36'], facecolors = 'none',  edgecolors = 'k', s = sizes['RASA36'], linewidth = linewidth*0.5, label = 'RASA36')
legend2_4 = axes[0].scatter(0, 0, marker = markers['LSGT'], facecolors = 'none',  edgecolors = 'k', s = sizes['LSGT'], linewidth = linewidth*0.5, label = 'LSGT')
legend1 = axes[0].legend(handles=[legend1_1[0]], loc=2)
legend2 = axes[0].legend(handles = [legend2_1, legend2_2, legend2_3, legend2_4], loc = 4)
axes[0].add_artist(legend1)

# text
text1 = axes[0].text(59529.3315+ 8.5, -17 - 5.7, f'U{offsets["U"]}', c = colors['U'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 4.3, f'B{offsets["B"]}', c = colors['B'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 3.2, f'g{offsets["g"]}', c = colors['g'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 1.8, f'  V', c = colors['V'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 0.1, f'r+{offsets["r"]}', c = colors['r'])
text1 = axes[0].text(59529.3315+ 8.5, -17 + 2.5, f'i+{offsets["i"]}', c = colors['i'])



from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
for observatory in ['KCT','LSGT','RASA36']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
for observatory in ['LasCumbres1m']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    
axes[1].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
amplitudelist = dict(U = 35.7756453, B = 101.575783, g = 189.019474, V = 333.903886, r = 319.233812, i = 174.071520)
alphalist = dict(U = 3.14959759, B = 3.01872108, g = 2.74397461, V = 2.38187287, r = 2.42068873, i = 2.4828142)

for filter_ in colors.keys():
    exptime = 59527.3002
    phase_range = np.arange(exptime, 59540, 0.1)
    amp = amplitudelist[filter_]
    alpha= alphalist[filter_]
    flux_model = fireball_model(phase_range, amp, exptime, alpha)
    mag_model = flux_to_mag(flux_model, zp = 25)
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
axes[1].plot(phase_range, mag_U-mag_B, c = 'cyan', linestyle = '--', linewidth = 1)
axes[1].plot(phase_range, mag_B-mag_V, c = 'b', linestyle = '--', linewidth = 1)
axes[1].plot(phase_range, mag_g-mag_r, c = 'r', linestyle = '--', linewidth = 1)

axes[1].set_xlim(59529.3315 - 2.5, 59529.3315+10.1)
axes[1].set_ylim(-1, 1.5)
axes[1].set_xlabel('Days since first detection')
axes[1].set_ylabel('Color')

    
axes[1].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )

#label
axes[1].set_xlim(59529.3315 - 2.5, 59529.3315+10.1)
axes[1].set_ylim(-1, 1.5)
axes[1].set_xlabel('Days since first detection')
axes[1].set_yticks(np.arange(-1, 1.6, 0.5), np.arange(-1, 1.6, 0.5))
axes[1].set_ylabel('Color')

# Text
text1 = axes[1].text(59529.3315+ 8.5, 0.8, f'U-B', c = colors['U'])
text1 = axes[1].text(59529.3315+ 8.5, -0.7, f'B-V', c = 'b')
text1 = axes[1].text(59529.3315+ 8.5, 0.1, f'g-r', c = 'r')

plt.subplots_adjust(
    wspace=0, # the width of the padding between subplots 
    hspace=0) # the height of the padding between subplots

#%%



#%% LC +Color [Early with fireball+CEI]
#%% Helper
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
def fireball_filter(params, time, filter_):
    exptime = params['exptime_fireball']
    amp = params[f'amp_{filter_}']
    alpha= params[f'alpha_{filter_}']
    return fireball_model(np.array(time), amp, exptime, alpha)

from companion_interaction_K10 import CompanionInteractionK10
def simulate_CEI_model(rstar : float,
                       wdmass : float,
                       v9 : float,
                       t_range : np.array = np.arange(0.1, 15, 0.1),
                       filterset :str = 'UBVRIugri',
                       model : bool = False
                       ):
    Comp = CompanionInteractionK10(rstar = rstar, m_wd = wdmass, v9 = v9, commonangle = False)
    Comp_tbl, _, _ = Comp.calc_magnitude(t_range, filterset = filterset)
    if model:
        return Comp_tbl, Comp
    return Comp_tbl

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
        spl, _ = interpolate_spline(phase_CEI + exptime_CEI, mag_CEI, show = False, smooth = smooth)
        spl_dict[filter_] = spl
    return spl_dict
def fit_both(rstar,
             wdmass,
             v9 : float = 1,
             
             fit_filterset = 'BVgri',
             fit_start_mjd : int = 59529,
             fit_end_mjd : int = 59538,
             fit_method = 'leastsq',
             simulate : bool = True
             ):

    # Input
    fit_idx = [filter_ in fit_filterset for filter_ in tbl_obs['filter']]
    fit_tbl = tbl_obs[fit_idx]
    early_fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
    early_fit_tbl.sort('obsdate')
    early_fit_tbl['flux'] = mag_to_flux(early_fit_tbl['mag'])
    early_fit_tbl['e_flux'] = 0.05*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5#early_fit_tbl['e_mag']*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    filter_tbls = early_fit_tbl.group_by('filter').groups
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
        fit_params_CEI_FB.add(f'alpha_{filter_}', value = 2, min = 1, max = 4)
        fit_params_CEI_FB.add(f'amplitude_{filter_}', value = 1000, min = 1, max = 500000)
    
    def chisq_both(params, x_fit, y_fit, e_y_fit, filter_key, 
                   model_CEI,
                   ):
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
            CEI_flux = mag_to_flux(spl_CEI(mjd)+DM)
            both_flux = fireball_flux + CEI_flux
            chisq_singlefilter = (((obs_flux - both_flux)/obs_fluxerr)**2)
            #chisq_singlefilter = (np.abs((obs_flux - both_flux)))

            chisq_allfilter.append(chisq_singlefilter)
        #print(np.sum(np.concatenate(chisq_allfilter)))
        #print(vals)
        return np.concatenate(chisq_allfilter)
    
    # Fitting
    phase_max = fit_end_mjd - fit_start_mjd
    t_range = np.arange(0.1, phase_max, 0.1)
    if simulate:
        model_CEI = simulate_CEI_model(t_range = t_range, rstar = rstar, wdmass= wdmass, v9 = v9, filterset = fit_filterset)
    else: 
        model_CEI = simulate_CEI_model(t_range = t_range, rstar = rstar, wdmass= wdmass, v9 = v9, filterset = fit_filterset)
        
    out = minimize(chisq_both, fit_params_CEI_FB, args = (x_fit, y_fit, e_y_fit, filter_key, model_CEI), method = fit_method)
    
    return out
#%% Code
from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
result_tbl = ascii.read('./CEI_fire_Result_1.txt', format = 'fixed_width')
result_tbl.sort('chisq')

#plt.figure(figsize = (5,5), dpi = 400)
#plt.gca().invert_yaxis()
fig, axes = plt.subplots(
    nrows=2, ncols=1, 
    sharex=True, # sharing properties among x axes
    figsize=(5, 8),
    dpi = 500,
    gridspec_kw={'height_ratios': [3, 1]})
#axes[0].set_title('Companion-ejecta + Power law')
for filter_ in colors.keys():
    for observatory in ['KCT','LSGT','RASA36']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5)
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            axes[0].errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
        if len(show_data_ul) > 0:
            for x,y in zip(show_data_ul['obsdate'], show_data_ul['mag']+offsets[filter_]-DM):
                axes[0].arrow(x=x,y=y,dx=0,dy=0.1, linewidth = 0.5, head_width = 0.25, head_length = 0.3, color = colors[filter_], shape = 'full', alpha = 0.5)
                axes[0].arrow(x=x-0.125,y=y-0.01,dx=0.25,dy=0,linewidth = 0.5, head_width = 0, color = colors[filter_], alpha = 0.5)            
    for observatory in ['LasCumbres1m']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5, alpha  = 1)
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            axes[0].errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    axes[0].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
axes[0].set_ylabel('Absolute Magnitude [AB] + Offset')
axes[0].set_xlim(59529.3315 - 2, 59529.3315+10.1)
axes[0].set_ylim(23-DM, 8-DM)
axes[0].set_yticks(np.arange(-24, -8, 2), np.arange(-24, -8, 2) )

# legend
legend1_1 = axes[0].plot(0, 0, c = 'k', linestyle = '--', linewidth = 1, label = 'Power-law rise')
legend1_2 = axes[0].plot(0, 0, c = 'k', linestyle = ':', linewidth = 1, label = 'Companion-ejecta interaction')
legend1_3 = axes[0].plot(0, 0, c = 'k', linestyle = '-', linewidth = 1, label = 'Combined (Power + Companion)')
legend2_1 = axes[0].scatter(0, 0, marker = markers['LasCumbres1m'], facecolors = 'none',  edgecolors = 'k', s = sizes['LasCumbres1m'], linewidth = linewidth*0.5, label = 'H22')
legend2_2 = axes[0].scatter(0, 0, marker = markers['KCT'], facecolors = 'none',  edgecolors = 'k', s = sizes['KCT'], linewidth = linewidth*0.5, label = 'KCT')
legend2_3 = axes[0].scatter(0, 0, marker = markers['RASA36'], facecolors = 'none',  edgecolors = 'k', s = sizes['RASA36'], linewidth = linewidth*0.5, label = 'RASA36')
legend2_4 = axes[0].scatter(0, 0, marker = markers['LSGT'], facecolors = 'none',  edgecolors = 'k', s = sizes['LSGT'], linewidth = linewidth*0.5, label = 'LSGT')
legend1 = axes[0].legend(handles=[legend1_1[0], legend1_2[0], legend1_3[0]], loc=2)
legend2 = axes[0].legend(handles = [legend2_1, legend2_2, legend2_3, legend2_4], loc = 4)
axes[0].add_artist(legend1)

# text
text1 = axes[0].text(59529.3315+ 8.5, -17 - 5.7, f'U{offsets["U"]}', c = colors['U'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 4.3, f'B{offsets["B"]}', c = colors['B'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 3.2, f'g{offsets["g"]}', c = colors['g'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 1.8, f'  V', c = colors['V'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 0.1, f'r+{offsets["r"]}', c = colors['r'])
text1 = axes[0].text(59529.3315+ 8.5, -17 + 2.5, f'i+{offsets["i"]}', c = colors['i'])


result_values = result_tbl[8]
exptime_CEI = result_values['exptime_FB']

exptime_FB = result_values['exptime_CEI']
color_key, offset_key, _, _, label_key = load_filt_keys()
#plt.gca().invert_yaxis()
phase_min_FB = np.max([59526, result_values['exptime_FB']-0.5])
phase_min_CEI = np.max([59526, result_values['exptime_CEI']-0.5])
phase_range_FB = np.arange(phase_min_FB-1, 59540, 0.1)
phase_range_CEI = np.arange(phase_min_CEI-1, 59540, 0.1)
phase_range_CEI = np.arange(np.min([phase_min_FB, phase_min_CEI]), 59540, 0.1)

CEI_spl = simulate_CEI_model(rstar = result_values['rstar'], wdmass = result_values['wdmass'], v9 = result_values['v9'])
spl_allfilt_CEI = get_CEI_spline(CEI_spl, exptime_CEI = result_values['exptime_CEI'])

for filter_ in colors.keys():
    #exptime_FB = out.params[f'exptime_{filter_}']
    amp = result_values[f'amplitude_{filter_}']
    alpha= result_values[f'alpha_{filter_}']
    tbl_filter = observed_data.get_filt_data(observed_data.data)[filter_]
    tbl_filter.sort('obsdate')
    tbl_UL = tbl_filter[tbl_filter['status'] =='UL']
    tbl_obs = tbl_filter[tbl_filter['status'] =='detected']
    flux_FB = fireball_model(time = phase_range_CEI, amplitude = amp, alpha = alpha, exptime = exptime_FB)
    flux_FB = mag_to_flux(flux_to_mag(flux_FB, zp = ZP)-DM)
    spl_CEI = spl_allfilt_CEI[filter_]
    flux_CEI = mag_to_flux(spl_CEI(phase_range_CEI))
    flux_FB[np.isnan(flux_FB)] = 0
    flux_both = flux_FB+ flux_CEI
    mag_model = flux_to_mag(flux_FB, zp = ZP)
    mag_DOM = flux_to_mag(flux_CEI, zp = ZP)
    mag_both = flux_to_mag(flux_both, zp = ZP)
    axes[0].plot(phase_range_CEI, mag_model + offset_key[filter_], c = color_key[filter_], label = rf'[{label_key[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
    axes[0].plot(phase_range_CEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= ':', linewidth = 1)
    axes[0].plot(phase_range_CEI, mag_both + offset_key[filter_], c = color_key[filter_], linestyle= '-', linewidth = 1)
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



# legends
#sorted_keys = {k: v for k, v in sorted(labels.items(), key=lambda item: item[1])}
sorted_keys = {k: v for k, v in labels.items()}
rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
name_rows = [labels[clr] for clr in sorted_keys]
columns = [axes[0].plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in list(name_telescopes.keys())]
name_columns = list(name_telescopes.values())
#axes[0].legend(columns, name_columns, loc=4, ncol = 1, fontsize = 9)

from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
for observatory in ['KCT','LSGT','RASA36']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 0.3)
        axes[1].errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.3)
        axes[1].errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.3)
        axes[1].errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
for observatory in ['LasCumbres1m']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 0.3)
        axes[1].errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.3)
        axes[1].errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.3)
        axes[1].errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    
axes[1].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
axes[1].plot(phase_range_CEI, mag_U_both-mag_B_both, c = 'cyan', linestyle = '-', linewidth = 1)
axes[1].plot(phase_range_CEI, mag_B_both-mag_V_both, c = 'b', linestyle = '-', linewidth = 1)

#label
#axes[1].plot(phase_range_CEI, mag_g_both-mag_r_both, c = 'k', linestyle = '-', linewidth = 1, label = 'Combined model')
axes[1].plot(phase_range_CEI, mag_g_both-mag_r_both, c = 'r', linestyle = '-', linewidth = 1)
axes[1].set_xlim(59529.3315 - 2.5, 59529.3315+10.1)
axes[1].set_ylim(-1, 1.5)
axes[1].set_xlabel('Days since first detection')
axes[1].set_yticks(np.arange(-1, 1.6, 0.5), np.arange(-1, 1.6, 0.5))
axes[1].set_ylabel('Color')

# Text
text1 = axes[1].text(59529.3315+ 8.5, 0.8, f'U-B', c = colors['U'])
text1 = axes[1].text(59529.3315+ 8.5, -0.7, f'B-V', c = 'b')
text1 = axes[1].text(59529.3315+ 8.5, 0.1, f'g-r', c = 'r')


plt.subplots_adjust(
    wspace=0, # the width of the padding between subplots 
    hspace=0) # the height of the padding between subplots


residual = dict()
phase = dict()
error = dict()
for filter_ in colors.keys():
    data_filt = data_detected[(data_detected['filter'] == filter_)] 
    data_filt = data_filt[data_filt['obsdate']<59540]
    if len(data_filt )> 0:
    #exptime_FB = out.params[f'exptime_{filter_}']
        amp = result_values[f'amplitude_{filter_}']
        alpha= result_values[f'alpha_{filter_}']
        flux_FB = fireball_model(time = data_filt['obsdate'], amplitude = amp, alpha = alpha, exptime = exptime_FB)
        flux_FB = mag_to_flux(flux_to_mag(flux_FB, zp = ZP)-DM)
        flux_FB[np.isnan(flux_FB)] = 0
        spl_CEI = spl_allfilt_CEI[filter_]
        flux_CEI = mag_to_flux(spl_CEI(data_filt['obsdate']))
        flux_both = flux_FB+ flux_CEI
        mag_model = flux_to_mag(flux_FB, zp = ZP)
        mag_DOM = flux_to_mag(flux_CEI, zp = ZP)
        mag_both = flux_to_mag(flux_both, zp = ZP)
        residual[filter_] = data_filt['mag']-DM - mag_both
        phase[filter_] = data_filt['obsdate']
        error[filter_] = data_filt['e_mag']
#%% Residual
for filter_ in residual.keys():
    residual_ = residual[filter_]
    phase_ = phase[filter_]
    error_ = error[filter_]
    plt.figure(figsize = (6,2), dpi = 300)
    plt.axhline(0, linestyle = '--', c= 'k')
    plt.scatter(phase_, residual_, c = colors[filter_])
    #plt.errorbar(phase_, residual_, yerr = error_, fmt = 'None', ecolor = 'k')
    plt.ylim(-1, 1)
    plt.xlim(59529.3315-2.5, 59529.3315+10.1)



#%%



#%% LC +Color [Early with fireball+DEI]
#%% Helper
def load_DEI_model(E_exp = 1.0,
                   M_ej = 1.0,
                   kappa = 0.05,
                   
                   M_dom = 0.1,
                   v_dom = 5e3,
                   f_dom = 0.1,
                   
                   t_delay = 5e3,
                   f_comp = 1.5,
                   home_dir : str = '/Users/hhchoi1022/Gitrepo/Research/Supernova/DOM_model/'
                   ):
    
    tbl = ascii.read(f'{home_dir}{E_exp}_{M_ej}_{kappa}_{t_delay}_{f_comp}_{M_dom}_{v_dom}_{f_dom}.dat', format = 'fixed_width')
    return tbl

def simulate_DEI_model(t_range : np.array = np.arange(0.1, 15, 0.1),
                       E_exp = 1.0,
                       M_ej = 1.0,
                       kappa = 0.05,
                   
                       M_dom = 0.1,
                       v_dom = 5e3,
                       f_dom = 0.1,
                       
                       t_delay = 5e3,
                       f_comp = 1.5,
                       ):
    model = DOMInteractionL17(E_exp = E_exp, M_ej = M_ej, kappa = kappa, M_dom = M_dom, V_dom = v_dom, f_dom = f_dom, t_delay = t_delay, f_comp = f_comp)
    tbl, _, _ = model.calc_magnitude(t_range, 'UBVRIugri')
    return tbl

def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    np.nan_to_num(flux, copy=False, nan=0.0001)
    return flux

def fireball_filter(params, time, filter_):
    exptime = params['exptime_fireball']
    amp = params[f'amp_{filter_}']
    alpha= params[f'alpha_{filter_}']
    return fireball_model(np.array(time), amp, exptime, alpha)
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
        spl, _ = interpolate_spline(phase_DEI + exptime_DEI, mag_DEI, show = False, smooth = smooth)
        spl_dict[filter_] = spl
    return spl_dict
def fit_both(E_exp,
             M_ej,
             kappa,
             M_dom,
             v_dom,
             f_dom,
             t_delay,
             f_comp,
             
             fit_filterset = 'BVgr',
             fit_start_mjd : int = 59529,
             fit_end_mjd : int = 59539,
             fit_method = 'leastsq',
             simulate : bool = False
             ):

    # Input
    fit_idx = [filter_ in fit_filterset for filter_ in tbl_obs['filter']]
    fit_tbl = tbl_obs[fit_idx]
    early_fit_tbl = fit_tbl[(fit_tbl['obsdate'] > fit_start_mjd)&(fit_tbl['obsdate'] < fit_end_mjd)]
    early_fit_tbl.sort('obsdate')
    early_fit_tbl['flux'] = mag_to_flux(early_fit_tbl['mag'])
    early_fit_tbl['e_flux'] =  0.05*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5#early_fit_tbl['e_mag']*mag_to_flux(early_fit_tbl['mag'])*2.303/2.5
    filter_tbls = early_fit_tbl.group_by('filter').groups
    filter_key = filter_tbls.keys['filter']
    fit_table = {filter_:filter_tbl for filter_, filter_tbl in zip(filter_key, filter_tbls)}
    x_fit = [np.array((fit_table[filter_]['obsdate'].tolist())) for filter_ in filter_key]
    y_fit = [np.array((fit_table[filter_]['flux'].tolist())) for filter_ in filter_key]
    e_y_fit = [np.array((fit_table[filter_]['e_flux'].tolist())) for filter_ in filter_key]

    # Parameters
    fit_params_DEI_FB = Parameters()
    fit_params_DEI_FB.add('exptime_DEI', value = 59528.5, min = 59528.3, max = 59535)
    fit_params_DEI_FB.add('exptime_FB', value = 59528.5, min = 59528.3, max = 59535)
    for filter_ in filter_key:
        fit_params_DEI_FB.add(f'alpha_{filter_}', value = 2, min = 0, max = 4)
        fit_params_DEI_FB.add(f'amplitude_{filter_}', value = 3000, min = 1, max = 500000)
    
    def chisq_both(params, x_fit, y_fit, e_y_fit, filter_key, 
                   model_DEI,
                   ):
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
            DEI_flux = mag_to_flux(spl_DEI(mjd)+DM)
            both_flux = fireball_flux + DEI_flux
            chisq_singlefilter = (((obs_flux - both_flux)/obs_fluxerr)**2)
            #chisq_singlefilter = (np.abs((obs_flux - both_flux)))

            chisq_allfilter.append(chisq_singlefilter)
        #print(np.sum(np.concatenate(chisq_allfilter)))
        #print(vals)
        return np.concatenate(chisq_allfilter)
    
    # Fitting
    phase_max = fit_end_mjd - fit_start_mjd
    t_range = np.arange(0.1, phase_max, 0.1)
    if simulate:
        model_DEI = simulate_DEI_model(t_range = t_range, E_exp = E_exp, M_ej = M_ej, kappa = kappa, M_dom = M_dom, v_dom = v_dom, f_dom = f_dom, t_delay = t_delay, f_comp = f_comp)
    else: 
        model_DEI = load_DEI_model(E_exp = E_exp, M_ej = M_ej, kappa = kappa, M_dom = M_dom, v_dom = v_dom, f_dom = f_dom, t_delay = t_delay, f_comp = f_comp)
        
    out = minimize(chisq_both, fit_params_DEI_FB, args = (x_fit, y_fit, e_y_fit, filter_key, model_DEI), method = fit_method)
    
    return out

#%% Code

from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
result_tbl = ascii.read('./DOM_fire_Result_3.txt', format = 'fixed_width')
result_tbl.sort('chisq')
fig, axes = plt.subplots(
    nrows=2, ncols=1, 
    sharex=True, # sharing properties among x axes
    figsize=(5, 8),
    dpi = 500,
    gridspec_kw={'height_ratios': [3, 1]})
#axes[0].set_title('DOM-ejecta + Power law')
for filter_ in colors.keys():
    for observatory in ['KCT','LSGT','RASA36']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5,)
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            axes[0].errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
        if len(show_data_ul) > 0:
            for x,y in zip(show_data_ul['obsdate'], show_data_ul['mag']+offsets[filter_]-DM):
                axes[0].arrow(x=x,y=y,dx=0,dy=0.1, linewidth = 0.5, head_width = 0.25, head_length = 0.3, color = colors[filter_], shape = 'full', alpha = 0.5)
                axes[0].arrow(x=x-0.125,y=y-0.01,dx=0.25,dy=0,linewidth = 0.5, head_width = 0, color = colors[filter_], alpha = 0.5)            
    for observatory in ['LasCumbres1m']:
        show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
        show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
        if len(show_data_detected) > 0:
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = sizes[observatory], linewidth = linewidth*0.5,  alpha  = 0.3)
            axes[0].scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 0.2)
            axes[0].errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_]-DM, show_data_detected['e_mag'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    axes[0].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
axes[0].set_ylabel('Absolute Magnitude [AB] + Offset')
axes[0].set_xlim(59529.3315 - 2, 59529.3315+10.1)
axes[0].set_ylim(23-DM, 8-DM)
axes[0].set_yticks(np.arange(-24, -8, 2), np.arange(-24, -8, 2) )

# legend
legend1_1 = axes[0].plot(0, 0, c = 'k', linestyle = '--', linewidth = 1, label = 'Power-law rise')
legend1_2 = axes[0].plot(0, 0, c = 'k', linestyle = ':', linewidth = 1, label = 'DOM-ejecta interaction')
legend1_3 = axes[0].plot(0, 0, c = 'k', linestyle = '-', linewidth = 1, label = 'Combined (Power + DOM)')
legend2_1 = axes[0].scatter(0, 0, marker = markers['LasCumbres1m'], facecolors = 'none',  edgecolors = 'k', s = sizes['LasCumbres1m'], linewidth = linewidth*0.5, label = 'H22')
legend2_2 = axes[0].scatter(0, 0, marker = markers['KCT'], facecolors = 'none',  edgecolors = 'k', s = sizes['KCT'], linewidth = linewidth*0.5, label = 'KCT')
legend2_3 = axes[0].scatter(0, 0, marker = markers['RASA36'], facecolors = 'none',  edgecolors = 'k', s = sizes['RASA36'], linewidth = linewidth*0.5, label = 'RASA36')
legend2_4 = axes[0].scatter(0, 0, marker = markers['LSGT'], facecolors = 'none',  edgecolors = 'k', s = sizes['LSGT'], linewidth = linewidth*0.5, label = 'LSGT')
legend1 = axes[0].legend(handles=[legend1_1[0], legend1_2[0], legend1_3[0]], loc=2)
legend2 = axes[0].legend(handles = [legend2_1, legend2_2, legend2_3, legend2_4], loc = 4)
axes[0].add_artist(legend1)

# text
text1 = axes[0].text(59529.3315+ 8.5, -17 - 5.7, f'U{offsets["U"]}', c = colors['U'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 4.3, f'B{offsets["B"]}', c = colors['B'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 3.2, f'g{offsets["g"]}', c = colors['g'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 1.8, f'  V', c = colors['V'])
text1 = axes[0].text(59529.3315+ 8.5, -17 - 0.1, f'r+{offsets["r"]}', c = colors['r'])
text1 = axes[0].text(59529.3315+ 8.5, -17 + 2.5, f'i+{offsets["i"]}', c = colors['i'])

result_values = result_tbl[0]
exptime_DEI = result_values['exptime_DEI']
exptime_FB = result_values['exptime_FB']
color_key, offset_key, _, _, label_key = load_filt_keys()

phase_min_FB = np.max([59526, result_values['exptime_FB']-0.5])
phase_min_DEI = np.max([59526, result_values['exptime_DEI']-0.5])
phase_range_FB = np.arange(phase_min_FB, 59540, 0.1)
phase_range_DEI = np.arange(phase_min_DEI, 59540, 0.1)

DOM_spl = simulate_DEI_model(E_exp = result_values['E_exp'], M_ej = result_values['M_ej'], kappa = result_values['kappa'], t_delay = result_values['t_delay'], f_comp = result_values['f_comp'], M_dom = result_values['M_dom'], v_dom = result_values['v_dom'], f_dom = result_values['f_dom'])
spl_allfilt_DEI = get_DEI_spline(DOM_spl, exptime_DEI = result_values['exptime_DEI'])


for filter_ in colors.keys():
    #exptime_FB = out.params[f'exptime_{filter_}']
    amp = result_values[f'amplitude_{filter_}']
    alpha= result_values[f'alpha_{filter_}']
    tbl_filter = observed_data.get_filt_data(tbl_obs)[filter_]
    tbl_filter.sort('obsdate')
    tbl_UL = tbl_filter[tbl_filter['status'] =='UL']
    tbl_obs = tbl_filter[tbl_filter['status'] =='detected']
    flux_FB = fireball_model(time = phase_range_DEI, amplitude = amp, alpha = alpha, exptime = exptime_FB)
    flux_FB = mag_to_flux(flux_to_mag(flux_FB, zp = ZP)-DM)
    flux_FB[np.isnan(flux_FB)] = 0
    spl_DEI = spl_allfilt_DEI[filter_]
    flux_DEI = mag_to_flux(spl_DEI(phase_range_DEI))
    flux_both = flux_FB+ flux_DEI
    mag_model = flux_to_mag(flux_FB, zp = ZP)
    mag_DOM = flux_to_mag(flux_DEI, zp = ZP)
    mag_both = flux_to_mag(flux_both, zp = ZP)
    #plt.axvline(exptime_FB, linestyle = '--', c='r')
    #plt.axvline(exptime_DEI, linestyle = '--', c='b')
    #plt.text(exptime_FB+0.1,20,'FB_EXPTIME = %.4f'%exptime_FB, c='r')
    #plt.text(exptime_DEI+0.1,21,'DEI_EXPTIME = %.4f'%exptime_DEI, c='b')
    axes[0].plot(phase_range_DEI, mag_model + offset_key[filter_], c = color_key[filter_], label = rf'[{label_key[filter_]}] $\alpha = {round(alpha,2)}$', linestyle= '--', linewidth = 1)
    axes[0].plot(phase_range_DEI, mag_DOM + offset_key[filter_], c = color_key[filter_], linestyle= ':', linewidth = 1)
    axes[0].plot(phase_range_DEI, mag_both + offset_key[filter_], c = color_key[filter_], linestyle= '-', linewidth = 1)
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



# legends
#sorted_keys = {k: v for k, v in sorted(labels.items(), key=lambda item: item[1])}
sorted_keys = {k: v for k, v in labels.items()}
rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
name_rows = [labels[clr] for clr in sorted_keys]
columns = [axes[0].plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in list(name_telescopes.keys())]
name_columns = list(name_telescopes.values())
#axes[0].legend(columns, name_columns, loc=4, ncol = 1, fontsize = 9)

from HHsupport_analysis import flux_to_mag
def fireball_model(time, amplitude, exptime, alpha):
    flux = amplitude * (time - exptime )**alpha
    return flux
for observatory in ['KCT','LSGT','RASA36']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
for observatory in ['LasCumbres1m']:
    show_data_detected = data_detected[(data_detected['observatory'] == observatory)]
    show_data_ul = data_ul[data_ul['observatory'] == observatory]
    filt_data = observed_data.get_filt_data(show_data_detected)
    if (len(filt_data['U']) > 0) & (len(filt_data['B'])>0):
        UB_tbl = match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 0.15)
        UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'], facecolors = colors['U'],  edgecolors = 'k', s =  sizes[observatory], marker = markers[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2'],  UB_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)
    if (len(filt_data['B']) > 0) & (len(filt_data['V'])>0):
        BV_tbl = match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 0.15)
        BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'], facecolors = 'b',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2'],  BV_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    if (len(filt_data['g']) > 0) & (len(filt_data['r'])>0):
        gr_tbl = match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 0.15)
        gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'none',  edgecolors = 'k', marker = markers[observatory], s = sizes[observatory], linewidth = linewidth*0.5)
        axes[1].scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], facecolors = 'r',  edgecolors = 'k', marker = markers[observatory], s =  sizes[observatory], linewidth = linewidth*0.5, alpha = 1)
        axes[1].errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = linewidth*0.5, capsize = 0, c = 'k', capthick = linewidth*0.5)

    
axes[1].set_xticks(np.min(data_detected['obsdate']) + np.arange(-40, 200, 2), np.arange(-40, 200, 2) )
axes[1].plot(phase_range_DEI, mag_U_both-mag_B_both, c = 'cyan', linestyle = '-', linewidth = 1)
axes[1].plot(phase_range_DEI, mag_B_both-mag_V_both, c = 'b', linestyle = '-', linewidth = 1)
#axes[1].plot(phase_range_DEI, mag_g_both-mag_r_both, c = 'k', linestyle = '--', linewidth = 1, label = 'Combined model')
axes[1].plot(phase_range_DEI, mag_g_both-mag_r_both, c = 'r', linestyle = '-', linewidth = 1)

axes[1].set_xlim(59529.3315 - 2.5, 59529.3315+10.1)
axes[1].set_ylim(-1, 1.5)
axes[1].set_xlabel('Days since first detection')
axes[1].set_yticks(np.arange(-1, 1.6, 0.5), np.arange(-1, 1.6, 0.5))
axes[1].set_ylabel('Color')

# Text
text1 = axes[1].text(59529.3315+ 8.5, 0.8, f'U-B', c = colors['U'])
text1 = axes[1].text(59529.3315+ 8.5, -0.7, f'B-V', c = 'b')
text1 = axes[1].text(59529.3315+ 8.5, 0.1, f'g-r', c = 'r')

#plt.legend(loc = 3)
plt.subplots_adjust(
    wspace=0, # the width of the padding between subplots 
    hspace=0) # the height of the padding between subplots

residual = dict()
phase = dict()
for filter_ in colors.keys():
    data_filt = data_detected[(data_detected['filter'] == filter_)] 
    data_filt = data_filt[data_filt['obsdate']<59540]
    if len(data_filt )> 0:
    #exptime_FB = out.params[f'exptime_{filter_}']
        amp = result_values[f'amplitude_{filter_}']
        alpha= result_values[f'alpha_{filter_}']
        flux_FB = fireball_model(time = data_filt['obsdate'], amplitude = amp, alpha = alpha, exptime = exptime_FB)
        flux_FB = mag_to_flux(flux_to_mag(flux_FB, zp = ZP)-DM)
        flux_FB[np.isnan(flux_FB)] = 0
        spl_DEI = spl_allfilt_DEI[filter_]
        flux_DEI = mag_to_flux(spl_DEI(data_filt['obsdate']))
        flux_both = flux_FB+ flux_DEI
        mag_model = flux_to_mag(flux_FB, zp = ZP)
        mag_DOM = flux_to_mag(flux_DEI, zp = ZP)
        mag_both = flux_to_mag(flux_both, zp = ZP)
        residual[filter_] = data_filt['mag']-DM - mag_both
        phase[filter_] = data_filt['obsdate']
#%% Residual
for filter_ in residual.keys():
    residual_ = residual[filter_]
    phase_ = phase[filter_]
    plt.figure(figsize = (6,2), dpi = 300)
    plt.axhline(0, linestyle = '--', c= 'k')
    plt.scatter(phase_, residual_, c = colors[filter_])
    plt.ylim(-1, 1)
    plt.xlim(59529.3315-2.5, 59529.3315+10.1)
