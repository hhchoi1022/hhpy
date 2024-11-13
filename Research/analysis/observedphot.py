
#%%
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from Research.helper import Helper
from astropy.time import Time
#%%

class ObservedPhot:
    """
    ObservedData class for photometric observed data

    Parameters
    ----------
    data_tbl : str
        The input data table with the columns of:
        - 'obsdate' : float or int
            Observation date in MJD
        - 'filter' : str
            Filter name (e.g., 'g', 'r', 'i', 'B', 'V', and etc.)
        - 'mag' : float or int
            Apparent magnitude
        - 'e_mag' : float or int
            Error in apparent magnitude
        - 'observatory' : str
            Observatory code (e.g., 'KGEO', 'LCO', and etc.)
        - 'status' : 'detected' or 'UL'
            Detection status: 
                -- 'detected' for detected observation, 
                -- 'UL' for the upper limit observation

    target : str, optional
        Target name, by default 'SN2021aefx'
    observer : str, optional
        Name of observer, by default None
    MW_extinction_corrected : bool, optional
        Flag for MW extinction correction. If True, object
        is corrected for MW extinction according to Schlegel, 
        Finkbeiner, & Davis (1998) using the Schlafly & 
        Finkbeiner (2011) recalibration. 
        By default, `MW_extinction_corrected` is False.
    Host_extinction_corrected : bool, optional
        Flag for host extinction correction. By default, 
        `Host_extinction_corrected` is False

    Methods
    -------
    get_observatory()
        Returns observatory name as a sorted list
    get_all_filter()
        Returns filter name as a sorted list
    get_defined_filter()
        Returns filter name as a sorted list if there's at least one detection
    get_data_detected(key='detected')
        Returns a table with rows that matched the input key
    get_data_ul(key='UL')
        Returns a table with rows that matched the input key
    get_marker_observatory()
        Returns a dictionary for observatory marker
    get_color_filter()
        Returns a dictionary for color filter
    get_label_filter()
        Returns a dictionary for filter label 
    get_offset_filter()
        Returns a dictionary for filter offset
    get_filt_data(data_tbl, filters='UBVugri')
        Returns a dictionary of tables with rows that matched the input filter
    exclude_observatory(observatory)
        Exclude row(s) from the data with given observatory
    exclude_filter(filter_)
        Exclude row(s) from the data with given filter
    show_lightcurve(dpi=300, figsize=(8,6), scatter_size=30,
                     scatter_linewidth=0.5, errorbar_linewidth=0.5,
                     errorbar_capsize=3, line=False, label=False,
                     label_location='upper right', color_BV=True,
                     color_gr=True)
        Plot a lightcurve with B-V, and g-r color

    """
    def __init__(self,
                 data_tbl : str,
                 target = 'SN2021aefx'
                 ):
        
        self.helper = Helper()
        self.target = target
        self.data = data_tbl
        
    def __repr__(self) -> str:
        return (f'(target={self.target},'
                f' observatory={self.get_observatory()},'
                f' filter={self.get_defined_filter()}')

    def get_observatory(self):
        return sorted(list(set(self.data['observatory'])))
    
    def get_all_filter(self):
        filt_array = np.array(list(self.helper.load_filt_keys()[1].keys()))
        for filt_ in list(self.helper.load_filt_keys()[1].keys()):
           if not filt_ in list(set(self.data['filter'])):
                filt_array = np.delete(filt_array, np.where(filt_array == filt_))
        return list(filt_array)
        #return sorted(list(set(self.data['filter'])))
    
    def get_defined_filter(self):
        return sorted(list(self.get_color_filter().keys()))
    
    def get_data_detected(self):
        return self.data[np.array([True if str(value).upper() == 'TRUE' else False for value in self.data['detected']])]
    
    def get_data_ul(self):
        return self.data[~ np.array([True if str(value).upper() == 'TRUE' else False for value in self.data['detected']])]
    
    def get_marker_observatory(self):
        observatories = self.get_observatory()
        markers = self.helper.load_marker_keys(len(observatories))
        marker_keys = dict()
        for observatory, marker in zip(observatories, markers):
            marker_keys[observatory] = marker
        return marker_keys
        
    def get_color_filter(self):
        filters = self.get_all_filter()
        colors, _, _, _, _ = self.helper.load_filt_keys(filters)
        return colors
    
    def get_label_filter(self):
        filters = self.get_all_filter()
        _, _, _, _, labels = self.helper.load_filt_keys(filters)
        return labels
    
    def get_offset_filter(self):
        filters = self.get_all_filter()
        _, offsets, _, _, _ = self.helper.load_filt_keys(filters)
        return offsets
    
    def get_filt_data(self,
                      filters : str = None):
        all_filt_data = dict()
        if filters == None:
            filters = self.get_defined_filter()
        for filter_ in filters:
            all_filt_data[filter_] = self.data[self.data['filter'] == filter_]
        return all_filt_data
    
    def exclude_observatory(self,
                            observatory : str):
        removed_tbl = self.helper.remove_rows_table(self.data, column_key = 'observatory', remove_keys = observatory)
        self.data = removed_tbl
    
    def exclude_filter(self,
                       filter_ : str):        
        removed_tbl = self.helper.remove_rows_table(self.data, column_key = 'filter', remove_keys = filter_)
        self.data = removed_tbl
    
    def show_lightcurve(self,
                        dpi = 300,
                        figsize = (12,8),
                        phase_format = 'datetime',
                        phase_min = None,
                        phase_max = None,
                        phase_binsize : int = 5,
                        mag_min = None,
                        mag_max = None,
                        mag_binsize : float = 2,
                        scatter_size = 20,
                        scatter_linewidth = 0.02,
                        errorbar_linewidth = 0.02,
                        errorbar_capsize = 0.02,
                        scatter_alpha : float = 0.4,
                        line : bool = False,
                        label : bool = False,
                        label_location = 'upper right',
                        color_UB : bool = False,
                        color_BV : bool = False,
                        color_ug : bool = False,
                        color_gr : bool = False,
                        UL : bool = False,
                        UL_linewidth_ver : float = 0.03,
                        UL_linewidth_hor : float = 0.03,
                        UL_linelength_ver : float = 0.01,
                        UL_linelength_hor : float = 0.02,
                        UL_headlength : float = 0.015,
                        UL_headwidth : float = 0.03,
                        UL_alpha : float = 0.5):
        import matplotlib.patches as mpatches
        
        # Set x axis
        if phase_min is None:
            phase_min = np.min(self.get_data_detected()['obsdate'])-10
        if phase_max is None:
            phase_max = np.max(self.get_data_detected()['obsdate'])+10
        if phase_format.upper() == 'DATETIME':
            dates = [Time(t, format='mjd').datetime.strftime('%Y-%m-%d') for t in np.arange(phase_min, phase_max, phase_binsize)]
        else:
            dates = np.arange(phase_min, phase_max, phase_binsize).astype(int)
        if mag_min is None:
            mag_min = int(np.min(self.get_data_detected()['mag']))-1
        if mag_max is None:
            mag_max = int(np.max(self.get_data_detected()['mag']))+1
        mags = np.arange(mag_min, mag_max, mag_binsize)

        # UL
        data_detected = self.get_data_detected()
        data_ul = self.get_data_ul()
        colors = self.get_color_filter()
        labels = self.get_label_filter()
        offsets = self.get_offset_filter()
        markers = self.get_marker_observatory()
        
        if color_BV | color_gr | color_UB:
            plt.subplots_adjust(hspace=0)
            
            plt.tight_layout(h_pad=-1, w_pad=-1) 
            plt.gca().invert_yaxis()
            gs = gridspec.GridSpec(nrows = 2, ncols =1, height_ratios= [10, 3], width_ratios = [6])
            
            #LC
            ax1 = plt.subplot(gs[0])
            for filter_ in colors.keys():
                for observatory in self.get_observatory():
                    show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
                    show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
                    if len(show_data_detected) > 0:
                        if line:
                            ax1.plot(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], color = colors[filter_], linewidth = scatter_linewidth)
                        else:
                            ax1.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], show_data_detected['e_mag'] , fmt = 'none', elinewidth = errorbar_linewidth, capsize = errorbar_capsize, c = 'k', capthick = errorbar_linewidth, zorder = 1)
                            if observatory == 'KCT':
                                ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, label = labels[filter_], zorder = 5)
                                ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 5)
                            else:
                                ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, label = labels[filter_], zorder = 2)
                                ax1.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 3)
                    if UL:
                        if len(show_data_ul) > 0:
                            for x,y in zip(show_data_ul['obsdate'], show_data_ul['depth_5sig']+offsets[filter_]):
                                ax1.arrow(x=x,y=y,dx=0,dy=UL_linelength_ver, linewidth = UL_linewidth_ver, head_width = UL_headwidth, head_length = UL_headlength, color = colors[filter_], shape = 'full', alpha = UL_alpha)
                                ax1.arrow(x=x-UL_linelength_hor,y=y,dx=2*UL_linelength_hor,dy=0,linewidth = UL_linewidth_hor, head_width = 0, color = colors[filter_], alpha = UL_alpha)
            ax1.set_ylabel('Apparent magnitude [AB]')
            ax1.invert_yaxis()
            ax1.set_xlim(phase_min, phase_max)
            ax1.set_yticks(mags, mags)
            ax1.set_ylim(mag_max, mag_min)
            
            #Color
            ax2 = plt.subplot(gs[1], sharex = ax1)
            
            filt_data = self.get_filt_data(data_detected)
            if color_UB:
                UB_tbl = self.helper.match_table(filt_data['U'], filt_data['B'], key = 'obsdate', tolerance = 1)
                if len(UB_tbl) > 0:
                    UB_tbl['e_color'] = (UB_tbl['e_mag_1'] **2 + UB_tbl['e_mag_2'] **2)**(1/2)
                    ax2.errorbar(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2']-0.5,  UB_tbl['e_color'] , fmt = 'none', elinewidth = errorbar_linewidth, capsize = errorbar_capsize, c='k', capthick = errorbar_linewidth, zorder = 1)
                    ax2.scatter(UB_tbl['obsdate_1'], UB_tbl['mag_1']-UB_tbl['mag_2']-0.5, label = 'U-B-0.5', facecolors = 'cyan',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 2)
            if color_BV:
                BV_tbl = self.helper.match_table(filt_data['B'], filt_data['V'], key = 'obsdate', tolerance = 1)
                if len(BV_tbl) > 0:
                    BV_tbl['e_color'] = (BV_tbl['e_mag_1'] **2 + BV_tbl['e_mag_2'] **2)**(1/2)
                    ax2.errorbar(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2']+0.5,  BV_tbl['e_color'] , fmt = 'none', elinewidth = errorbar_linewidth, capsize = errorbar_capsize, c='k', capthick = errorbar_linewidth, zorder = 1)
                    ax2.scatter(BV_tbl['obsdate_1'], BV_tbl['mag_1']-BV_tbl['mag_2']+0.5, label = 'B-V+0.5', facecolors = 'b',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 2)
            if color_gr:
                gr_tbl = self.helper.match_table(filt_data['g'], filt_data['r'], key = 'obsdate', tolerance = 1)
                if len(gr_tbl) > 0:
                    gr_tbl['e_color'] = (gr_tbl['e_mag_1'] **2 + gr_tbl['e_mag_2'] **2)**(1/2)
                    ax2.errorbar(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'],  gr_tbl['e_color'] , fmt = 'none', elinewidth = errorbar_linewidth, capsize = errorbar_capsize, c = 'k', capthick = errorbar_linewidth, zorder = 1)
                    ax2.scatter(gr_tbl['obsdate_1'], gr_tbl['mag_1']-gr_tbl['mag_2'], label = 'g-r', facecolors = 'g',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 2)
            if color_ug:
                ug_tbl = self.helper.match_table(filt_data['u'], filt_data['g'], key = 'obsdate', tolerance = 1)
                if len(ug_tbl) > 0:
                    ug_tbl['e_color'] = (ug_tbl['e_mag_1'] **2 + ug_tbl['e_mag_2'] **2)**(1/2)
                    ax2.errorbar(ug_tbl['obsdate_1'], ug_tbl['mag_1']-ug_tbl['mag_2'],  ug_tbl['e_color'] , fmt = 'none', elinewidth = errorbar_linewidth, capsize = errorbar_capsize, c = 'k', capthick = errorbar_linewidth, zorder = 1)
                    ax2.scatter(ug_tbl['obsdate_1'], ug_tbl['mag_1']-ug_tbl['mag_2'], label = 'u-g', facecolors = 'magenta',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 2)

            plt.xlabel('Phase [days]')
            ax2.set_ylabel('Color')
            ax2.set_xticks(np.arange(phase_min, phase_max, phase_binsize), dates, rotation = 45)
            ax2.set_ylim(-1.5, 2)
            ax2.set_xlim(phase_min, phase_max)
            
            # legends
            sorted_keys = {k: v for k, v in labels.items()}
            rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
            name_rows = [labels[clr] for clr in sorted_keys]
            columns = [plt.plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in self.get_observatory()]
            name_columns = self.get_observatory()
            if label:
                ax1.legend(rows + columns, name_rows + name_columns, loc=label_location, ncol = 2)
                ax2.legend(loc = label_location)
            
            return ax1, ax2
        
        else:
            #LC
            plt.gca().invert_yaxis()
            for filter_ in colors.keys():
                for observatory in self.get_observatory():
                    show_data_detected = data_detected[(data_detected['filter'] == filter_) & (data_detected['observatory'] == observatory)]
                    show_data_ul = data_ul[(data_ul['filter'] == filter_) & (data_ul['observatory'] == observatory)]
                    if len(show_data_detected) > 0:
                        if line:
                            plt.plot(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], color = colors[filter_], linewidth = scatter_linewidth)
                        else:
                            plt.errorbar(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], show_data_detected['e_mag'] , fmt = 'none', elinewidth = errorbar_linewidth, capsize = errorbar_capsize, c = 'k', capthick = errorbar_linewidth, zorder = 1)
                            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = 'none',  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, label = labels[filter_], zorder = 2)
                            plt.scatter(show_data_detected['obsdate'], show_data_detected['mag']+offsets[filter_], marker = markers[observatory], facecolors = colors[filter_],  edgecolors = 'k', s = scatter_size, linewidth = scatter_linewidth, alpha = scatter_alpha, zorder = 3)
                    if UL:
                        if len(show_data_ul) > 0:
                            for x,y in zip(show_data_ul['obsdate'], show_data_ul['depth_5sig']+offsets[filter_]):
                                plt.arrow(x=x,y=y,dx=0,dy=UL_linelength_ver, linewidth = UL_linewidth_ver, head_width = UL_headwidth, head_length = UL_headlength, color = colors[filter_], shape = 'full', alpha = UL_alpha)
                                plt.arrow(x=x-UL_linelength_hor,y=y,dx=2*UL_linelength_hor,dy=0,linewidth = UL_linewidth_hor, head_width = 0, color = colors[filter_], alpha = UL_alpha)
            plt.xticks(np.arange(phase_min, phase_max, phase_binsize), dates, rotation = 45)
            plt.xlim(phase_min,phase_max)
            plt.yticks(mags, mags)
            plt.ylim(mag_max, mag_min)
            # legends
            #sorted_keys = {k: v for k, v in sorted(labels.items(), key=lambda item: item[1])}
            sorted_keys = {k: v for k, v in labels.items()}
            rows = [mpatches.Patch(color=colors[clr]) for clr in sorted_keys]
            name_rows = [labels[clr] for clr in sorted_keys]
            columns = [plt.plot([],[], markers[observatory], markerfacecolor ='w', markeredgecolor='k')[0] for observatory in self.get_observatory()]
            name_columns = self.get_observatory()
            if label:
                plt.legend(rows + columns, name_rows + name_columns, loc=label_location, ncol = 2)
            
            plt.xlabel('Phase [days]')
            plt.ylabel('Apparent Magnitude [AB]')
            
        


# %%
if __name__ =='__main__':    
    plt.figure(dpi = 400, figsize =(9,7))
    filepath_phot = '/data1/supernova_rawdata/SN2021aefx/photometry/all_phot_MW_dereddening_Host_dereddening.dat'
    filepath_spec = '/data1/supernova_rawdata/SN2021aefx/photometry/all_spec_MW_dereddening_Host_dereddening.dat'

    #filepath_all = '/data7/yunyi/temp_supernova/Gitrepo/Research/analysis/all_phot_MW_dereddening_Host_dereddening.dat'

    tbl_phot = ascii.read(filepath_phot, format = 'fixed_width')
    tbl_spec = ascii.read(filepath_spec, format = 'fixed_width')
    obs_phot = ObservedPhot(tbl_phot)
    obs_spec = ObservedPhot(tbl_spec)
    #observed_data.exclude_observatory('Swope')

    ax1, ax2 = obs_phot.show_lightcurve(UL = True, scatter_alpha = 1, UL_headlength=0.2, UL_headwidth=2,  UL_linelength_ver=0.3, scatter_size= 30, errorbar_capsize=0, UL_linewidth_hor=0.5, UL_linewidth_ver=0.8, UL_linelength_hor=1, label =True, color_BV = True, color_gr = True, color_UB = True, color_ug = False)
    #ax3, ax4 = obs_spec.show_lightcurve(UL = True, scatter_alpha = 1, UL_headlength=0.2, UL_headwidth=2,  UL_linelength_ver=0.3, scatter_size= 30, errorbar_capsize=0, UL_linewidth_hor=0.5, UL_linewidth_ver=0.8, UL_linelength_hor=1, label =True, color_BV = True, color_gr = True, color_UB = True, color_ug = False)
    ax1.set_ylim(22, 6)
    plt.xlim(59525, 59540)
    plt.show()
# %%
if __name__ == '__main__':
    from astropy.table import vstack
    filepath_imsng = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/all_IMSNG.dat'
    filepath_atlas = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/ATLAS_new.dat'
    filepath_7DT = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/all_7DT.dat'
    
    tbl_IMSNG = ascii.read(filepath_imsng, format = 'fixed_width')
    tbl_ATLAS = ascii.read(filepath_atlas, format = 'fixed_width')
    tbl_7DT = ascii.read(filepath_7DT, format = 'fixed_width')
    tbl_all = vstack([tbl_IMSNG, tbl_ATLAS, tbl_7DT])
    obs_all = ObservedPhot(tbl_all)
    obs_all.show_lightcurve(dpi = 300,
                            figsize = (12,8),
                            phase_format = 'datetime',
                            phase_min = 60180,
                            phase_max = 60390,
                            phase_binsize  = 50,
                            mag_min = 12,
                            mag_max = 22,
                            mag_binsize  = 2,
                            scatter_size = 60,
                            scatter_linewidth = 0.4,
                            errorbar_linewidth = 0.02,
                            errorbar_capsize = 0.02,
                            scatter_alpha = 1,
                            line = False,
                            label = True,
                            label_location = 'upper right',
                            color_UB = False,
                            color_BV = False,
                            color_ug = False,
                            color_gr = False,
                            UL = True,
                            UL_linewidth_ver = 1,
                            UL_linewidth_hor = 1,
                            UL_linelength_ver = 0.5,
                            UL_linelength_hor = 2.5,
                            UL_headlength = 0.2,
                            UL_headwidth = 3,
                            UL_alpha = 0.5)
                    

    
    
    
        
        
        
# %% Tutorial
if __name__ == '__main__':
    
    #filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostmwextin3.10.dat'
    #filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostmwextin3.10.dat'
    #filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostmwextin3.10.dat'
    
    #filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_hostextin3.10.dat'
    #filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_hostextin3.10.dat'
    #filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_hostextin3.10.dat'
    
    filepath_1 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Hosseinzadeh2022/Hosseinzadeh2022_noextin.dat'
    filepath_2 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/IMSNG/IMSNG_noextin.dat'
    filepath_3 = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Ashall2022/Ashall2022_noextin.dat'
    filepath_all = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/lightcurve/Alldata_No_cor.dat'
    
    tbl1 = ascii.read(filepath_1, format = 'fixed_width')
    tbl2 = ascii.read(filepath_2, format = 'fixed_width')
    tbl3 = ascii.read(filepath_3, format = 'fixed_width')
    tbl_all = ascii.read(filepath_all, format = 'fixed_width')
    from wiserepSpectrum import WiseRepSpectrum
    spec = WiseRepSpectrum(specfilekey = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/spectrum/WISeREP/ascii/*FTS*')
    spec_phot = spec.phot
    #spec_phot.remove_row(index = 0)
    
    tbl4 = spec.phot_formatted(spec_phot)
#%%
if __name__ =='__main__':    
    #tbl3 = ascii.read(filepath_3, format = 'fixed_width')
    from astropy.table import vstack
    tbl_obs = vstack([tbl_all])
    observed_data = ObservedPhot(tbl_obs, MW_extinction_corrected= True, Host_extinction_corrected= True)
    observed_data.get_color_filter()
    plt.figure(dpi = 500, figsize = (8,6))
    observed_spec = ObservedPhot(tbl4)
    observed_spec.show_lightcurve(scatter_size = 80)
    observed_data.show_lightcurve(label = True, label_location=0, UL= True, UL_headwidth=0.3, UL_linelength_ver=0.5, UL_alpha = 0.3, color_BV = False, color_gr = False, phase_binsize= 5)
    
    plt.xlim(59525, 59540)
    plt.ylim(24, 5)
    plt.figure(dpi = 300)
    from matplotlib import cm
    cmap = cm.get_cmap('hsv')
    num_spec = 4
    for i in range(num_spec):
        idx = list(spec.spec.keys())[i]
        spec.show_spec_date(idx, show_flux_unit= 'fnu', normalize = False, normalize_cenwl=6500, color = cmap((i+1)/num_spec), label = f'{i+1} th spectrum')
        plt.legend()
        