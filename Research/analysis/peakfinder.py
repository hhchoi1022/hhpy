
#%%
from Research.analysis.observedphot import ObservedPhot
import numpy as np
import matplotlib.pyplot as plt
#%%
class PeakFinder:
    
    def __init__(self, obs_phot : ObservedPhot):
        self.observedphot = obs_phot
        self.data_by_filter = self.observedphot.get_filt_data()
        self.data_for_fitting = self.data_by_filter
    
    def find_peak_for_all_filters(self,
                                  poly_fit_degree = 4):
        for filter_, data in self.data_by_filter.items():
            self.define_range_mjd_for_peak(mjd_min = data[np.argmin(data['mag'])]['obsdate']-10, mjd_max = data[np.argmin(data['mag'])]['obsdate']+10)
            plt.figure(figsize=(10, 6))
            plt.title('Filter: %s'%filter_)
            self.find_peak(mjdlist = np.array(self.data_for_fitting[filter_]['obsdate']), maglist = np.array(self.data_for_fitting[filter_]['mag']), magerrlist = np.array(self.data_for_fitting[filter_]['e_mag']), poly_fit_degree = 4)
    
    def define_range_mjd_for_peak(self,
                                  mjd_min : float = None,
                                  mjd_max : float = None
                                  ):
        if (mjd_min is None) or (mjd_max is None):
            mjd_peak = np.mean([tbl_filter[np.argmin(tbl_filter['mag'])]['obsdate'] for tbl_filter in list(self.data_by_filter.values())])
            mjd_min = mjd_peak - 10
            mjd_max = mjd_peak + 10
        for filter_, data in self.data_by_filter.items():
            self.data_for_fitting[filter_] = data[(data['obsdate'] > mjd_min) & (data['obsdate'] < mjd_max)]
    
    def find_peak(self,
                  mjdlist,
                  maglist,
                  magerrlist = None,
                  poly_fit_degree = 4):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.optimize import fmin

        # Scale time data
        time_min, time_max = mjdlist.min(), mjdlist.max()
        time_shifted = mjdlist - time_min

        # Fit the polynomial to scaled time data
        if magerrlist is not None:
            w = 1 / np.array(magerrlist)
        coefficients = np.polyfit(time_shifted, maglist, poly_fit_degree, w = None)
        poly_function = np.poly1d(coefficients)

        # Plot
        mjd_axis = np.linspace(0, np.max(time_shifted), 500)
        mag_axis = poly_function(mjd_axis)

        # Find the peak time in scaled coordinates
        mag_model = poly_function(mjd_axis)
        peak_idx = np.argmin(mag_model)
        peak_mag = mag_model[peak_idx]
        peak_time = mjd_axis[peak_idx] + time_min
        
        plt.plot(mjdlist, maglist, 'o', label='Data', markersize=8)
        plt.plot(mjd_axis + time_min, mag_axis, '-', label='Fitted Polynomial (Scaled)', linewidth=2)
        if magerrlist is not None:
            plt.errorbar(mjdlist, maglist, yerr = magerrlist, fmt = 'None', elinewidth = 1, ecolor = 'k', capsize = 2)

        plt.plot(peak_time, peak_mag, 'ro', label=f'Peak ({peak_time:.2f}, {peak_mag:.2f})', markersize=10)

        plt.xlabel('Obsdate (MJD)')
        plt.ylabel('Magnitude (AB)')
        plt.legend()
        plt.grid(True)
        plt.show()

        
    
# %%
from astropy.io import ascii
tbl = ascii.read('/mnt/data1/supernova_rawdata/SN2021aefx/photometry/all_phot_MW_dereddening_Host_dereddening.dat', format = 'fixed_width')
obs = ObservedPhot(tbl, target = 'SN2021aefx')
obs.exclude_observatory(['DLT40','LasCumbres0.4m','Swift'])
# %%
P = PeakFinder(obs)
#%%
P.find_peak_for_all_filters()
# %%
