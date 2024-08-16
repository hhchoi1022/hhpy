
#%%
import glob
from astropy.table import Table
from datetime import datetime
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.units as u
import numpy as np
from astropy.time import Time
from timeseriesspectrum import TimeSeriesSpectrum
#%%
class WiseRepSpectrum(TimeSeriesSpectrum):
    '''
    #### Timeseries observed spectrum module for WISeREP data
    Paramters
    =========
    1. specfilekey : str = abs path of the spectrum files (only ascii txt format)
    
    Methods
    =========
    1. get_spec_date(time) = return Spectrum instance at the time (relative to the first spectrum observing date)
    2. show_spec_date(time, 
                      show_flux_unit : str = 'fnu',
                      smooth_factor : int = 1,
                      normalize : bool = False,
                      normalize_cenwl : int = 6700,
                      normalize_wlwidth : int = 20) 
                      
                      = show the spectrum at the time (relative to the first spectrum observing date)
    
    Notes
    =========
                
    '''
    def __init__(self,
                 specfilekey : str,
                 format_ = 'ascii'):

        self._format = format_
        self._specfilekey = specfilekey
        self.filelist, self._phaselist = self._get_filelist()
        phase_array, wavelength_array, flamb_array = self._timeSeriesSpectrum()
        TimeSeriesSpectrum.__init__(self, time_array = phase_array, wavelength_array= wavelength_array, flux_array= flamb_array, flux_unit= 'flamb')

    
    def _get_obsdate(self,
                    filename : str,
                    date_pattern : str = '(\d\d\d\d)-(\d\d)-(\d\d)_(\d\d)-(\d\d)',
                    date_pattern_2 : str = '(245\d\d\d\d.\d\d\d)'):
        import re
        try:
            year, month, day, hour, minute = re.search(date_pattern, filename).groups()
        except:
            jd, = re.search(date_pattern_2, filename).groups()
            dt = Time(jd, format= 'jd').datetime
            year, month, day, hour, minute = dt.year, dt.month, dt.day, dt.hour, dt.minute
        return int(year), int(month), int(day), int(hour), int(minute)

    def _get_filelist(self):
        return self._sort_filelist(glob.glob(self._specfilekey))
    
    def _get_spec_data(self,
                      filename : str):
        data_tbl = ascii.read(filename)
        if len(data_tbl.keys()) == 2:
            data_tbl.rename_columns(data_tbl.keys(), ['wavelength[AA]', 'flamb'])
        else:
            data_tbl.rename_columns(data_tbl.keys(), ['wavelength[AA]', 'flamb', 'e_flamb'])
        return data_tbl
    
    def _sort_filelist(self,
                       filelist):
        
        datelist = []
        for specfile in filelist:
            pattern = '(\d\d\d\d)-(\d\d)-(\d\d)_(\d\d)-(\d\d)'
            year, month, day, hour, minute = self._get_obsdate(specfile, date_pattern= pattern)
            date = datetime(year = year, month = month, day= day, hour = hour, minute = minute)
            datelist.append(date)

        tmp_tbl = Table()
        tmp_tbl['date'] = datelist
        tmp_tbl['file'] = filelist
        tmp_tbl.sort('date')
        phaselist = [round(Time(date, format = 'datetime').mjd,3) for date in tmp_tbl['date']]
        sorted_filelist = tmp_tbl['file']
        return sorted_filelist, phaselist

    def _timeSeriesSpectrum(self):
        wavelengthlist = []
        flamblist = []
        for file_ in self.filelist:
            data_tbl = self._get_spec_data(file_)
            wavelength = list(data_tbl['wavelength[AA]'])
            flamb = list(data_tbl['flamb'])
            wavelengthlist.append(wavelength)
            flamblist.append(flamb)
        phase_array = np.array(self._phaselist)
        wavelength_array = np.array(wavelengthlist, dtype = 'object')
        flamb_array = np.array(flamblist, dtype = 'object')
        return phase_array, wavelength_array, flamb_array
        
            
# %%

#%%
if __name__ == '__main__':
    specfolder = '/Users/hhchoi1022/Gitrepo/data/SN2018oh/observation/spectrum/WISeREP/ascii/*'
    
    A = WiseRepSpectrum(specfolder)
    #%%
    for phase in A.phase[:2]:
        
        A.show_spec_date(phase, color = None, normalize = True, label = phase, show_flux_unit= 'jy', smooth_factor= 11)
        plt.xlim(3000, 9000)
        #plt.ylim(-1, 5)
        plt.legend()
    # %%
