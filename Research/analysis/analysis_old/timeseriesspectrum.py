
#%%

from astropy.io import ascii
import matplotlib.pyplot as plt
import glob
import astropy.units as u
from astropy.table import Table
from spectrum import Spectrum
import numpy as np
#%%

class TimeSeriesSpectrum:
    '''
    Paramters
    =========
    1. time_array : array = time array, must be 1 dimensional
    2. wavelength_array : array = time x wavelength array, must be 2 dimensional. For example, if 500 wave coordinates with 3 days, 3x500 array must be input shape
    3. flux_array : array = time x flux array, must have same shape as wavelength array
    4. flux_unit : str = select among [fnu, flamb, jy], input flux unit
    
    Methods
    =========
    1. get_spec_date(time) = return Spectrum instance at the time
    2. show_spec_date(time, 
                      show_flux_unit : str = 'fnu',
                      smooth_factor : int = 1,
                      normalize : bool = False,
                      normalize_cenwl : int = 6700,
                      normalize_wlwidth : int = 20) 
                      
                      = show the spectrum at the time
                
    '''
    
    def __init__(self,
                 time_array,
                 wavelength_array,
                 flux_array,
                 flux_unit : str = 'fnu',
                 phot_filterset = 'UBVRIugri'
                 ):
        
        self._flux_array = flux_array
        self._flux_unit = flux_unit
        self._wavelength_array = wavelength_array
        
        self.phase = time_array
        self.data = self._data()
        self.spec = self._all_spec()
        self.phot = self._all_phot(filterset  =phot_filterset)
        
        if not (any([length == len(self.phase) for length in self._flux_array.shape])) & (self._flux_array.shape == self._wavelength_array.shape):
            raise ValueError('Dimension not matching Dim(Phase):%d, Dim(wl):%dx%d, Dim(flux):%dx%d'%(len(self.phase), self._wavelength_array.shape[0], self._wavelength_array.shape[1], self._flux_array.shape[0], self._flux_array.shape[1]))
    
    def _data(self):
        all_data = dict()
        for phase, wavelength, flux in zip(self.phase, self._wavelength_array, self._flux_array):
            all_data[phase] = dict()
            all_data[phase]['wavelength'] = wavelength
            all_data[phase]['flux'] = flux
        return all_data
    
    def _all_spec(self):
        all_spectrum = dict()
        for phase in self.phase:
            wavelength = self.data[phase]['wavelength']
            flux = self.data[phase]['flux']
            spec = Spectrum(wavelength= wavelength, flux = flux, flux_unit = self._flux_unit)
            all_spectrum[phase] = spec
        return all_spectrum
    
    def _all_phot(self,
                  filterset : str = 'UBVRIugri'):
        phot_tbl = Table()
        phot_tbl['MJD'] = list(self.spec.keys())
        for filter_ in filterset:
            phot_tbl[filter_] = 0.0
        for i, spec in enumerate(list(self.spec.values())):
            photometry = spec.photometry
            for filter_ in filterset:
                phot_tbl[filter_][i] = photometry[filter_]
        return phot_tbl
    
    def phot_formatted(self,
                       phot_tbl,
                       filterset : str ='UBVRIugri'):
        formatted_tbl = Table(names = ['obsdate','filter','mag','e_mag','observatory','status','[mag_cor-mag_orign]'], dtype =[float, str, float, float, str, str, float])
        header = phot_tbl.colnames
        for data in phot_tbl:
            for filter_ in filterset:
                if filter_ in header:
                    formatted_tbl.add_row([data['MJD'], filter_, data[filter_], 0.00, 'Spectrum', 'detected', 0.0])
        return formatted_tbl

    
    def get_spec_date(self,
                      time):
        exact_time = self.phase[np.argmin(np.abs(time-self.phase))]
        if not exact_time == time:
            print('time %.4f is not found in TimeSeriesSpectrum %.4f is selected instead.' %(time, exact_time))
        return self.spec[exact_time]

    def show_spec_date(self,
                       time,
                       show_flux_unit : str = 'fnu',
                       smooth_factor : int = 1,
                       normalize : bool = False,
                       normalize_cenwl : int = 6700,
                       normalize_wlwidth : int = 20,
                       label : str = '',
                       color : str = None
                       ):
        spec = self.get_spec_date(time = time)
        spec.show(show_flux_unit= show_flux_unit, smooth_factor= smooth_factor, normalize= normalize, normalize_cenwl= normalize_cenwl, normalize_wlwidth= normalize_wlwidth, label = label, color = color)


# %% Tutorial
if __name__ == '__main__':
    from HHsupport_analysis import read_Polin2019_spec
    specfolder = '/Users/hhchoi1022/Gitrepo/Data/IaSNe_Model/spectrum/Polin/ddet_Polin2019/*.h5'
    file_ = glob.glob(specfolder)[30]
    Lnu, Llamb, lamb, time, mu =  read_Polin2019_spec(file_)
    lamb = np.tile(lamb, [len(time),1])
    #lamb = lamb * np.ones(Lnu.shape)
    A = TimeSeriesSpectrum(time_array = time, wavelength_array=lamb, flux_array=Lnu)
#%%
if __name__ == '__main__':
    speckey = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/spectrum/WISeREP/ascii/*.txt'
    data = ascii.read(glob.glob(specfolder)[0])
    from wiserepSpectrum import WiseRepSpectrum
    A  = WiseRepSpectrum(specfilekey = speckey)
#%%
    A.show_spectrum(A.filelist[0:5], fnu = False, normalize = False, logscale = True, median_filter = True, median_filtersize = 11)

# %%
