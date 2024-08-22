
#%%

from astropy.io import ascii
import matplotlib.pyplot as plt
import glob
import astropy.units as u
from astropy.table import Table
import numpy as np

from Research.spectroscopy import Spectrum
from Research.spectroscopy import SpectroscopyFile
#%%

class TimeSeriesSpectrum:
    def __init__(self,
                 specfilelist : str,
                 flux_unit : str = 'flamb'
                 ):
        self._filelist = specfilelist
        self._obsdate_data, self._wavelength_data, self._flux_data, self._headers, self._filelist = self._read_specfilelist()
        if flux_unit.upper() in ['FNU', 'FLAMB', 'JY']:
           self._flux_unit = flux_unit
        else:
           raise ValueError('%s is not registered as spectral flux unit'%flux_unit)
        
        self.spectrum = self._get_all_spec()
    @property
    def data(self):
        spec_datalist = list()
        for obsdate, wavelength, flux in zip(self._obsdate_data, self._wavelength_data, self._flux_data):
            spec_data = dict()
            spec_data['obsdate'] = obsdate
            spec_data['wavelength'] = wavelength
            spec_data['flux'] = flux
            spec_datalist.append(spec_data)
        return spec_datalist

    def _read_specfilelist(self):
        obsdate_data = list()
        wavelength_data = list()
        flux_data = list()
        headerlist = list()
        filenamelist = list()
        for file_ in self._filelist:
            specfile = SpectroscopyFile(file_)
            obsdate_data.append(specfile.obsdate)
            wavelength_data.append(specfile.wavelength)
            flux_data.append(specfile.flux)
            headerlist.append(specfile.header)
            filenamelist.append(file_)
        return obsdate_data, wavelength_data, flux_data, headerlist, filenamelist
    
    def _get_all_spec(self):
        all_spectrum = list()
        for _obsdate_data, wavelength, flux in zip(self._obsdate_data, self._wavelength_data, self._flux_data):
            spec = Spectrum(wavelength= wavelength, flux = flux, flux_unit = self._flux_unit)
            all_spectrum.append(spec)
        return all_spectrum
    
    def photometry(self,
                   filterset : str = 'UBVRIugri',
                   save : bool = False):
        phot_tbl = Table()
        phot_tbl['obsdate'] = self._obsdate_data
        for filter_ in filterset:
            phot_tbl[filter_] = 0.0
            phot_tbl['observatory'] = '                 '
            phot_tbl['file'] = '                                                                                                                                                                                                                                                                                       '
        for i, (spec,hdr,file_) in enumerate(zip(self.spectrum,self._headers, self._filelist)):
            # set observatory 
            if 'OBSERVER' in hdr.keys():
                observatory = hdr['OBSERVER']
            if 'INSTRUME' in hdr.keys():
                observatory = hdr['INSTRUME']
            if 'OBSERVAT' in hdr.keys():
                observatory = hdr['OBSERVAT']
            photometry = spec.photometry(filterset = filterset)
            for filter_ in filterset:
                phot_tbl[filter_][i] = photometry[filter_]
                phot_tbl['observatory'][i] = observatory
                phot_tbl['file'][i] = file_
            
        formatted_tbl = Table(names = ['obsdate', 'mag', 'e_mag', 'magsys', 'filter', 'depth_5sig', 'zp', 'observatory','detected', 'filename'], dtype =[float, float, float, str, str, float, float, str, bool, str])
        header = phot_tbl.colnames
        for data in phot_tbl:
            for filter_ in filterset:
                if filter_ in 'UBVRI':
                    formatted_tbl.add_row([data['obsdate'], data[filter_], 0.05, 'Vega', filter_ , None, None, data['observatory'], True, data['file']])
                else:
                    formatted_tbl.add_row([data['obsdate'], data[filter_], 0.05, 'AB', filter_ , None, None, data['observatory'], True, data['file']])
        if save:
            formatted_tbl.write('timeseriesspec.synphot', format = 'ascii.fixed_width', overwrite = True)
            print('timeseriesspec.synphot is saved.')
        return formatted_tbl
    
    def get_spec_date(self,
                      time):
        closest_idx = np.argmin(np.abs(time-np.array(self._obsdate_data)))
        return self.spectrum[closest_idx], self._obsdate_data[closest_idx]

    def show_spec_date(self,
                       time,
                       show_flux_unit : str = 'flamb',
                       smooth_factor : int = 1,
                       normalize : bool = False,
                       normalize_cenwl : int = 6700,
                       normalize_wlwidth : int = 20,
                       label : str = '',
                       color : str = None
                       ):
        spec, obsdate = self.get_spec_date(time = time)
        spec.show(show_flux_unit= show_flux_unit, smooth_factor= smooth_factor, normalize= normalize, normalize_cenwl= normalize_cenwl, normalize_wlwidth= normalize_wlwidth, label = obsdate, color = color)


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
    speckey = '/data1/supernova_rawdata/SN2021aefx/spectroscopy/WISeREP/ascii/*'
    filelist = glob.glob(speckey)
    A = TimeSeriesSpectrum(specfilelist = filelist, flux_unit = 'flamb')
    A.photometry(save = True)
#%%
    A.show_spec_date(59533, normalize = True)
    A.show_spec_date(59535, normalize = True)
    A.show_spec_date(59536, normalize = True)

    plt.legend()

    #A.show_spec_date(59536, normalize = True)

# %%
