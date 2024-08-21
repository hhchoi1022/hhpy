

#%%
import mosfit
import numpy as np
import json
from astropy.table import Table

from Research.spectroscopy import TimeSeriesSpectrum
#%%

#%%
class mosfitData:
    
    def __init__(self,
                 target : str = 'SN2021aefx',
                 save : bool = True,
                 save_path : str = './',
                 photometry : bool = True,
                 spectra : bool = True):
        self.target = target
        self.data = self.find_target(target)
        self._save = save
        self._path = save_path
        if photometry:
            self.photdata = self._photdata()
        if spectra:
            self.specdata = self._specdata()
            
    
    def find_target(self,
                    target):
        fetcher = mosfit.fetcher.Fetcher()
        fetched = fetcher.fetch(target)[0]
        path = fetched['path']
        with open(path, 'r') as file_:
            data = json.load(file_) 
        self.data = data[target]
        return data[target]
    
    def _photdata(self):
        if not 'photometry' in list(self.data.keys()):
            raise MemoryError(f'No photometric data for {self.target}')
        header_key = ['obsdate', 'mag', 'e_mag', 'magsys', 'filter', 'depth_5sig', 'zp', 'observatory', 'detected']
        dtype_key = ['float64', 'float64', 'float64', 'str', 'str', 'float64', 'float64', 'str', 'bool']
        phot_tbl = Table(names = header_key, dtype = dtype_key)

        for photometry in self.data['photometry']:
            ul = 99
            e_ul = 99
            if 'time' in photometry.keys():
                obsdate = photometry['time']
            else:
                obsdate = 0

            if 'magnitude' in photometry.keys():
                mag = photometry['magnitude']
            else:
                mag = 99
                
            if 'e_magnitude' in photometry.keys():
                e_magnitude = photometry['e_magnitude']
            else:
                e_magnitude = 99
                                                
            if 'system' in photometry.keys():
                magsys = photometry['system']
            else:
                magsys = 'None'
                
            if 'band' in photometry.keys():
                filter_ = photometry['band']
            else:
                filter_ = 'None'
            
            if 'upperlimit' in photometry.keys():
                detected = False
                ul = mag
            else:
                detected = True
                ul = None
                                
            if 'telescope' in photometry.keys():
                telescope = photometry['telescope']
            else:
                telescope = 'None'

            if telescope == '':
                if 'instrument' in photometry.keys():
                    telescope = photometry['instrument']
                else:
                    telescope = 'None'
            
            if telescope == '':
                if 'observer' in photometry.keys():
                    telescope = photometry['observer']
                else:
                    telescope = 'None'
            
            phot_data = [obsdate, mag, e_magnitude, magsys, filter_, ul, None, telescope, detected]
            try:
                phot_tbl.add_row(phot_data)
            except:
                pass
        if self._save:
            os.makedirs(f'{self._path}mostfit/{self.target}', exist_ok = True)  
            phot_tbl.write(f'{self._path}mostfit/{self.target}/{self.target}.phot', format = 'ascii.fixed_width')
        return phot_tbl

    def _specdata(self):
        if not 'spectra' in list(self.data.keys()):
            raise MemoryError(f'No spectroscopic data for {self.target}')
        wavelengthlist = []
        flamblist = []
        mjdlist = []
        for spectrometry in self.data['spectra']:
            if not spectrometry['u_fluxes'] =='Uncalibrated':
                try:
                    filename = f'{spectrometry["filename"].split(".")[0]}_{spectrometry["time"]}.spec'
                except:
                    filename = f'{spectrometry["time"]}.spec'
                data = np.array(spectrometry['data']).astype(float)
                if 'u_errors' in spectrometry.keys():
                    wavelength = data[:,0]
                    flux = data[:,1]
                    fluxerr = data[:,2]
                    spec_tbl = Table()
                    spec_tbl['wavelength[AA]'] = wavelength
                    spec_tbl['flamb'] = flux
                    spec_tbl['e_flamb'] = fluxerr
                    spec_tbl['f_unit'] = spectrometry["u_fluxes"]
                    spec_tbl['lamb_unit'] = spectrometry["u_wavelengths"]
                else:
                    wavelength = data[:,0]
                    flux = data[:,1]
                    spec_tbl = Table()
                    spec_tbl['wavelength[AA]'] = wavelength
                    spec_tbl['flamb'] = flux
                    spec_tbl['f_unit'] = spectrometry["u_fluxes"]
                    spec_tbl['lamb_unit'] = spectrometry["u_wavelengths"]
                mjd = np.array(spectrometry['time']).astype(float)
                mjdlist.append(mjd)
                wavelengthlist.append(wavelength)
                flamblist.append(flux)
            if self._save:
                os.makedirs(f'{self._path}mostfit/{self.target}', exist_ok = True)  
                spec_tbl.write(f'{self._path}mostfit/{self.target}/{filename}', format = 'ascii.fixed_width')
        
        phase_array = np.array(mjdlist)
        wavelength_array = np.array(wavelengthlist, dtype = 'object')
        flamb_array = np.array(flamblist, dtype = 'object')
        
        timespec = TimeSeriesSpectrum(time_array = phase_array, wavelength_array= wavelength_array, flux_array= flamb_array, flux_unit= 'flamb')
        return timespec
                
                
                
            
                
                
                
        
            
            
        

                

    

# %%
if __name__ == '__main__':
    A = mosfitData(target = 'SN2017cbv', spectra = True, save = True)
# %%
