

#%%
import mosfit
import numpy as np
import json
from astropy.table import Table
from timeseriesspectrum import TimeSeriesSpectrum
#%%

#%%
class LoadData_mosfit:
    
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
        header_key = ['obsdate', 'mag', 'e_mag', 'magsys', 'filter', 'UL5_4', 'e_UL5_4', 'observatory', 'status', 'unit_time']
        dtype_key = ['float64', 'float64', 'float64', 'str', 'str', 'float64', 'float64', 'str', 'str', 'str']
        phot_tbl = Table(names = header_key, dtype = dtype_key)

        for photometry in self.data['photometry']:
            ul = 99
            e_ul = 99
            if 'time' in photometry.keys():
                obsdate = photometry['time']
            else:
                obsdate = 0
            
            if 'u_time' in photometry.keys():
                u_time = photometry['u_time']
            else:
                u_time = 'None'
                
            if 'band' in photometry.keys():
                filter_ = photometry['band']
            else:
                filter_ = 'None'
                
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

            if 'upperlimit' in photometry.keys():
                status = photometry['upperlimit']
                if status:
                    status = 'UL'
                    ul = mag
                    e_ul = e_magnitude
                else:
                    status = 'detected'

            else:
                status = 'detected'
                
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
            
            phot_data = [obsdate, mag, e_magnitude, magsys, filter_, ul, e_ul, telescope, status, u_time]
            try:
                phot_tbl.add_row(phot_data)
            except:
                pass
        if self._save:
            phot_tbl.write(f'{self._path}{self.target}_phot.dat', format = 'ascii.fixed_width')
        return phot_tbl

    def _specdata(self):
        if not 'spectra' in list(self.data.keys()):
            raise MemoryError(f'No spectroscopic data for {self.target}')
        wavelengthlist = []
        flamblist = []
        mjdlist = []
        for spectrum in self.data['spectra']:
            if not spectrum['u_fluxes'] =='Uncalibrated':
                filename = spectrum['filename']
                data = np.array(spectrum['data']).astype(float)
                if 'u_errors' in spectrum.keys():
                    wavelength = data[:,0]
                    flux = data[:,1]
                    fluxerr = data[:,2]
                    spec_tbl = Table()
                    spec_tbl['wavelength[AA]'] = wavelength
                    spec_tbl['flamb'] = flux
                    spec_tbl['e_flamb'] = fluxerr
                else:
                    wavelength = data[:,0]
                    flux = data[:,1]
                    spec_tbl = Table()
                    spec_tbl['wavelength[AA]'] = wavelength
                    spec_tbl['flamb'] = flux
                mjd = np.array(spectrum['time']).astype(float)
                mjdlist.append(mjd)
                wavelengthlist.append(wavelength)
                flamblist.append(flux)
                
                save_filename = f'[{mjd}]+{filename}+.dat'
            if self._save:
                spec_tbl.write(f'{self._path}{save_filename}', format = 'ascii.fixed_width')
        
        phase_array = np.array(mjdlist)
        wavelength_array = np.array(wavelengthlist, dtype = 'object')
        flamb_array = np.array(flamblist, dtype = 'object')
        
        timespec = TimeSeriesSpectrum(time_array = phase_array, wavelength_array= wavelength_array, flux_array= flamb_array, flux_unit= 'flamb')
        return timespec
                
                
                
            
                
                
                
        
            
            
        

                

    

# %%
if __name__ == '__main__':
    A = LoadData_mosfit(target = 'SN2011fe', spectra = False, save = True)
    #%%
    import matplotlib.pyplot as plt
    for phase in A.specdata.phase[1:5]:
        A.specdata.show_spec_date(phase, normalize = True, label = phase)
    plt.legend()
# %%
