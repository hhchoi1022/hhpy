


#%%
from astropy.io import ascii
#%%
from timeseriesspectrum import TimeSeriesSpectrum
from HHsupport_analysis import read_Polin2019_spec, lflux_to_nuflux, fnu_to_mag
import glob, os
import re
import numpy as np
from astropy.table import Table
import pyphot
#%%


class HesmaSpectrum:
    
    def __init__(self,
                 specfilekey : str = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/spectrum/HESMA/*/*.dat'):
        self._specfilekey = specfilekey
        self._authorlist = self._get_authorlist()
        self.filelist = self._get_filelist()
        self.model_path = self._get_all_models()
        self.models = list(self.model_path.keys())
    
    def _get_filelist(self):
        return glob.glob(self._specfilekey)
    
    def _get_authorlist(self):
        return list(set([os.path.basename(os.path.dirname(file_)) for file_ in glob.glob(self._specfilekey)]))
    
    def _get_author_modellist(self,
                              author : str):
        model_key = os.path.dirname(os.path.dirname(self._specfilekey))+'/'+author+'/*.dat'
        model_path = glob.glob(model_key)
        model_name = [os.path.basename(path) for path in model_path]
        return model_path, model_name
    
    def _get_all_models(self):
        all_models = dict()
        for author in self._authorlist:
            all_models[author] = dict()
            model_path, model_name = self._get_author_modellist(author)
            all_models[author] = model_path
        return all_models
    
    def _get_filename(self,
                      wdmass : float,
                      he_shellmass : float):
        pattern = f'_{wdmass}_{he_shellmass}_'
        filename = np.array(self._filelist)[[pattern in file_ for file_ in self._filelist]]
        if len(filename) == 1:
            return filename[0]
        else:
            raise ValueError(f'WDmass:{wdmass}, He_shellmass{he_shellmass} not exist in the database')
    
    def timeSeriesSpectrum(self,
                           filename : str):
        data_tbl = ascii.read(filename, header_start = 0, data_start = 1 )
        wavelength = np.array(data_tbl['0.00'])
        data_tbl.remove_column('0.00')
        phase = np.array(list(data_tbl.columns)).astype(float)
        flux_array = data_tbl.to_pandas().to_numpy().transpose()
        wavelength = wavelength * np.ones(flux_array.shape)
        spec = TimeSeriesSpectrum(time_array = phase, wavelength_array= wavelength, flux_array = flux_array, flux_unit = 'flamb')
        spec.file = filename
        return spec

        
#%%
if __name__ == '__main__':
    import glob, os
    specfile_key = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/spectrum/HESMA/*/*.dat'
    A = HesmaSpectrum(specfile_key)

# %%
