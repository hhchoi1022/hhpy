
#%%
from timeseriesspectrum import TimeSeriesSpectrum
from HHsupport_analysis import read_Polin2019_spec
import glob
import re
import numpy as np
from astropy.table import Table
#%%


class PolinSpectrum:
    
    def __init__(self,
                 specfilekey : str):
        self._specfilekey = specfilekey
        self._filelist = self._get_filelist()
        self.cases = self._get_all_case()
    
    def _get_filelist(self):
        return glob.glob(self._specfilekey)
    
    def _get_pattern(self,
                     filename : str,
                     pattern : str = '_(\d+.\d+)_(\d+.\d+)'):
        
        wdmass, hemass = re.search(pattern, filename).groups()
        return float(wdmass), float(hemass)
    
    def _get_filename(self,
                      wdmass : float,
                      he_shellmass : float):
        pattern = f'_{wdmass}_{he_shellmass}_'
        filename = np.array(self._filelist)[[pattern in file_ for file_ in self._filelist]]
        if len(filename) == 1:
            return filename[0]
        else:
            raise ValueError(f'WDmass:{wdmass}, He_shellmass{he_shellmass} not exist in the database')
    
    def _get_all_case(self):
        wdmasslist = []
        hemasslist = []
        for file_ in self._filelist:
            wdmass, hemass = self._get_pattern(file_)
            wdmasslist.append(wdmass)
            hemasslist.append(hemass)
        all_case_tbl = Table()
        all_case_tbl['wdmass'] = wdmasslist
        all_case_tbl['heshellmass'] = hemasslist
        all_case_tbl.sort(['wdmass', 'heshellmass'])
        return all_case_tbl

    def timeSeriesSpectrum(self,
                           wdmass : float,
                           he_shellmass : float):
        filename = self._get_filename(wdmass = wdmass, he_shellmass= he_shellmass)
        Lnu, Llamb, lamb, time, mu = read_Polin2019_spec(filename = filename)
        lamb = lamb * np.ones(Lnu.shape)
        timespec = TimeSeriesSpectrum(time_array= time, wavelength_array= lamb, flux_array= Lnu, flux_unit='fnu')
        return timespec
        
        

    
    
    
        
        
#%%
if __name__ == '__main__':
    specfolder =  '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/spectrum/Polin/ddet_Polin2019/*.h5'
    A = PolinSpectrum(specfolder)

# %%
