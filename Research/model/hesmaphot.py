#%%
from hesmaspectrum import HesmaSpectrum
import os, glob
from astropy.io import ascii
import pyphot
from HHsupport_analysis import load_filt_keys
from astropy.table import Table
from pyphot import unit
import numpy as np
from convert_AB_Vega import ABVegaMagnitude
from astropy.table import vstack
#%%
class HesmaPhot:
    
    def __init__(self,
                 photfilekey : str = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/lightcurve/HESMA/*/*.dat',
                 specfilekey : str = '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/spectrum/HESMA/*/*.dat'):
        self._photfilekey = photfilekey
        self._specfilekey = specfilekey
        self._spec = HesmaSpectrum(specfilekey)
        self._authorlist = self._get_all_authorlist()
        self.specfilelist = self._spec.filelist
        self.photfilelist = self._get_all_filelist()
        self.model_path = self._get_all_models()
        self.models = list(self.model_path.keys())
        self.vegaconverter = ABVegaMagnitude.convert_vega_to_AB
        
    
    def _get_all_filelist(self):
        return glob.glob(self._photfilekey)
    
    def _get_all_authorlist(self):
        return list(set([os.path.basename(os.path.dirname(file_)) for file_ in glob.glob(self._photfilekey)]))
    
    def _get_all_author_modellist(self,
                                  author : str):
        model_key = os.path.dirname(os.path.dirname(self._photfilekey))+'/'+author+'/*.dat'
        model_path = glob.glob(model_key)
        model_name = [os.path.basename(path) for path in model_path]
        return model_path, model_name
    
    def _get_all_models(self):
        all_models = dict()
        for author in self._authorlist:
            all_models[author] = dict()
            model_path, model_name = self._get_all_author_modellist(author)
            all_models[author] = model_path
        return all_models
    
    def _get_author(self,
                    filename : str):
        return os.path.basename(os.path.dirname(filename))
    
    def _get_specmodel(self,
                       filename : str):
        return os.path.basename(filename).split('spectra')[0]
    
    def _get_photmodel(self,
                       filename : str):
        return os.path.basename(filename).split('lightcurve')[0]
        
        
    def synthphot(self,
                  specfilename : str,
                  filterkey : str = 'UBVRIugriz',
                  early : bool = False):
        model_key = self._get_specmodel(specfilename)
        author_key = self._get_author(specfilename)
        earlyfile = glob.glob(f'{os.path.dirname(os.path.dirname(self._photfilekey))}/{author_key}/{model_key}*early*.dat')
        early_tbl = Table()
        if len(earlyfile) > 0:
            early_tbl = ascii.read(earlyfile[0], format = 'fixed_width')
            early_tbl.rename_column('t','phase')
            '''            
            for filter_ in early_tbl.keys()[1:]:
                mag = early_tbl[filter_]
                ABmag = self.vegaconverter(self, mag, filter_ = filter_)
                early_tbl[filter_] = ABmag
            '''
        else:
            pass

        spec = self._spec.timeSeriesSpectrum(specfilename)
        allphase = spec.phase 
        spec_data = spec.data
        synth_phot_tbl = Table()
        
        
        synth_phot_tbl['phase'] = allphase
        lib = pyphot.get_library()
        _, _, _, pyphot_key, _ = load_filt_keys(filterkey)
        for filter_ in filterkey:
            filt_pyphot = lib[pyphot_key[filter_]]
            magset = []
            for phase in allphase:
                spec = spec_data[phase]
                wl = spec['wavelength']
                flamb = spec['flux']
                flux = filt_pyphot.get_flux(wl*unit['AA'],flamb*unit['ergs/s/cm**2/AA'], axis = 1)
                mag = round(-2.5*np.log10(flux.value) - filt_pyphot.AB_zero_mag,4)
                magset.append(mag)
            synth_phot_tbl[f'{filter_}'] = magset
        if early:
            synth_phot_tbl = vstack([synth_phot_tbl, early_tbl])
        else:
            pass
        
        synth_phot_tbl.sort('phase')
        return synth_phot_tbl

        

# %%
A = HesmaPhot(specfilekey= '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/spectrum/HESMA/*/*.dat')
#%%

earlykey= '/Users/hhchoi1022/Gitrepo/data/SN2021aefx/model/lightcurve/Polin/*/*/original/*.mag'
early_files = glob.glob(earlykey)
from HHsupport_analysis import read_Polin2019
for early_file in early_files:
    '''    
    path =  os.path.dirname(early_file)+'/original/'
    os.makedirs(path, exist_ok= True)
    os.system(f'cp {early_file} {path}')
    '''
    path =  os.path.dirname(os.path.dirname(early_file))+'/'
    dat = read_Polin2019(early_file)
    for filter_ in dat.keys()[1:]:
        if filter_ in 'UBVRI':
            ABmag = ABVegaMagnitude.convert_vega_to_AB(ABVegaMagnitude, magnitude_vega = dat[filter_], filter_ = filter_)
        else:
            ABmag = dat[filter_]
        dat[filter_] = np.array(ABmag).round(3)
    dat.write(path + os.path.basename(early_file), overwrite = True, format= 'ascii')

# %%
