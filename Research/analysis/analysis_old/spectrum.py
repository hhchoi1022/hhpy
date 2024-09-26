#%%


from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import median_smooth, extract_region

from astropy.io import ascii
import astropy.units as u

from pyphot import unit
import matplotlib.pyplot as plt
import numpy as np
#%%
class Spectrum:
    
    def __init__(self,
                 wavelength : np.array,
                 flux : np.array,
                 flux_unit : str = 'fnu',
                 phot_filterset = 'UBVRIugri'):
        self._wavelength_data = wavelength
        self._flux_data = flux
        if flux_unit.upper() in ['FNU', 'FLAMB', 'JY']:
            self.flux_unit = flux_unit
        else:
            raise ValueError('%s is not registered as spectral flux unit'%flux_unit)
        self.spectrum = self._set_spectrum()
        self.wavelength = self.spectrum.wavelength
        self.flux = self.spectrum.flux
        self.fnu = self._convert_fnu()
        self.flamb = self._convert_flamb()
        self.fjy = self._convert_jy()
        self.photometry = self.get_photometry(filterset = phot_filterset)
    
    def _set_spectrum(self):
        wavelength = self._wavelength_data * u.AA
        if self.flux_unit.upper() == 'FNU':
            f_unit = 'erg cm-2 s-1 /Hz'
        elif self.flux_unit.upper() == 'FLAMB':
            f_unit = 'erg cm-2 s-1 /AA'
        elif self.flux_unit.upper() == 'JY':
            f_unit = 'jy'
        flux = self._flux_data * u.Unit(f_unit)
        spec = Spectrum1D(spectral_axis = wavelength, flux = flux)
        return spec
    
    def _convert_fnu(self):
        f_unit = 'erg cm-2 s-1 /Hz'
        return self.spectrum.new_flux_unit(u.Unit(f_unit)).flux
    
    def _convert_flamb(self):
        f_unit = 'erg cm-2 s-1 /AA'
        return self.spectrum.new_flux_unit(u.Unit(f_unit)).flux
    
    def _convert_jy(self):
        f_unit = 'Jy'
        return self.spectrum.new_flux_unit(u.Unit(f_unit)).flux
    
    def show(self,
             show_flux_unit : str = 'fnu',
             smooth_factor : int = 1,
             normalize : bool = False,
             normalize_cenwl : int = 6700,
             normalize_wlwidth : int = 20,
             label : str = '',
             color : str = 'k'):
        if not show_flux_unit.upper() in ['FNU', 'FLAMB', 'JY']:
            raise ValueError('%s is not registered as spectral flux unit'%show_flux_unit)
        if show_flux_unit.upper() == 'FNU':
            wavelength = self.wavelength
            wavelength_label = self.wavelength.unit.to_string()
            flux = self.fnu
            flux_label = self.fnu.unit.to_string()
        elif show_flux_unit.upper() == 'FLAMB':
            wavelength = self.wavelength
            wavelength_label = self.wavelength.unit.to_string()
            flux = self.flamb
            flux_label = self.flamb.unit.to_string()
        elif show_flux_unit.upper() == 'JY':
            wavelength = self.wavelength
            wavelength_label = self.wavelength.unit.to_string()
            flux = self.fjy
            flux_label = self.fjy.unit.to_string()
            
        spec = Spectrum1D(spectral_axis = wavelength, flux = flux)
        if smooth_factor != 1:
            spec = median_smooth(spec, smooth_factor)
        if normalize:
            norm_region = SpectralRegion((normalize_cenwl-normalize_wlwidth//2)*u.AA, (normalize_cenwl+normalize_wlwidth//2)*u.AA)
            norm_spec = extract_region(spec, norm_region)
            spec /= np.mean(norm_spec.flux)
            flux_label = 'Normalized flux density'
        plt.step(spec.wavelength, spec.flux, label = label, c = color)
        plt.xlabel(wavelength_label)
        plt.ylabel(flux_label)
        
    def get_photometry(self,
                       filterset : str = 'UBVRIugri'):
        from HHsupport_analysis import load_filt_keys
        from astropy.table import Table
        import pyphot
        phot_tbl = dict()
        _, _, _, pyphot_key, _ = load_filt_keys(filterset)
        lib = pyphot.get_library()
        for filter_ in filterset:
            filt_pyphot = lib[pyphot_key[filter_]]
            flux = filt_pyphot.get_flux(self.wavelength.value * unit['AA'], self.flamb.value * unit['ergs/s/cm**2/AA'], axis = 1)
            mag = -2.5*np.log10(flux.value) - filt_pyphot.AB_zero_mag
            phot_tbl[filter_] = np.round(mag,5)
        return phot_tbl

        
        
        
    
    
    
        
    
    
        
# %%
if __name__ == '__main__':
    import glob
    specfilelist = glob.glob('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/spectrum/WISeREP/ascii/*.txt')
    ascii.read(specfilelist[0])
    from wiserepSpectrum import WiseRepSpectrum
    spec = WiseRepSpectrum('/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/spectrum/WISeREP/ascii/*.txt')
    #%%
    data_1 = ascii.read(glob.glob('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/spectrum/ascii/*.txt')[1])
    data_2 = ascii.read(glob.glob('/Users/hhchoi1022/Gitrepo/data/SN2021aefx/observation/spectrum/ascii/*.txt')[20])
    wl = data_1['col1']
    flux = data_1['col2']
    wl_2 = data_2['col1']
    flux_2 = data_2['col2']
    spec_1 = Spectrum(wl, flux, flux_unit = 'flamb')
    spec_2 = Spectrum(wl_2, flux_2, flux_unit = 'flamb')
    spec_1.show(show_flux_unit = 'jy', normalize= False, smooth_factor = 11)
    spec_2.show(show_flux_unit = 'jy', normalize= False, smooth_factor = 1)
    #%%