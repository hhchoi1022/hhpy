#%%


from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import median_smooth, extract_region

from astropy.io import ascii
import astropy.units as u

from pyphot import unit
import matplotlib.pyplot as plt
import numpy as np
from Research.helper import Helper
#%%
class Spectrum:
    
    def __init__(self,
                 wavelength : np.array,
                 flux : np.array,
                 flux_unit : str = 'fnu'):
        self.helper = Helper()
        self._wavelength_data = wavelength
        self._flux_data = flux
        if flux_unit.upper() in ['FNU', 'FLAMB', 'JY', 'MJY']:
            self.flux_unit = flux_unit
        else:
            raise ValueError('%s is not registered as spectral flux unit'%flux_unit)
        self.spectrum = self._set_spectrum()
        self.wavelength = self.spectrum.wavelength
        self.flux = self.spectrum.flux
        self.fnu = self._convert_fnu()
        self.flamb = self._convert_flamb()
        self.fjy = self._convert_jy()
        self.fmjy = self._convert_mjy()

    @property
    def data(self):
        spec_data = dict()
        spec_data['wavelength'] = self.wavelength
        spec_data['flux'] = self.flux
        return spec_data    

    def _set_spectrum(self):
        wavelength = self._wavelength_data * u.AA
        if self.flux_unit.upper() == 'FNU':
            f_unit = 'erg cm-2 s-1 /Hz'
        elif self.flux_unit.upper() == 'FLAMB':
            f_unit = 'erg cm-2 s-1 /AA'
        elif self.flux_unit.upper() == 'JY':
            f_unit = 'jy'
        elif self.flux_unit.upper() == 'MJY':
            f_unit = 'mJy'
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
    
    def _convert_mjy(self):
        f_unit = 'mJy'
        return self.spectrum.new_flux_unit(u.Unit(f_unit)).flux
    
    def show(self,
             show_flux_unit : str = 'fnu',
             smooth_factor : int = 1,
             redshift : float = 0,
             normalize : bool = False,
             normalize_cenwl : int = 6700,
             normalize_wlwidth : int = 20,
             linewidth : int = 1,
             linestyle : str = '-',
             label : str = '',
             color : str = 'k',
             log : bool = False,
             offset : float = 0,
             axis = None):
        if not show_flux_unit.upper() in ['FNU', 'FLAMB', 'JY', 'AB', 'MJY']:
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
        elif show_flux_unit.upper() == 'MJY':
            wavelength = self.wavelength
            wavelength_label = self.wavelength.unit.to_string()
            flux = self.fmjy
            flux_label = self.fmjy.unit.to_string()
        if show_flux_unit.upper() == 'AB':
            wavelength = self.wavelength
            wavelength_label = self.wavelength.unit.to_string()
            flux = self.fnu
            mag = -2.5*np.log10(flux.value) - 48.6
            flux = mag * u.AB
            flux_label = "AB"
            
        spec = Spectrum1D(spectral_axis = wavelength * (1+redshift), flux = flux)
        if smooth_factor != 1:
            spec = median_smooth(spec, smooth_factor)
        if normalize:
            norm_region = SpectralRegion((normalize_cenwl-normalize_wlwidth//2)*u.AA, (normalize_cenwl+normalize_wlwidth//2)*u.AA)
            norm_spec = extract_region(spec, norm_region)
            spec /= np.mean(norm_spec.flux)
            flux_label = 'Normalized flux density'
        if axis is None:
            if log:
                plt.step( np.log10(spec.wavelength.value), np.log10(spec.flux.value), label = label, c = color, linestyle = linestyle, linewidth = linewidth)
            else:
                plt.step(spec.wavelength, spec.flux + offset, label = label, c = color, linestyle = linestyle, linewidth = linewidth)
            plt.xlabel(wavelength_label)
            plt.ylabel(flux_label)
        else:
            if log:
                axis.step( np.log10(spec.wavelength.value), np.log10(spec.flux.value), label = label, c = color, linestyle = linestyle, linewidth = linewidth)
            else:
                axis.step(spec.wavelength, spec.flux + offset, label = label, c = color, linestyle = linestyle, linewidth = linewidth)
            axis.set_xlabel(wavelength_label)
            axis.set_ylabel(flux_label)
        
    def photometry(self,
                   filterset : str = 'UBVRIugri',
                   mag_type : str = 'AB',
                   visualize : bool = False,
                   visualize_unit : str = 'mag' # mag or flux
                   ):
        from astropy.table import Table
        import pyphot
        phot_tbl = dict()
        color_key, _, _, pyphot_key, _ = self.helper.load_filt_keys(filterset)
        lib = pyphot.get_library()

        if visualize:
            # Create a figure and axis for the spectrum
            fig, ax1 = plt.subplots(figsize = (10,6))
            ax1.set_xlabel('Wavelength (Å)')
            if visualize_unit == 'mag':
                self.show(show_flux_unit='AB', smooth_factor=11, log=False, color='k', label='Spectrum')
                ax1.set_ylabel('AB magnitude', color='k')
            else:
                self.show(show_flux_unit='flamb', smooth_factor=11, log=False, color='k', label='Spectrum')
                ax1.set_ylabel('Flux Density', color='k')   
            ax1.tick_params(axis='y', labelcolor='k')

            # Create a secondary axis for filter transmission
            ax2 = ax1.twinx()
            ax2.set_ylabel('Filter Transmission', color='k')
            ax2.tick_params(axis='y', labelcolor='k')

        for filter_ in filterset:
            filt_pyphot = lib[pyphot_key[filter_]]
            flux = filt_pyphot.get_flux(self.wavelength.value * unit['AA'], self.flamb.value * unit['ergs/s/cm**2/AA'], axis=1)
            if filter_ in 'UBVRI':
                mag = -2.5 * np.log10(flux.value) - filt_pyphot.Vega_zero_mag
            else:
                mag = -2.5 * np.log10(flux.value) - filt_pyphot.AB_zero_mag
            phot_tbl[filter_] = np.round(mag, 5)
            
            if visualize:
                # Check if the transmission is in percentage or fraction and normalize
                transmission = filt_pyphot.transmit
                if np.max(transmission) > 1:  # Assuming if max value is greater than 1, it is in percentage
                    transmission_normalized = transmission / 100.0
                else:  # It is already in fraction form
                    transmission_normalized = transmission
                    
                # Plot the photometric point on the spectrum axis
                if visualize_unit == 'mag':
                    ax1.scatter(filt_pyphot.lpivot.value, mag, label=f'{filter_} mag: {mag:.2f}', color=color_key[filter_])
                    ax1.invert_yaxis()
                else:
                    ax1.scatter(filt_pyphot.lpivot.value, flux.value, label=f'{filter_} flux: {flux.value:.2e}', color=color_key[filter_])

                # Plot the normalized filter transmission curve on the secondary axis
                ax2.plot(filt_pyphot.wavelength.value, transmission_normalized, alpha=0.5, color=color_key[filter_])
   
        if visualize:
            # Add legends and display the plot
            ax1.legend(loc='upper right', ncol = 1)
            plt.xlim(np.min(self.wavelength.value) -1000, np.max(self.wavelength.value) + 1000)
            plt.title('Spectrum with Filter Transmission Curves')
            plt.show()
    
        return phot_tbl
    
# %%
if __name__ == '__main__':
    from astropy.io import fits
    import glob
    from spectroscopyfile import SpectroscopyFile
    specfilelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/spectroscopy/WISeREP/ascii/*'))
    specfile = SpectroscopyFile(specfilelist[0])
    flux = specfile.flux
    wl = specfile.wavelength
    spec = Spectrum(wl, flux, flux_unit = 'flamb')
    #spec.show(show_flux_unit = 'flamb', normalize= False, smooth_factor = 11, log = False)
    synth_phot_tbl = spec.photometry(visualize = True, visualize_unit = 'flux')
    print('specfile.obsdate: ',specfile.obsdate)
    print('g-r: ',synth_phot_tbl['g'] - synth_phot_tbl['r'])
    print('B-V: ',synth_phot_tbl['B'] - synth_phot_tbl['V'])
# %%
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm  # Import the colormap

if __name__ == '__main__':
    plt.figure(dpi = 300)
    visualize_files =[specfilelist[-5],specfilelist[-6],specfilelist[-8],specfilelist[-10], specfilelist[-11],specfilelist[-12],specfilelist[-13]] #specfilelist[0,6]  # Choose the first 6 files to visualize
    visualize_files =[specfilelist[0],specfilelist[1], specfilelist[3], specfilelist[4]]
    num_files = len(visualize_files)  # Determine the number of files
    colormap = cm.get_cmap('jet', num_files)  # Choose a colormap and set the number of colors

    for i, file_ in enumerate(visualize_files):  # Loop over all files in specfilelist
        specfile = SpectroscopyFile(file_)
        flux = specfile.flux
        wl = specfile.wavelength
        spec = Spectrum(wl, flux, flux_unit='flamb')
        
        from Research.helper import AnalysisHelper
        Planck = AnalysisHelper().planck
        val = Planck(temperature=8000, wl_AA=wl)
        
        # Get a color from the colormap
        color = colormap(i)
        
        # If spec.show() supports a color parameter
        spec.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.04, normalize_cenwl=7500, color=color, label = rf'$t_0$ + {np.round(float(specfile.obsdate-59529.3318),2)}', offset = -2*i)
        plt.axvline(6300, linestyle = '--', c= 'r', linewidth = 0.5)
        plt.axvline(6364, linestyle = '--', c= 'r', linewidth = 0.5)
        plt.axvline(5875, linestyle = '--', c= 'b', linewidth = 0.5)
        plt.axvline(6678, linestyle = '--', c= 'b', linewidth = 0.5)
        plt.axvline(6563, linestyle = '--', c= 'b', linewidth = 0.5)
    plt.legend(ncol = 2)
    plt.legend()
    plt.show()
#%%
if __name__ == '__main__':
    file1 = specfilelist[0]
    file2 = specfilelist[1]
    normalize_cenwl = 7000
    normalize_wlwidth = 20
    bb_ratio = 30
    pw_ratio = 60
    
    specfile1 = SpectroscopyFile(file1)
    flux1 = specfile1.flux
    wl1 = specfile1.wavelength
    spec1 = Spectrum(wl1, flux1, flux_unit='flamb')
    color = colormap(0)
    spec1.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.0, normalize_cenwl=normalize_cenwl, color=color, label = rf'$t_0$ + {np.round(float(specfile1.obsdate-59529.3318),2)}')

    specfile2 = SpectroscopyFile(file2)
    flux2 = specfile2.flux
    wl2 = specfile2.wavelength
    spec2 = Spectrum(wl2, flux2, flux_unit='flamb')
    color = colormap(1)
    spec2.show(show_flux_unit='flamb', normalize=True, smooth_factor=11, log=False, redshift=0.0, normalize_cenwl=normalize_cenwl, color=color, label = rf'$t_0$ + {np.round(float(specfile2.obsdate-59529.3318),2)}')
    
    from Research.helper import AnalysisHelper
    Planck = AnalysisHelper().planck
    spec = Spectrum1D(spectral_axis = spec2.wavelength, flux = spec2.flamb)
    norm_region = SpectralRegion((normalize_cenwl-normalize_wlwidth//2)*u.AA, (normalize_cenwl+normalize_wlwidth//2)*u.AA)
    norm_spec = extract_region(spec, norm_region)
    spec /= np.mean(norm_spec.flux)
    bb = Planck(temperature= 20000, wl_AA=wl2)
    spec_bb = Spectrum1D(spectral_axis = bb['wl']*u.AA, flux = bb['flamb']*u.Unit('erg cm-2 s-1 /AA'))
    norm_region = SpectralRegion((normalize_cenwl-normalize_wlwidth//2)*u.AA, (normalize_cenwl+normalize_wlwidth//2)*u.AA)
    norm_spec_bb = extract_region(spec_bb, norm_region)
    spec_bb /= np.mean(norm_spec_bb.flux)
    flux3 = spec.flux.value * pw_ratio + spec_bb.flux.value * bb_ratio
    spec_comb = Spectrum1D(spectral_axis = bb['wl']*u.AA, flux = flux3*u.Unit('erg cm-2 s-1 /AA'))
    norm_region = SpectralRegion((normalize_cenwl-normalize_wlwidth//2)*u.AA, (normalize_cenwl+normalize_wlwidth//2)*u.AA)
    norm_spec_comb = extract_region(spec_comb, norm_region)
    spec_comb /= np.mean(norm_spec_comb.flux)
    plt.step(spec_comb.wavelength, spec_comb.flux)
    #spec3 = Spectrum(wl2, flux3, flux_unit='flamb')
    #spec3.show(show_flux_unit='flamb', normalize=False, smooth_factor=11, log=False, redshift=0.05, normalize_cenwl=7500, color='k', label = rf'$t_0$ + {np.round(float(specfile2.obsdate-59529.3318),2)} + BB')
    
    
    
    plt.legend()


# %%    
if __name__ == '__main__':  
    redshift = 0.005
    model = ascii.read('/home/hhchoi1022/syn++/files/es-0.98.1/example/SN2021aefx_211111_SALT.fit')
    specfile = SpectroscopyFile(file_)
    flux = specfile.flux
    wl = specfile.wavelength/(1+redshift)
    spec = Spectrum(wl, flux, flux_unit = 'flamb')
    #spec.show(show_flux_unit = 'flamb', normalize= True, smooth_factor = 11, log = False, normalize_cenwl=6700)
    Planck = AnalysisHelper().planck
    val  = Planck(temperature = 8000, wl_AA = wl)
    #plt.plot(val['wl'], + val['flamb']/np.mean(val['flamb'][np.where((np.arange(np.min(wl), np.max(wl), 1) > 6700-30) & (np.arange(np.min(wl), np.max(wl), 1) < 6700+30))]))   
    plt.axvline(3700)
    plt.axvline(5800)
    plt.axvline(8500)
    plt.axvline(7000)
    spec.spectrum
    norm_region = SpectralRegion((7000-30//2)*u.AA, (7000+30//2)*u.AA)
    norm_spec = extract_region(spec.spectrum, norm_region)
    spec.spectrum /= np.mean(norm_spec.flux)
    fit_flux = spec.spectrum.flux
    fit_wl = spec.spectrum.wavelength 
    fit_fluxerr = np.ones_like(fit_flux) * 0.1
    plt.plot(model['col1'], model['col2'])
    from astropy.table import Table
    tbl = Table()
    tbl['wl'] = np.round(fit_wl,2)
    tbl['flux'] = fit_flux
    tbl['fluxerr'] = fit_fluxerr
    spec_new = Spectrum(list(tbl['wl']), list(tbl['flux']), flux_unit = 'flamb')
    spec_new.show(show_flux_unit = 'flamb', normalize= False, smooth_factor = 11, log = False)
    tbl.write('SN2021aefx_211111_SALT.dat', format='ascii.basic', overwrite=True)
    tbl_cbv = ascii.read('/home/hhchoi1022/Desktop/Gitrepo/Research/spectroscopy/mostfit/SN2017cbv/57822.6863542.spec', format = 'fixed_width')
    wl_cbv, flux_cbv = tbl_cbv['wavelength[AA]'], tbl_cbv['flamb']
    #Spectrum(wl_cbv, flux_cbv, flux_unit = 'flamb').show(show_flux_unit = 'flamb', normalize= True, smooth_factor = 11, log = False,  normalize_cenwl=7500)
# %%
