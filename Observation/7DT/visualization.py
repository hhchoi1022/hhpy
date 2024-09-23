#%%
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import glob
#%%
class SDTVisualizer:
    def __init__(self, ascii_files, filters):
        self.ascii_files = ascii_files  # List of ASCII files
        self.filters = filters          # List of filter names (e.g., m400, m425, ..., m875)
        self.dataframes = dict()            # To store data from ASCII files
        
    def load_data(self):
        # Load each ASCII file into a dataframe and store it
        for file, filter_ in zip(self.ascii_files, self.filters):
            df = pd.read_csv(file, delim_whitespace=True)  # assuming space-delimited ASCII file
            self.dataframes[filter_] = df
    
    def find_nearest_source(self, ra, dec, radius_arcsec=2.0):
        # Match sources based on coordinates with a given radius
        matched_sources = []
        matched_filters = []
        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
        
        for filter_, df in self.dataframes.items():
            source_coords = SkyCoord(ra=df['ALPHA_J2000']*u.degree, dec=df['DELTA_J2000']*u.degree)
            sep = coord.separation(source_coords).arcsec
            idx_match = np.argmin(sep)
            if sep[idx_match] < radius_arcsec:
                matched_sources.append(df.iloc[idx_match])  # Add the first matched source
                matched_filters.append(filter_)
                
        return matched_filters, matched_sources
    
    def get_photometry(self, filters, sources):
        wavelengths = []
        magnitudes = []
        magnitude_errors = []
        
        for filter_, source in zip(filters, sources):
            mag_key = f'MAG_APER_2_{filter_}'
            magerr_key = f'MAGERR_APER_2_{filter_}'
            
            mag = source[mag_key] if mag_key in source else np.nan
            magerr = source[magerr_key] if magerr_key in source else np.nan
            
            # Convert filter name to wavelength (assuming 'm400' means 4000Å)
            wavelength = int(filter_[1:]) * 10
            
            wavelengths.append(wavelength)
            magnitudes.append(mag)
            magnitude_errors.append(magerr)
        
        return np.array(wavelengths), np.array(magnitudes), np.array(magnitude_errors)
    
    def plot_photometry(self, wavelengths, magnitudes, magnitude_errors):
        plt.errorbar(wavelengths, magnitudes, yerr=magnitude_errors, fmt='o', label='Photometry', color='blue')
        plt.gca().invert_yaxis()  # Magnitudes are inverted (brighter objects have lower magnitudes)
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Magnitude')
        plt.title('GECKO24l')
        plt.legend()
        plt.show()
#%%
# Example usage:
# List of ASCII files (replace with actual file paths)

ascii_files = glob.glob('/data1/7DT_obsdata/T04646/*.cat')
filters = [re.findall(r'_(m\d{3})_', ascii_files[i])[0] for i in range(len(ascii_files))]
# Initialize the class with the files and filters
photometry = SDTVisualizer(ascii_files, filters)

# Load data from the files
photometry.load_data()

# Coordinates of the target source (in degrees)
ra_target =  110.4148614  # Example RA
dec_target =  -49.2397919   # Example Dec
#ra_target =  109.723912  # Example RA
#dec_target =  -49.5148821   # Example Dec

# Matching radius in arcseconds
radius_arcsec = 5

# Find the source at the given coordinates
matched_filters, matched_sources = photometry.find_nearest_source(ra_target, dec_target, radius_arcsec)

# Extract the photometry for the matched source
wavelengths, magnitudes, magnitude_errors = photometry.get_photometry(matched_filters, matched_sources)

# Plot the photometric spectrum
photometry.plot_photometry(wavelengths, magnitudes, magnitude_errors)

# %%
