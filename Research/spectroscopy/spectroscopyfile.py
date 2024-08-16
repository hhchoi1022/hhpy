#%%
import os
import datetime
from astropy.time import Time
class SpectroscopyFile:
    """
    A class to represent and read spectroscopy files in either FITS or ASCII format.

    Attributes
    ----------
    filename : str
        Path to the spectroscopy file.
    format : str
        Format of the file ('fits' or 'ascii').
    header : dict
        Header information from the file.
    wavelength : numpy.ndarray
        Array of wavelength values.
    flux : numpy.ndarray
        Array of flux values.

    Methods
    -------
    __init__(filename)
        Initializes the SpectroscopyFile object and determines the file format.
    __repr__()
        Returns a string representation of the SpectroscopyFile object.
    read_fits()
        Reads the FITS file format and extracts wavelength and flux data.
    read_ascii()
        Reads the ASCII file format and extracts wavelength and flux data.
    """

    def __init__(self, filename: str):
        """
        Initializes the SpectroscopyFile object and determines the file format.

        Parameters
        ----------
        filename : str
            Path to the spectroscopy file.
        
        Raises
        ------
        ValueError
            If the file format is not recognized.
        """
        self.filename = filename
        try:
            self.read_fits()
            self.format = 'fits'
        except:
            pass
        try:
            self.read_ascii()
            self.format = 'ascii'
        except:
            pass
        if not hasattr(self, 'format'):
            raise ValueError('File format not recognized')

    def __repr__(self):
        """
        Returns a string representation of the SpectroscopyFile object.

        Returns
        -------
        str
            A string representation of the SpectroscopyFile object.
        """
        return f"SpectroscopyFile(Filename = {os.path.basename(self.filename)}, Format = {self.format})"

    def read_fits(self):
        """
        Reads the FITS file format and extracts wavelength and flux data.

        This method sets the header, wavelength, and flux attributes for the FITS file.

        Raises
        ------
        Exception
            If there is an error in reading the FITS file.
        """
        from astropy.io import fits
        import numpy as np

        # Open the FITS file
        hdulist = fits.open(self.filename)
        
        # Access the flux data
        flux = hdulist[0].data

        # Extract the header
        header = hdulist[0].header

        # Extract WCS information from the header
        crval1 = header['CRVAL1']  # Reference value (starting wavelength)
        cdelt1 = header['CDELT1']  # Wavelength increment per pixel
        crpix1 = header['CRPIX1']  # Reference pixel (usually 1)

        # Calculate the wavelength array
        n_pixels = flux.size
        wavelength = crval1 + cdelt1 * (np.arange(n_pixels) + 1 - crpix1)

        # Close the FITS file
        hdulist.close()

        # Set the attributes
        self.header = header
        self.wavelength = wavelength
        self.flux = flux
        self.obsdate = self._get_obsdate_from_header()
        if not self.obsdate:
            self.obsdate = self._get_obsdate_from_filename(self.filename)

    def read_ascii(self):
        """
        Reads the ASCII file format and extracts wavelength and flux data.

        This method sets the header, wavelength, and flux attributes for the ASCII file.

        Raises
        ------
        Exception
            If there is an error in reading the ASCII file.
        """
        import numpy as np
        
        with open(self.filename, 'r') as file:
            lines = file.readlines()
            
        # Separate the header and data
        header_lines = []
        data_lines = []
        for line in lines:
            if line.startswith('#'):
                header_lines.append(line.strip())
            else:
                data_lines.append(line.strip())

        # Parse the header (if needed, for now we are just reading it)
        header = {}
        for line in header_lines:
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().lstrip('# ')
                value = value.split('/')[0].strip().replace("'", "").replace(" ", "")
                header[key] = value

        # Read the data
        data = np.loadtxt(data_lines)
        wavelength = data[:, 0]
        flux = data[:, 1]

        # Set the attributes
        self.header = header
        self.wavelength = wavelength
        self.flux = flux
        self.obsdate = self._get_obsdate_from_header()
        if not self.obsdate:
            self.obsdate = self._get_obsdate_from_filename(self.filename)
    
    def _get_obsdate_from_header(self):
        obsdate = None
        try:
            obsdate = Time(self.header['DATE-OBS']).mjd
        except:
            pass
        try:
            obsdate = Time(self.header['JD'], format = 'jd').mjd
        except:
            pass
        return obsdate
        
    def _get_obsdate_from_filename(self,
                                   filename : str,
                                   date_pattern : str = '(\d\d\d\d)-(\d\d)-(\d\d)_(\d\d)-(\d\d)',
                                   date_pattern_2 : str = '(245\d\d\d\d.\d\d\d)'):
        import re
        try:
            year, month, day, hour, minute = re.search(date_pattern, filename).groups()
        except:
            jd, = re.search(date_pattern_2, filename).groups()
            dt = Time(jd, format= 'jd').datetime
            year, month, day, hour, minute = dt.year, dt.month, dt.day, dt.hour, dt.minute
        dt = datetime.datetime(year = int(year), month = int(month), day= int(day), hour = int(hour), minute = int(minute))
        obsdate = Time(dt).mjd
        return obsdate