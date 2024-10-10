
# %%
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table
from astropy.wcs import WCS

from tqdm import tqdm
import os
import numpy as np
from astropy.table import unique
import inspect
import subprocess
import re

# %%

import signal
import functools

class TimeoutError(Exception):
    pass

class ActionFailedError(Exception):
    pass

def timeout(seconds=10, error_message="Function call timed out"):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Set the timeout signal
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                # Disable the alarm after the function completes
                signal.alarm(0)
            return result
        return wrapper
    return decorator

class PhotometryHelper():

    @property
    def photpath(self):
        # Get the file where the class is defined
        file_path = inspect.getfile(PhotometryHelper)

        # Convert the file path to an absolute path using os.path.abspath
        absolute_path = os.path.abspath(file_path)

        path_dir = os.path.join(os.path.dirname(absolute_path),'../photometry')

        return path_dir

    @property
    def scamppath(self):
        return os.path.join(self.photpath, 'scamp')
    
    @property  
    def swarppath(self):
        return os.path.join(self.photpath, 'swarp')
    
    @property
    def sexpath(self):
        return os.path.join(self.photpath, 'sextractor')
 
    def __repr__(self):
        methods = [f'PhotometryHelper.{name}()\n' for name, method in inspect.getmembers(
            PhotometryHelper, predicate=inspect.isfunction) if not name.startswith('_')]
        txt = '[Methods]\n'+''.join(methods)
        return txt
    # Load information

    def get_imginfo(self, filelist, keywords=['jd', 'group', 'filter', 'exptime', 'object', 'telescop', 'instrume', 'ra', 'dec', 'xbinning', 'ybinning', 'imagetyp']):
        '''
        parameters
        ----------
        1. filelist: (list) filelist of the fits images

        2. keywords : (list) list of the kwargs in the header

        returns 
        -------
        1. result_tbl : astropy.table
                        all information in the header of the filelist

        notes 
        -----
        "comment" and "history" are removed in the result_tbl due to too long length of the information 
        -----
        '''
        from ccdproc import ImageFileCollection
        from astropy.table import vstack
        from astropy.table import join 

        filedirs = set([os.path.dirname(filename) for filename in filelist])
        match_tbl = Table([filelist],names = ['file'])
        result_tbl = Table()
        for filedir in filedirs:
            files = [os.path.basename(file) for file in filelist if os.path.dirname(file)==filedir]
            iminfo = ImageFileCollection(location = filedir, keywords=keywords, filenames = files).summary
            absfiles = []
            for file in iminfo['file']:
                absfiles.append(filedir + '/'+ file)
            iminfo['file'] = absfiles
            result_tbl = vstack([iminfo, result_tbl])
        result_tbl = join(match_tbl, result_tbl, keys = 'file', join_type = 'inner')
        return result_tbl

    def get_telinfo(self, telescope: str = None, ccd: str = None, readoutmode: str = None, key_observatory='obs', key_ccd='ccd', key_mode = 'mode', obsinfo_file=None):
        '''
        parameters
        ----------
        1. observatory : str
                        observatory name to search in observaroy information file
        2. ccd : optional, str
                        ccd name of the observatroy (None)
        3. readoutmode : optional, str
                        readoutmode for RASA36 ccd [High, Merge] (High) 
        4. key_observatory : str
                        the key for searching observatory in observaroy information file (obs)
        5. key_ccd : str
                        the key for seraching ccd in observaroy information file (ccd)
        6. obsinfo_file : str
                        observaroy information file
        returns 
        -------
        1. obsinfo : astropy.table
                        The information of the observatory/ccd 

        notes 
        -----
        -----
        '''

        if obsinfo_file == None:
            obsinfo_file = os.path.join(self.photpath, 'CCD.dat')

        def output_valid_func(output): return len(output) == 1

        all_obsinfo = ascii.read(obsinfo_file, format='fixed_width')
        if telescope == None:
            telescope = input(
                f'Choose the Telescope : {set(all_obsinfo[key_observatory])}')
        obs_info = all_obsinfo[all_obsinfo[key_observatory] == telescope]
        if not telescope in all_obsinfo[key_observatory]:
            raise AttributeError(
                f'{telescope} information not exist.\n available :{set(all_obsinfo[key_observatory])}')

        if ccd != None:
            obs_info = obs_info[obs_info[key_ccd] == ccd]
            if output_valid_func(obs_info):
                return obs_info[0]
        
        if readoutmode != None:
            if readoutmode.upper() == 'MERGE':
                obs_info = obs_info[obs_info['mode'] == 'MERGE']
            if readoutmode.upper() == 'HIGH':
                obs_info = obs_info[obs_info['mode'] == 'HIGH']
            if readoutmode.upper() == 'LOW':
                obs_info = obs_info[obs_info['mode'] == 'HIGH']
            if output_valid_func(obs_info):
                return obs_info[0]
        
        # Other observatories (Multiple CCDs)
        else:
            if len(set(obs_info[key_ccd])) > 1:
                ccd = input(
                    f'Multiple CCDs of the observatory is found. Choose the CCD : {set(obs_info[key_ccd])}')
                obs_info = all_obsinfo[(all_obsinfo[key_observatory] == telescope) & (all_obsinfo[key_ccd] == ccd)]
                
            if len(set(obs_info[key_mode])) > 1:
                mode = input(
                    f'Multiple modes of the observatory is found. Choose the readout mode : {set(obs_info[key_mode])}')
                obs_info = all_obsinfo[(all_obsinfo[key_observatory] == telescope) & (all_obsinfo[key_mode] == mode)]
            
            if output_valid_func(obs_info):
                return obs_info[0]

        if not output_valid_func(obs_info):
            raise AttributeError(
                f'{ccd} information not exist in {telescope}.\n available CCD name:{list(set(all_obsinfo[all_obsinfo[key_observatory] == telescope][key_ccd]))}')

    def get_sexconfig(self):
        '''
        parameters
        ----------

        returns 
        -------
        default_config : dict
                        default configuration for Source Extractor
        notes 
        -------
        -----
        '''
        # source extractor configuration parameters
        conf_param = dict(
            # Default configuration file for Source Extractor 2.25.0
            # EB 2018-02-08
            #

            # -------------------------------- Catalog ------------------------------------

            CATALOG_NAME='./sextractor/temp.cat',       # name of the output catalog
            CATALOG_TYPE='ASCII_HEAD',     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
            # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
            PARAMETERS_NAME='default.param',  # name of the file containing catalog contents

            # ------------------------------- Extraction ----------------------------------

            # CCD (linear) or PHOTO (with gamma correction)
            DETECT_TYPE='CCD',
            DETECT_MINAREA=5,              # min. # of pixels above threshold
            # max. # of pixels above threshold (0=unlimited)
            DETECT_MAXAREA=0,
            # threshold type: RELATIVE (in sigmas)
            THRESH_TYPE='RELATIVE',
            # or ABSOLUTE (in ADUs)
            DETECT_THRESH=3,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
            ANALYSIS_THRESH=3,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

            FILTER='Y',              # apply filter for detection (Y or N)?
            FILTER_NAME='default.conv',   # name of the file containing the filter

            DEBLEND_NTHRESH=64,             # Number of deblending sub-thresholds
            DEBLEND_MINCONT=0.0001,          # Minimum contrast parameter for deblending

            CLEAN='N',              # Clean spurious detections? (Y or N)?
            CLEAN_PARAM=1.0,            # Cleaning efficiency

            MASK_TYPE='CORRECT',        # type of detection MASKing: can be one of
            # NONE, BLANK or CORRECT

            # -------------------------------- WEIGHTing ----------------------------------

            WEIGHT_TYPE='NONE',           # type of WEIGHTing: NONE, BACKGROUND,
            # MAP_RMS, MAP_VAR or MAP_WEIGHT
            # Rescale input weights/variances (Y/N)?
            RESCALE_WEIGHTS='N',
            WEIGHT_IMAGE='weight.fits',    # weight-map filename
            # modulate gain (E/ADU) with weights? (Y/N)
            WEIGHT_GAIN='N',

            # ------------------------------ Photometry -----------------------------------

            # MAG_APER aperture diameter(s) in pixels
            PHOT_APERTURES='6',
            PHOT_AUTOPARAMS='2.5,3.5',       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
            PHOT_PETROPARAMS='2.0,3.5',       # MAG_PETRO parameters: <Petrosian_fact>,
            # <min_radius>
            PHOT_AUTOAPERS='0.0,0.0',        # <estimation>,<measurement> minimum apertures
            # for MAG_AUTO and MAG_PETRO
            # flux fraction[s] used for FLUX_RADIUS
            PHOT_FLUXFRAC=0.5,

            # level (in ADUs) at which arises saturation
            SATUR_LEVEL=50000.0,
            # keyword for saturation level (in ADUs)
            SATUR_KEY='SATURATE',

            MAG_ZEROPOINT=0.0,            # magnitude zero-point
            # gamma of emulsion (for photographic scans)
            MAG_GAMMA=4.0,
            GAIN=0.0,            # detector gain in e-/ADU
            GAIN_KEY='GAIN',           # keyword for detector gain in e-/ADU
            # size of pixel in arcsec (0=use FITS WCS info)
            PIXEL_SCALE=1.0,

            # ------------------------- Star/Galaxy Separation ----------------------------

            SEEING_FWHM=5,            # stellar FWHM in arcsec
            STARNNW_NAME='default.nnw',    # Neural-Network_Weight table filename

            # ------------------------------ Background -----------------------------------

            BACK_TYPE='AUTO',           # AUTO or MANUAL
            BACK_VALUE=0.0,            # Default background value in MANUAL mode
            BACK_SIZE=64,             # Background mesh: <size> or <width>,<height>
            BACK_FILTERSIZE=3,              # Background filter: <size> or <width>,<height>

            BACKPHOTO_TYPE='GLOBAL',         # can be GLOBAL or LOCAL
            BACKPHOTO_THICK=24,             # thickness of the background LOCAL annulus
            BACK_FILTTHRESH=0.0,            # Threshold above which the background-
            # map filter operates

            # ------------------------------ Check Image ----------------------------------

            # can be NONE, BACKGROUND, BACKGROUND_RMS,
            CHECKIMAGE_TYPE='SEGMENTATION',
            # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
            # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
            # or APERTURES
            CHECKIMAGE_NAME='check.fits',     # Filename for the check-image

            # --------------------- Memory (change with caution!) -------------------------

            MEMORY_OBJSTACK=3000,           # number of objects in stack
            MEMORY_PIXSTACK=300000,         # number of pixels in stack
            MEMORY_BUFSIZE=1024,           # number of lines in buffer

            # --------------------------- Experimental Stuff -----------------------------

            PSF_NAME='default.psf',    # File containing the PSF model
            PSF_NMAX=1,              # Max.number of PSFs fitted simultaneously
            PATTERN_TYPE='RINGS-HARMONIC',  # can RINGS-QUADPOLE, RINGS-OCTOPOLE,
            # RINGS-HARMONICS or GAUSS-LAGUERRE
            SOM_NAME='default.som'    # File containing Self-Organizing Map weights
        )
        return conf_param

    def load_sexconfig(self, sexconfig: str) -> dict:
        config_dict = {}

        with open(sexconfig, 'r') as file:
            for line in file:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                # Split the line into key and value
                key_value = line.split(maxsplit=1)
                if len(key_value) == 2:
                    key, value = key_value
                    # Remove inline comments
                    value = value.split('#', 1)[0].strip()
                    # Attempt to convert value to appropriate type
                    try:
                        # Handle lists
                        if ',' in value:
                            value = [float(v) if '.' in v else int(v)
                                     for v in value.split(',')]
                        else:
                            # Convert to float if possible
                            value = float(
                                value) if '.' in value else int(value)
                    except ValueError:
                        # Keep as string if conversion fails
                        pass
                    config_dict[key] = value
        return config_dict

    # Calculation

    def to_skycoord(self, ra, dec, frame: str = 'icrs'):
        import astropy.units as u
        '''
		parameters
		----------
		1. ra : str or float
				Right ascension in diverse format(see notes)
		2. dec : str or float
				Declination in diverse format(see notes)
		
		returns 
		-------
		1. skycoord : SkyCoord
		
		notes 
		-----
		Current supported formats
				1. 15h32m10s, 50d15m01s
				2. 15 32 10, 50 15 01
				3. 15:32:10, 50:15:01
				4. 230.8875, 50.5369
		-----
		'''
        ra = str(ra)
        dec = str(dec)
        if (':' in ra) & (':' in dec):
            skycoord = SkyCoord(ra=ra, dec=dec, unit=(
                u.hourangle, u.deg), frame=frame)
        elif ('h' in ra) & ('d' in dec):
            skycoord = SkyCoord(ra=ra, dec=dec, unit=(
                u.hourangle, u.deg), frame=frame)
        elif (' ' in ra) & (' ' in dec):
            skycoord = SkyCoord(ra=ra, dec=dec, unit=(
                u.hourangle, u.deg), frame=frame)
        else:
            skycoord = SkyCoord(ra=ra, dec=dec, unit=(
                u.deg, u.deg), frame=frame)
        return skycoord

    def bn_median(self, masked_array, axis=None):
        """

        parameters
        ----------
        masked_array : `numpy.ma.masked_array`
                        Array of which to find the median.
        axis : optional, int 
                        Axis along which to perform the median. Default is to find the median of
                        the flattened array.

        returns
        ----------

        notes
        ----------
        Source code from Gregory S.H. Paek
        Perform fast median on masked array
        ----------
        """

        import numpy as np
        import bottleneck as bn
        data = masked_array.filled(fill_value=np.NaN)
        med = bn.nanmedian(data, axis=axis)
        # construct a masked array result, setting the mask from any NaN entries
        return np.ma.array(med, mask=np.isnan(med))

    # Table operation

    def cross_match(self, obj_catalog, sky_catalog, max_distance_second=5):
        '''
        parameters
        ----------
        1. obj_catalog : SkyCoord
                        source coordinates to be matched 
        2. sky_catalog : SkyCoord
                        coordinates of reference stars
        3. max_distance_second : float
                        separation distance in seconds which indicates tolerance for matching (5)

        returns 
        -------
        1. matched_object_idx : np.array
                        matched object catalog indices 
        2. matched_catalog_idx : np.array
                        matched sky catalog indices
        3. no_matched_object_idx : np.array
                        not matched object catalog indices

        notes 
        -----
        To extract matched objects from each catalog, 
        USE : obj_catalog[matched_object_idx] & sky_catalog[matched_catalog_idx]
        -----
        '''
        from astropy.coordinates import match_coordinates_sky
        closest_ids, closest_dists, closest_dist3d = match_coordinates_sky(obj_catalog, sky_catalog)
        max_distance = max_distance_second/3600
        matched_object_idx = []
        matched_catalog_idx = []
        no_matched_object_idx = []
        for i in range(len(closest_dists)):
            if closest_dists.value[i] < max_distance:
                matched_object_idx.append(i)
                matched_catalog_idx.append(closest_ids[i])
            else:
                no_matched_object_idx.append(i)
        return matched_object_idx, matched_catalog_idx, no_matched_object_idx

    def group_table(self, tbl: Table, key: str, tolerance: float = 0.1):
        '''
        parameters
        ----------
        group components in the {table} with the difference of the {key} smaller than the {tolerance}

        returns 
        -------
        1. groupped tables

        notes 
        -----

        -----
        '''
        from astropy.table import Table
        from astropy.table import vstack
        i = 0
        table = tbl.copy()
        table['group'] = 0
        groupped_tbl = Table()
        while len(table) >= 1:
            group_idx = (np.abs(table[0][key] - table[key]) < tolerance)
            group_tbl = table[group_idx]
            group_tbl['group'] = i
            remove_idx = np.where(group_idx == True)
            table.remove_rows(remove_idx)
            groupped_tbl = vstack([group_tbl, groupped_tbl])
            i += 1

        return groupped_tbl

    def match_table(self, tbl1, tbl2, key, tolerance=0.01):
        '''
        parameters
        ----------
        {two tables} to combine with the difference of the {key} smaller than the {tolerance}

        returns 
        -------
        1. combined table
        2. phase

        notes 
        -----
        Combined table have both columns of original tables. 
        They are horizontally combined in the order of tbl1, tbl2
        -----
        '''

        from astropy.table import vstack, hstack

        matched_tbl = Table()
        for obs in tbl1:
            ol_idx = (np.abs(obs[key] - tbl2[key]) < tolerance)
            if True in ol_idx:
                closest_idx = np.argmin(np.abs(obs[key]-tbl2[key]))
                compare_tbl = tbl2[closest_idx]
                # join(obs, compare_tbl, keys = 'observatory', join_type = 'outer')
                compare_tbl = hstack([obs, compare_tbl])
                matched_tbl = vstack([matched_tbl, compare_tbl])

        return matched_tbl

    def binning_table(self, tbl, key, tolerance=0.01):
        '''
        Parameters
        ----------
        tbl : Astropy.Table
                The input table to be binned.
        key : str
                The column name to apply the binning on.
        tolerance : float, optional
                The tolerance within which to bin the rows. Default is 0.01.

        Returns
        -------
        Astropy.Table
                The binned table with duplicates removed based on the specified tolerance.
        '''
        import pandas as pd
        table = tbl.to_pandas()
        # Sort the table by the key for efficient processing
        table = table.sort_values(by=key).reset_index(drop=True)

        binned_rows = []
        start_idx = 0

        while start_idx < len(table):
            end_idx = start_idx
            while (end_idx < len(table)) and (table[key].iloc[end_idx] - table[key].iloc[start_idx] < tolerance):
                end_idx += 1

            compare_table = table.iloc[start_idx:end_idx]

            # Aggregate the values within the tolerance range
            row = []
            for col in table.columns:
                if pd.api.types.is_numeric_dtype(table[col]):
                    result_val = round(compare_table[col].mean(), 4)
                else:
                    result_val = compare_table[col].iloc[0]
                row.append(result_val)

            binned_rows.append(row)
            start_idx = end_idx

        binned_table = pd.DataFrame(binned_rows, columns=tbl.columns)
        binned_tbl = Table().from_pandas(binned_table)
        return binned_tbl

    def remove_rows_table(self, tbl, column_key, remove_keys):
        '''
        Parameters
        ----------
        tbl : astropy.table.Table
                The input table from which rows need to be removed.
        column_key : str
                The column name based on which rows will be removed.
        remove_keys : str or list
                The value or list of values to be removed from the specified column in the table.

        Returns
        -------
        astropy.table.Table
                The table with specified rows removed.

        Notes
        -----
        This function removes rows from the input table where the values in the specified column
        match any of the values in `remove_keys`. `remove_keys` can be a single value (string)
        or a list of values. The function modifies the table in place and returns the modified table.
        -----
        '''
        if isinstance(remove_keys, str):
            remove_mask = tbl[column_key] == remove_keys
            remove_idx = np.where(remove_mask == True)
            tbl.remove_rows(remove_idx)
        else:
            for remove_key in remove_keys:
                remove_mask = tbl[column_key] == remove_key
                remove_idx = np.where(remove_mask == True)
                tbl.remove_rows(remove_idx)
        return tbl

 # Image processing
 
    def calculate_rotang(self, target_img, update_header : bool = False):
        from astropy.io import fits
        from astropy.wcs import WCS
        import numpy as np
        #fits_file = filelist[]
        # Load the FITS file with the astrometry solution
        #fits_file = im.target_image  # Replace with your actual solved FITS file
        hdul = fits.open(target_img)
        wcs = WCS(hdul[0].header)

        # Define a pixel at the center of the image
        ny, nx = hdul[0].data.shape
        center_pixel = [nx // 2, ny // 2]

        # Get the pixel coordinates offset along the y-axis (to simulate up direction)
        north_pixel = [center_pixel[0], center_pixel[1] + 1]

        # Convert these pixel positions to celestial coordinates (RA, Dec)
        center_coord = wcs.pixel_to_world(*center_pixel)
        north_coord = wcs.pixel_to_world(*north_pixel)

        # Calculate the position angle (PA) between the north pixel and the center pixel
        delta_ra = np.deg2rad(north_coord.ra.deg - center_coord.ra.deg) * np.cos(np.deg2rad(center_coord.dec.deg))
        delta_dec = np.deg2rad(north_coord.dec.deg - center_coord.dec.deg)

        pa_radians = np.arctan2(delta_ra, delta_dec)
        pa_degrees = np.rad2deg(pa_radians)

        # Adjust the angle to get the position angle from north
        if pa_degrees < 0:
            pa_degrees += 360
        
        if update_header:
            hdul[0].header['ROTANG'] = pa_degrees
            hdul.writeto(target_img, overwrite=True)

        print(f"Camera rotation angle (Position Angle) toward North: {pa_degrees:.2f} degrees")
        

    def cutout_img(self, target_img, size=0.9, prefix='cut_', xcenter=None, ycenter=None):
        '''
        parameters
        ----------
        1. target_img : str
                        {abspath} of the target image
        2. size : float or int
                        (float) ratio of the cut image (0.9)
                        (int) size of the cut image in pixel 
        3. prefix : str
                        prefix of the output image
        4. xcenter : optional, int or (float or str)
                        (int) pixel coordinate of the center
                        (float or str) RA coordinate of the center
        4. ycenter : optional, int or (float or str)
                        (int) pixel coordinate of the center
                        (float or str) Dec coordinate of the center
        returns 
        -------
        1. outputname : str
                        {abspath} of the cutout image

        notes 
        -----
        if xcenter, ycenter == None, cutout image will be centered to the center of the original image
        -----
        '''
        from astropy.wcs import WCS
        from astropy.nddata import Cutout2D

        hdu = fits.open(target_img)[0]
        wcs = WCS(hdu.header)
        if size < 1:
            size = size*int(len(hdu.data))
        if (xcenter == None) & (ycenter == None):
            xcenter, ycenter = len(hdu.data)//2, len(hdu.data)//2
        if (type(xcenter) != int) & (type(ycenter) != int):
            center_coords = self.to_skycoord(xcenter, ycenter)
            cutouted = Cutout2D(
                data=hdu.data, position=center_coords, size=size, wcs=wcs)
        else:
            cutouted = Cutout2D(data=hdu.data, position=(
                xcenter, ycenter), size=size, wcs=wcs)
        cutouted_hdu = hdu
        cutouted_hdu.data = cutouted.data
        cutouted_hdu.header.update(cutouted.wcs.to_header())
        cutouted_hdu.header['NAXIS1'] = int(size)
        cutouted_hdu.header['NAXIS2'] = int(size)
        outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
        cutouted_hdu.writeto(outputname, overwrite=True)
        return outputname

    def align_img(self, target_img, reference_img, prefix='align_', detection_sigma = 5):
        """

        parameters
        ----------
        1. target_img : str
                        absolute path of a target image 
        2. reference_img : str
                        absolute path of a reference image
                        The image as a reference frame. All images will be aligned on the basis of this image
        3. prefix = str
                        prefix of the aligned iamge

        returns
        ----------
        1. outputname : str
                        absolute path of the aligned image
        notes
        ----------
        ----------
        """

        import astroalign as aa
        from ccdproc import CCDData
        from astropy.wcs import WCS

        tgt_data, tgt_hdr = fits.getdata(target_img, header=True)
        ref_data, ref_hdr = fits.getdata(reference_img, header=True)
        
        ref_wcs = WCS(ref_hdr)
        wcs_hdr = ref_wcs.to_header(relax = True)
        wcs_hdr.remove('DATE-OBS', ignore_missing = True)
        wcs_hdr.remove('MJD-OBS', ignore_missing = True)
        wcs_hdr.remove('RADESYS', ignore_missing = True)
        wcs_hdr.remove('EQUINOX', ignore_missing = True)
        tgt_hdr.update(wcs_hdr)
        tgt_data = tgt_data.byteswap().newbyteorder()
        ref_data = ref_data.byteswap().newbyteorder()
        try:
            aligned_data, footprint = aa.register(tgt_data, ref_data, fill_value=0, detection_sigma= detection_sigma, max_control_points=30)
            aligned_tgt = CCDData(aligned_data, header=tgt_hdr, unit='adu')
            outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
            fits.writeto(outputname, aligned_tgt.data, aligned_tgt.header, overwrite=True)
            return outputname
        except:
            raise ActionFailedError('Failed to align the image. Check the image quality and the detection_sigma value.')

    def combine_img(self,
                    filelist,
                    clip='extrema',
                    combine='median',
                    scale='multiply',
                    prefix='com_',
                    zp_key='ZP5_1',

                    # Clipping
                    clip_sigma_low=2,
                    clip_sigma_high=5,
                    clip_minmax_min=3,
                    clip_minmax_max=3,
                    clip_extrema_nlow=1,
                    clip_extrema_nhigh=1,
                    ):
        '''
        parameters
        ----------
        1. filelist : list or np.array 
                        filelist to be combined
        2. clip : str
                        method for clipping [None, minmax, sigma, extrema] (sigma)
        3. combine : str
                        method for combining [mean, median, sum] (median)
        4. scale : str
                        method for scaling [None, zero, multiply] (zero)
        5. prefix : str
                        prefix of the combined image

        2.1. clip_sigma_low : optional, int
                        Threshold for rejecting pixels that deviate below the baseline value.
        2.2. clip_sigma_high : optional, int
                        Threshold for rejecting pixels that deviate above the baseline value.    
        2.3. clip_minmax_min : optional, int
                        If not None, all pixels with values below min_clip will be masked.
        2.4. clip_minmax_max : optional, int
                        If not None, all pixels with values above min_clip will be masked.
        2.5. clip_extrema_nlow : optional, int
                        If not None, the number of low values to reject from the combination.
        2.6. clip_extrema_nhigh : optional, int
                        If not None, the number of high values to reject from the combination.

        returns 
        -------
        1. outputname : str
                        absolute path of the combined image

        notes 
        -----
        For more information : https://ccdproc.readthedocs.io/en/latest/image_combination.html
        -----
        '''
        from ccdproc import CCDData
        from ccdproc import Combiner

        if len(filelist) <3:
            clip = None
            print('Number of filelist is lower than the minimum. Skip clipping process... \n')
            
        print('Combining... \n')
        print(60*'=')

        ccdlist = []
        hdrlist = []
        hdr = fits.getheader(filelist[0])
        for file in filelist:
            ccd = CCDData.read(file, unit = 'adu')
            hdr = fits.getheader(file)
            ccdlist.append(ccd)
            hdrlist.append(hdr)
        for i, file in enumerate(filelist):
            hdr[f'COMBIM{i+1}'] = os.path.basename(file)
        # 여기 다시 zp - zp_ref 해서 스케일링 하는걸로 변경
        hdr['NCOMBINE'] = int(len(filelist))
        if 'JD' in hdr.keys():
            hdr['JD'] = Time(np.mean([inim.header['JD'] for inim in ccdlist]), format='jd').value
        if 'DATE-OBS' in hdr.keys():
            hdr['DATE-OBS'] = Time(np.mean([Time(inim.header['DATE-OBS']).jd for inim in ccdlist]), format='jd').isot
        hdr['EXPTIME'] = float(np.sum([inim.header['EXPTIME'] for inim in ccdlist]))

        combiner = Combiner(ccdlist, dtype=np.float32)
        if scale == 'multiply':
            zp_median = np.median([inim.header[zp_key] for inim in ccdlist])
            for i, inim in enumerate(ccdlist):
                zp = inim.header[zp_key]
                zp_diff = zp - zp_median
                ccdlist[i].data = ccdlist[i].data * 100 ** (-zp_diff/5)
            
        elif (scale == 'zero'):
            averages = [np.mean(ccddata) for ccddata in ccdlist]
            delvalues = averages - averages[0]
            for i, delvalue in enumerate(delvalues):
                ccdlist[i].data = ccdlist[i].data-delvalue
                combiner = Combiner(ccdlist, dtype=np.float32)

        # Clipping
        if clip == 'minmax':
            combiner.minmax_clipping(min_clip=clip_minmax_min, max_clip=clip_minmax_max)
        if clip == 'sigma':
            combiner.sigma_clipping(low_thresh=clip_sigma_low, high_thresh=clip_sigma_high, func=np.ma.median)
        if clip == 'extrema':
            combiner.clip_extrema(nlow=clip_extrema_nlow, nhigh=clip_extrema_nhigh)

        # Combining
        if combine == 'median':
            combined = combiner.median_combine(median_func=self.bn_median)
        if combine == 'mean':
            combined = combiner.average_combine()
        if combine == 'sum':
            combined = combiner.sum_combine()

        combined.header = hdr

        outputname = f'{os.path.dirname(filelist[0])}/{prefix}{os.path.basename(filelist[0])}'
        if (len(filelist) == 1):
            ccd.header = hdr
            ccd.write(outputname, overwrite=True, format='fits')
        else:
            combined.write(outputname, overwrite=True, format='fits')
        init_mean, init_std = np.mean(ccdlist[0].data), np.std(ccdlist[0].data)
        fin_mean, fin_std = np.mean(combined.data), np.std(combined.data)

        print('Combine complete \n')
        print('Combine information')
        print(60*'=')
        print(f'Ncombine = {len(filelist)}')
        print(f'method   = {clip}(clipping), {combine}(combining)')
        print(f'mean     = {round(init_mean,3)} >>> {round(fin_mean,3)}')
        print(f'std      = {round(init_std,3)} >>> {round(fin_std,3)}')
        print(f'image path = {outputname}')
        return outputname

    def subtract_img(self,
                     target_img,
                     reference_img,
                     prefix='sub_',
                     method='hotpants',
                     # hotpants config
                     iu=60000,
                     il=-100000,
                     tu=600000000,
                     tl=-100000,
                     v=0,
                     ng='3 3 1.0 2 0.7 1 0.4'
                     ):
        '''
        parameters
        ----------
        1. target_img : str
                {abspath} of the target image to be subtracted
        2. reference_img : str
                {abspath} of the reference image 
        3. prefix : str
                prefix of the output image (sub_)
        4. method : str
                method for subtraction (hotpants)
        5. iu : int
                upper valid data count, image (60000)
        6. tu : int
                upper valid data count, template (600000000)
        7. tl : int
                lower valid data count, template (-100000)
        8. v : int
                level of verbosity, 0-2 (0)
        9. ng : str
                'ngauss degree0 sigma0 .. degreeN sigmaN'
                : ngauss = number of gaussians which compose kernel (3)
                : degree = degree of polynomial associated with gaussian # (3 2 1)
                : sigma  = width of gaussian # (1.0 0.7 0.4)
        returns 
        -------

        notes 
        -----
        For more information : https://github.com/acbecker/hotpants
        -----
        '''
        outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
        if method == 'hotpants':
            os.system(
                f'hotpants -c t -n i -inim {target_img} -tmplim {reference_img} -outim {outputname} -iu {iu} -il {il} -tu {tu} -tl {tl} -v {v} -ng {ng} > .out && rm -rf .out')
        return outputname

    # Program running
    @timeout(seconds = 15)
    def run_astrometry(self,
                       image, 
                       sex_configfile : str,
                       ra : float = None,
                       dec : float = None,
                       radius : float = None,
                       scalelow : float = 0.6, 
                       scalehigh : float = 0.8, 
                       overwrite : bool = False,
                       remove : bool = True
                       ):
        """
        1. Description
        : Solving WCS coordinates using Astrometry.net software. For better performance in especially B band images, --use-sextractor mode is added. This mode needs SExtractor configuration files. So please posit configuration files for your working directory. cpulimit 300 is also added to prevent too long processing time for bad images.
        : scalelow and scalehigh for the range of pixscale estimation

        2. History
        2018.03    Created by G.Lim.
        2018.12.18 Edited by G.Lim. SExtractor mode is added.
        2018.12.21 Edited by G.Lim. Define SAO_astrometry function.
        2020.03.01 --backend-config is added to have the system find INDEX files.
        2021.12.29 Edited by HH.Choi.  
        """
        import os,sys
        import glob
        import subprocess
        import numpy as np
        
        """
        Running the Astrometry process with options to pass RA/Dec and a timeout.
        """
        try:
            print('Start Astrometry process...=====================')
            # Set up directories and copy configuration files
            current_dir = os.getcwd()
            sex_dir = self.sexpath
            image_dir = os.path.dirname(image)
            os.chdir(sex_dir)
            os.system(f'cp {sex_configfile} {sex_dir}/*.param {sex_dir}/*.conv {sex_dir}/*.nnw {image_dir}')
            
            os.chdir(image_dir)
            print(f'Solving WCS using Astrometry with RA/Dec of {ra}/{dec} and radius of {radius} arcmin')

            # Building the command string
            if overwrite:
                new_filename = os.path.join(image_dir,os.path.basename(image))
                com = f'solve-field {image} --cpulimit 60 --overwrite --use-source-extractor --source-extractor-config {sex_configfile} --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low {str(scalelow)} --scale-high {str(scalehigh)} --no-remove-lines --uniformize 0 --no-plots --new-fits {new_filename} --temp-dir .'
            else:
                new_filename = os.path.join(image_dir,"a"+os.path.basename(image))
                com = f'solve-field {image} --cpulimit 60 --use-source-extractor --source-extractor-config {sex_configfile} --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low {str(scalelow)} --scale-high {str(scalehigh)} --no-remove-lines --uniformize 0 --no-plots --new-fits {new_filename} --temp-dir .'
            
            if ra is not None and dec is not None:
                com += f' --ra {ra} --dec {dec}'
            if radius is not None:
                com += f' --radius {radius}'
            
            print(f'Command to run: {com}')
            
            # Use subprocess.run with timeout
            result = subprocess.run(com, shell=True, timeout=900, check=True, text=True, capture_output=True)
            print(result.stdout)
            print(result.stderr)

            orinum = subprocess.check_output(f'ls C*.fits | wc -l', shell=True)
            resnum = subprocess.check_output(f'ls a*.fits | wc -l', shell=True)
            print(f"From {str(orinum[:-1])} files, {str(resnum[:-1])} files are solved.")

            # Clean up
            if remove:
                os.system(f'rm tmp* astrometry* *.conv default.nnw *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match check.fits *.param {os.path.basename(sex_configfile)}')
            print('Astrometry process finished=====================')
            return new_filename

        except subprocess.TimeoutExpired:
            print(f"The astrometry process exceeded the timeout limit.")
            return None
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the astrometry process: {e}")
            return None

    def run_sextractor(self, image, sex_configfile, sex_params: dict = None, return_result: bool = True, print_progress : bool = True):
        """
        Parameters
        ----------
        1. image : str
                Absolute path of the target image.
        2. sex_params : dict
                Configuration parameters in dict format. Can be loaded by load_sexconfig().
        3. sex_configfile : str
                Path to the SExtractor configuration file.
        4. return_result : bool
                If True, returns the result as an astropy table.

        Returns
        -------
        1. result : astropy.table.Table or str
                    Source extractor result as a table or the catalog file path.

        Notes
        -------
        This method runs SExtractor on the specified image using the provided configuration and parameters.
        """
        if print_progress:
            print('Start SExtractor process...=====================')

        # Switch to the SExtractor directory
        current_path = os.getcwd()
        os.chdir(self.sexpath)
        
        # Load and apply SExtractor parameters
        all_params = self.load_sexconfig(sex_configfile)
        sexparams_str = ''

        if sex_params:
            for key, value in sex_params.items():
                sexparams_str += f'-{key} {value} '
                all_params[key] = value

        # Command to run SExtractor
        command = f"source-extractor {image} -c {sex_configfile} {sexparams_str}"

        try:
            # Run the SExtractor command using subprocess.run
            subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if print_progress:
                print("SExtractor process finished=====================")

            if return_result:
                # Read the catalog produced by SExtractor
                sexresult = ascii.read(all_params['CATALOG_NAME'])
                os.chdir(current_path)
                return sexresult
            else:
                return all_params['CATALOG_NAME']
        except:
            if print_progress:
                print(f"Error during SExtractor execution")
            os.chdir(current_path)
            return None

        
    def run_scamp(self, filelist : str or list, sex_configfile : str, scamp_configfile : str = None, update_files : bool = True, print_progress : bool = True):
        
        if isinstance(filelist, str):
            filelist = [filelist]
        
        if print_progress:
            print(f'Start SCAMP process on {len(filelist)} images...=====================')
        sex_output_images = dict()
        for image in tqdm(filelist, desc='Running Source extractor...'):
            sex_params = dict()
            sex_params['CATALOG_NAME'] = f"{self.scamppath}/result/{os.path.basename(image).split('.')[0]}.sexcat"
            sex_params['PARAMETERS_NAME'] = f'{self.sexpath}/scamp.param'
            output_file = self.run_sextractor(image = image, sex_configfile = sex_configfile, sex_params = sex_params, return_result = False, print_progress = False)
            sex_output_images[image] = output_file
        
        if scamp_configfile is None:
            scamp_configfile = os.path.join(self.scamppath, 'default.scamp')
        
        # Sort out if the sextractor failed
        sex_output_images = {key: value for key, value in sex_output_images.items() if value is not None}
        scamp_output_images = {key: value.replace('.sexcat', '.head') for key, value in sex_output_images.items()}

        all_images_str = ' '.join(sex_output_images.values())
        command = f'scamp {all_images_str} -c {scamp_configfile}'
        
        try:
            current_path = os.getcwd()
            os.chdir(os.path.join(self.scamppath,'result'))
            # Run the SExtractor command using subprocess.run
            subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if print_progress:
                print("SCAMP process finished=====================")

            if update_files:
                def sanitize_header(header: fits.Header) -> fits.Header:
                    """
                    Sanitize a FITS header by removing or cleaning non-ASCII and non-printable characters.
                    
                    Parameters:
                    header (fits.Header): The FITS header to be sanitized.
                    
                    Returns:
                    fits.Header: The sanitized header.
                    """
                    sanitized_header = fits.Header()
                    
                    # Loop through each header card and sanitize it
                    for card in header.cards:
                        key, value, comment = card
                        if isinstance(value, str):
                            # Remove non-ASCII characters from the value
                            value = re.sub(r'[^\x20-\x7E]+', '', value)
                        
                        # Add sanitized card to the new header
                        sanitized_header[key] = (value, comment)
                    
                    return sanitized_header

                def update_fits_with_head(image_file: str, head_file: str):
                    """
                    Update the WCS and other relevant header information in a FITS file using a SCAMP-generated .head file.
                    
                    Parameters:
                    image_file (str): Path to the FITS image file to be updated.
                    head_file (str): Path to the SCAMP-generated .head file with updated WCS and other parameters.
                    """
                    # Read the header from the .head file
                    with open(head_file, 'r') as head:
                        head_content = head.read()

                    # Convert the head file content to an astropy header object
                    head_header = fits.Header.fromstring(head_content, sep='\n')
                    
                    # Sanitize the header to remove non-ASCII characters
                    head_header = sanitize_header(head_header)

                    # Open the FITS image and update its header with WCS information from the .head file
                    with fits.open(image_file, mode='update') as hdu:
                        # Update the header of the primary HDU (index 0) with the sanitized WCS-related keys
                        hdu[0].header.update(head_header)
                        hdu.flush()  # Write changes to disk

                    print(f"Updated WCS and relevant header information for {image_file} using {head_file}")

                
                for image, header in scamp_output_images.items():
                    update_fits_with_head(image, header)
            else:
                return scamp_output_images.values()
        except:
            print(f"Error during SCAMP execution")
            return
        finally:
            os.chdir(current_path)

    def run_ds9(self, filelist, shell: str = '/bin/bash'):
        import subprocess
        import numpy as np
        '''
		parameters
		----------
		1. filelist : str or list or np.array
				(str) abspath of the fitsfile for visualization
				(list or np.array) the list of abspath of the fitsfiles for visualization
		
		returns 
		-------
		
		notes 
		-----
		-----
		'''
        ds9_options = "-scalemode zscale -scale lock yes -frame lock image "
        names = ""
        if (type(filelist) == str) | (type(filelist) == np.str_):
            names = filelist
        else:
            for file in filelist:
                names += file+" "
        ds9_command = "ds9 "+ds9_options+names+" &"
        print('Running "'+ds9_command+'" in the terminal...')
        sp = subprocess.Popen([shell, "-i", "-c", ds9_command])
        sp.communicate()
        # os.system(ds9_command)

    def to_regions(self, reg_ra, reg_dec, reg_size: float = 5.0, output_file_name: str = 'regions.reg'):
        from regions import CircleSkyRegion, write_ds9
        import astropy.units as u
        # Check if ra and dec are single float values or lists
        if isinstance(reg_ra, float) and isinstance(reg_dec, float):
            ra_list = [reg_ra]
            dec_list = [reg_dec]
        elif isinstance(reg_ra, list) and isinstance(reg_dec, list):
            ra_list = reg_ra
            dec_list = reg_dec
        else:
            ra_list = reg_ra
            dec_list = reg_dec

        regions = []
        for ra, dec in zip(ra_list, dec_list):
            center = SkyCoord(ra, dec, unit='deg', frame='icrs')
            radius = reg_size * u.arcsec  # Example radius
            region = CircleSkyRegion(center, radius)
            regions.append(region)
        # Write the regions to a DS9 region file
        output_file_path = os.path.join(self.photpath, output_file_name)
        write_ds9(regions, output_file_path)
        return output_file_path
    
    def visualize_image(self, filename : str):
        from astropy.visualization import ImageNormalize, ZScaleInterval
        import matplotlib.pyplot as plt
        data = fits.getdata(filename)
        zscale_interval = ZScaleInterval()
        norm = ImageNormalize(data, interval=zscale_interval)

        # Plot the normalized image
        plt.imshow(data, cmap='gray', norm=norm, origin='lower')
        plt.colorbar()
        plt.title(f'{os.path.basename(filename)}')
        plt.show()
        

# %%
if __name__ == '__main__':
    import glob
    A = PhotometryHelper()
    file_ = '/data1/supernova_rawdata/SN2023rve/analysis/RASA36/reference_image/com_align_Calib-RASA36-NGC1097-20210719-091118-r-60.fits'
    #file1 = '/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/r/align_com_align_cutoutmost_Calib-KCT_STX16803-NGC1097-20230801-075323-r-120.fits'
    #file2 = '/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/reference_image/cut_Ref-KCT_STX16803-NGC1097-r-5400.com.fits'
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_HIGH.scampconfig'
    filelist = glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/RASA36/reference_image_tmp/cutout*.fits')
    #sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/KCT.config'
    #A.run_astrometry(image = file1, sex_configfile = sex_configfile)
    #A.run_scamp(filelist = file_, sex_configfile = sex_configfile)
    
    #sex_params = dict()
    #sex_params['CATALOG_NAME'] = f"{A.scamppath}/catalog/{os.path.basename(file_).split('.')[0]}.cat"
    #sex_params['PARAMETERS_NAME'] = f'{A.sexpath}/scamp.param'
    target_img = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/g/Calib-KCT_STX16803-NGC1097-20230927-062948-g-120.fits'
    reference_img = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/g/Calib-KCT_STX16803-NGC1097-20230927-063834-g-120.fits'
    A.visualize_image(target_img)
    A.visualize_image(reference_img)
    #A.align_img(target_img, reference_img)
# %%
