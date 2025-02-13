
#%%
from Research.helper import Helper
from Research.photometry import ScienceImage
from Research.photometry import CalibrationImage
from astropy.time import Time
from typing import List, Dict, Union
import ccdproc
from ccdproc import CCDData
import astropy.units as u
import numpy as np
from tqdm import tqdm
#%%
class Preprocess(Helper):
    
    def __init__(self):
        super().__init__()
    
    def correct_bdf(self, tgt_image: ScienceImage, bias_image: CalibrationImage, dark_image: CalibrationImage, flat_image: CalibrationImage,
                    remove : bool = False
                    ):
        if tgt_image.status.biascor['status']:
            tgt_image.logger.warning(f"BIAS correction already applied to {tgt_image.path}. BIAS correction is not applied.")
            raise RuntimeError(f"BIAS correction already applied to {tgt_image.path}")
        if tgt_image.status.darkcor['status']:
            tgt_image.logger.warning(f"DARK correction already applied to {tgt_image.path}. DARK correction is not applied.")
            raise RuntimeError(f"DARK correction already applied to {tgt_image.path}")
        if tgt_image.status.flatcor['status']:
            tgt_image.logger.warning(f"FLAT correction already applied to {tgt_image.path}. FLAT correction is not applied.")
            raise RuntimeError(f"FLAT correction already applied to {tgt_image.path}")
        
        # Convert input images to CCDData
        sci_ccddata = ccdproc.CCDData(data = tgt_image.data, meta = tgt_image.header, unit = 'adu')
        bias_ccddata = ccdproc.CCDData(data = bias_image.data, meta = bias_image.header, unit = 'adu')
        dark_ccddata = ccdproc.CCDData(data = dark_image.data, meta = dark_image.header, unit = 'adu')
        flat_ccddata = ccdproc.CCDData(data = flat_image.data, meta = flat_image.header, unit = 'adu')

        # Perform bias, dark, flat correction
        calib_data = self._correct_bdf(tgt_data = sci_ccddata, bias_data = bias_ccddata, dark_data = dark_ccddata, flat_data = flat_ccddata)        

        # Determine data types and convert to selected data type
        tgt_dtype = tgt_image.data.dtype
        bias_dtype = bias_image.data.dtype
        dark_dtype = dark_image.data.dtype
        flat_dtype = flat_image.data.dtype
        selected_dtype = np.promote_types(tgt_dtype, bias_dtype, dark_dtype, flat_dtype)
        calib_data.data = calib_data.data.astype(selected_dtype)
        
        # Add metadata
        calib_data.meta['BIASCOR'] = True
        calib_data.meta['BCORTIME'] = Time.now().isot
        calib_data.meta['BIASPATH'] = bias_image.path
        calib_data.meta['DARKCOR'] = True
        calib_data.meta['DCORTIME'] = Time.now().isot
        calib_data.meta['DARKPATH'] = dark_image.path
        calib_data.meta['FLATCOR'] = True
        calib_data.meta['FCORTIME'] = Time.now().isot
        calib_data.meta['FLATPATH'] = flat_image.path
        
        # Determine output filename
        dirpath = os.path.dirname(tgt_image.path)
        filename = os.path.basename(tgt_image.path)
        if 'cor_' in filename:
            filepath = os.path.join(dirpath, 'fdb' + filename)
        else:
            filepath = os.path.join(dirpath, 'fdbcor_' + filename)
        calib_data.write(filepath, overwrite = True)
        
        # Create new image object
        calib_image = type(tgt_image)(path  = filepath, telinfo = tgt_image.telinfo)
        calib_image.logger.info(f"BIAS, DARK, FLAT correction applied with {bias_image.path}, {dark_image.path}, {flat_image.path}")
        calib_image.status = tgt_image.status
        calib_image.update_status(process_name = 'biascor')
        calib_image.update_status(process_name = 'darkcor')
        calib_image.update_status(process_name = 'flatcor')
        
        # Log information
        tgt_image.logger.info(f"BIAS, DARK, FLAT correction applied: FILEPATH = {calib_image.path}")
        bias_image.logger.info(f"Used for BIAS correction: FILEPATH = {calib_image.path}")
        dark_image.logger.info(f"Used for DARK correction: FILEPATH = {calib_image.path}")
        flat_image.logger.info(f"Used for FLAT correction: FILEPATH = {calib_image.path}")
        if remove:
            os.remove(tgt_image.path)
            os.remove(tgt_image.loggerpath)
            os.remove(tgt_image.statuspath)
        return calib_image
    
    def _correct_bdf(self, tgt_data : CCDData, bias_data : CCDData, dark_data : CCDData, flat_data : CCDData):
        bcalib_data = self._correct_bias(tgt_data = tgt_data, bias_data = bias_data)
        dbcalib_data = self._correct_dark(tgt_data = bcalib_data, dark_data = dark_data)
        fdbcalib_data = self._correct_flat(tgt_data = dbcalib_data, flat_data = flat_data)
        return fdbcalib_data
    
    def correct_db(self, tgt_image: ScienceImage, bias_image: CalibrationImage, dark_image: CalibrationImage, 
                    remove : bool = False
                    ):
        if tgt_image.status.biascor['status']:
            tgt_image.logger.warning(f"BIAS correction already applied to {tgt_image.path}. BIAS correction is not applied.")
            raise RuntimeError(f"BIAS correction already applied to {tgt_image.path}")
        if tgt_image.status.darkcor['status']:
            tgt_image.logger.warning(f"DARK correction already applied to {tgt_image.path}. DARK correction is not applied.")
            raise RuntimeError(f"DARK correction already applied to {tgt_image.path}")
        
        # Convert input images to CCDData
        sci_ccddata = ccdproc.CCDData(data = tgt_image.data, meta = tgt_image.header, unit = 'adu')
        bias_ccddata = ccdproc.CCDData(data = bias_image.data, meta = bias_image.header, unit = 'adu')
        dark_ccddata = ccdproc.CCDData(data = dark_image.data, meta = dark_image.header, unit = 'adu')

        # Perform bias, dark correction
        calib_data = self._correct_bd(tgt_data = sci_ccddata, bias_data = bias_ccddata, dark_data = dark_ccddata)

        # Determine data types and convert to selected data type
        tgt_dtype = tgt_image.data.dtype
        bias_dtype = bias_image.data.dtype
        dark_dtype = dark_image.data.dtype
        selected_dtype = np.promote_types(tgt_dtype, bias_dtype, dark_dtype)
        calib_data.data = calib_data.data.astype(selected_dtype)
        
        # Add metadata
        calib_data.meta['BIASCOR'] = True
        calib_data.meta['BCORTIME'] = Time.now().isot
        calib_data.meta['BIASPATH'] = bias_image.path
        calib_data.meta['DARKCOR'] = True
        calib_data.meta['DCORTIME'] = Time.now().isot
        calib_data.meta['DARKPATH'] = dark_image.path
        
        # Determine output filename
        dirpath = os.path.dirname(tgt_image.path)
        filename = os.path.basename(tgt_image.path)
        if 'cor_' in filename:
            filepath = os.path.join(dirpath, 'db' + filename)
        else:
            filepath = os.path.join(dirpath, 'dbcor_' + filename)
        calib_data.write(filepath, overwrite = True)
        
        # Create new image object
        calib_image = type(tgt_image)(path  = filepath, telinfo = tgt_image.telinfo)
        calib_image.logger.info(f"BIAS, DARK correction applied with {bias_image.path}, {dark_image.path}")
        calib_image.status = tgt_image.status
        calib_image.update_status(process_name = 'biascor')
        calib_image.update_status(process_name = 'darkcor')
        
        # Log information
        tgt_image.logger.info(f"BIAS, DARK correction applied: FILEPATH = {calib_image.path}")
        bias_image.logger.info(f"Used for BIAS correction: FILEPATH = {calib_image.path}")
        dark_image.logger.info(f"Used for DARK correction: FILEPATH = {calib_image.path}")
        if remove:
            os.remove(tgt_image.path)
            os.remove(tgt_image.loggerpath)
            os.remove(tgt_image.statuspath)
        return calib_image
    
    def _correct_bd(self, tgt_data : CCDData, bias_data : CCDData, dark_data : CCDData):
        bcalib_data = self._correct_bias(tgt_data = tgt_data, bias_data = bias_data)
        dbcalib_data = self._correct_dark(tgt_data = bcalib_data, dark_data = dark_data)
        return dbcalib_data
        
    def correct_bias(self, tgt_image: ScienceImage or CalibrationImage, bias_image: CalibrationImage,
                     remove : bool = False
                     ):
        """ Corrects bias in the image """
        if tgt_image.status.biascor['status']:
            tgt_image.logger.warning(f"BIAS correction already applied to {tgt_image.path}. BIAS correction is not applied.")
            raise RuntimeError(f"BIAS correction already applied to {tgt_image.path}")
        
        # Convert input images to CCDData
        sci_ccddata = ccdproc.CCDData(data = tgt_image.data, meta = tgt_image.header, unit = 'adu')
        bias_ccddata = ccdproc.CCDData(data = bias_image.data, meta = bias_image.header, unit = 'adu')
        
        # Perform bias correction
        calib_data = self._correct_bias(tgt_data = sci_ccddata, bias_data = bias_ccddata)
        
        # Determine data types and convert to selected data type
        tgt_dtype = tgt_image.data.dtype
        bias_dtype = bias_image.data.dtype
        selected_dtype = np.promote_types(tgt_dtype, bias_dtype)
        calib_data.data = calib_data.data.astype(selected_dtype)

        # Add metadata
        calib_data.meta['BIASCOR'] = True
        calib_data.meta['BCORTIME'] = Time.now().isot
        calib_data.meta['BIASPATH'] = bias_image.path
        
        # Determine output filename
        dirpath = os.path.dirname(tgt_image.path)
        filename = os.path.basename(tgt_image.path)
        if 'cor_' in filename:
            filepath = os.path.join(dirpath, 'b' + filename)
        else:
            filepath = os.path.join(dirpath, 'bcor_' + filename)          
        calib_data.write(filepath, overwrite = True)
        
        # Create new image object
        calib_image = type(tgt_image)(path  = filepath, telinfo = tgt_image.telinfo)
        calib_image.logger.info(f"BIAS correction applied with {bias_image.path}")
        calib_image.status = tgt_image.status
        calib_image.update_status(process_name = 'biascor')
        
        # Log information
        tgt_image.logger.info(f"BIAS correction applied: FILEPATH = {calib_image.path}")
        bias_image.logger.info(f"Used for BIAS correction: FILEPATH = {calib_image.path}")
        if remove:
            os.remove(tgt_image.path)
            os.remove(tgt_image.loggerpath)
            os.remove(tgt_image.statuspath)
        return calib_image

    
    def _correct_bias(self, tgt_data : CCDData, bias_data : CCDData):
        calib_data = ccdproc.subtract_bias(tgt_data, bias_data)
        return calib_data
    
    def correct_dark(self, tgt_image: ScienceImage or CalibrationImage, dark_image: CalibrationImage, 
                     remove : bool = False
                     ):
        """ Corrects dark in the image """
        if tgt_image.status.darkcor['status']:
            tgt_image.logger.warning(f"DARK correction already applied to {tgt_image.path}. DARK correction is not applied.")
            raise RuntimeError(f"DARK correction already applied to {tgt_image.path}")
        
        # Convert input images to CCDData
        sci_ccddata = ccdproc.CCDData(data = tgt_image.data, meta = tgt_image.header, unit = 'adu')
        dark_ccddata = ccdproc.CCDData(data = dark_image.data, meta = dark_image.header, unit = 'adu')
        
        # Perform dark correction
        calib_data = self._correct_dark(tgt_data = sci_ccddata, dark_data = dark_ccddata)
        
        # Determine data types and convert to selected data type
        tgt_dtype = tgt_image.data.dtype
        dark_dtype = dark_image.data.dtype
        selected_dtype = np.promote_types(tgt_dtype, dark_dtype)
        calib_data.data = calib_data.data.astype(selected_dtype)
        
        # Add metadata
        calib_data.meta['DARKCOR'] = True
        calib_data.meta['DCORTIME'] = Time.now().isot
        calib_data.meta['DARKPATH'] = dark_image.path
        
        # Determine output filename
        dirpath = os.path.dirname(tgt_image.path)
        filename = os.path.basename(tgt_image.path)
        if 'cor_' in filename:
            filepath = os.path.join(dirpath, 'd' + filename)
        else:
            filepath = os.path.join(dirpath, 'dcor_' + filename)           
        calib_data.write(filepath, overwrite = True)
        
        # Create new image object
        calib_image = type(tgt_image)(path  = filepath, telinfo = tgt_image.telinfo)
        calib_image.logger.info(f"DARK correction applied with {dark_image.path}")
        calib_image.status = tgt_image.status
        calib_image.update_status(process_name = 'darkcor')
        
        # Log information
        tgt_image.logger.info(f"DARK correction applied: FILEPATH = {calib_image.path}")
        dark_image.logger.info(f"Used for DARK correction: FILEPATH = {calib_image.path}")
        if remove:
            os.remove(tgt_image.path)
            os.remove(tgt_image.loggerpath)
            os.remove(tgt_image.statuspath)
        return calib_image

    def _correct_dark(self, tgt_data : CCDData, dark_data : CCDData):
        calib_data = ccdproc.subtract_dark(tgt_data, dark_data, scale = True, exposure_time = 'EXPTIME', exposure_unit = u.second)
        return calib_data
    
    def correct_flat(self, tgt_image: ScienceImage, flat_image: CalibrationImage,
                     remove : bool = False
                     ):
        if tgt_image.status.flatcor['status']:
            tgt_image.logger.warning(f"FLAT correction already applied to {tgt_image.path}. FLAT correction is not applied.")
            raise RuntimeError(f"FLAT correction already applied to {tgt_image.path}")
        
        # Convert input images to CCDData
        sci_ccddata = ccdproc.CCDData(data = tgt_image.data, meta = tgt_image.header, unit = 'adu')
        flat_ccddata = ccdproc.CCDData(data = flat_image.data, meta = flat_image.header, unit = 'adu')
        
        # Perform flat correction
        calib_data = self._correct_flat(tgt_data = sci_ccddata, flat_data = flat_ccddata)
        
        # Determine data types and convert to selected data type
        tgt_dtype = tgt_image.data.dtype
        flat_dtype = flat_image.data.dtype
        selected_dtype = np.promote_types(tgt_dtype, flat_dtype)
        calib_data.data = calib_data.data.astype(selected_dtype)
        
        # Add metadata
        calib_data.meta['FLATCOR'] = True
        calib_data.meta['FCORTIME'] = Time.now().isot
        calib_data.meta['FLATPATH'] = flat_image.path
        
        # Determine output filename
        dirpath = os.path.dirname(tgt_image.path)
        filename = os.path.basename(tgt_image.path)
        if 'cor_' in filename:
            filepath = os.path.join(dirpath, 'f' + filename)
        else:
            filepath = os.path.join(dirpath, 'fcor_' + filename)
        calib_data.write(filepath, overwrite = True)
        
        # Create new image object
        calib_image = type(tgt_image)(path  = filepath, telinfo = tgt_image.telinfo)
        calib_image.logger.info(f"FLAT correction applied with {flat_image.path}")
        calib_image.status = tgt_image.status
        calib_image.update_status(process_name = 'flatcor')
        
        # Log information
        tgt_image.logger.info(f"FLAT correction applied: FILEPATH = {calib_image.path}")
        flat_image.logger.info(f"Used for FLAT correction: FILEPATH = {calib_image.path}")
        if remove:
            os.remove(tgt_image.path)
            os.remove(tgt_image.loggerpath)
            os.remove(tgt_image.statuspath)
        return calib_image
        
    def _correct_flat(self, tgt_data : CCDData, flat_data : CCDData):
        calib_data = ccdproc.flat_correct(tgt_data, flat_data)
        return calib_data
    
    def generate_master_frame(self, calib_imagelist : List[CalibrationImage], 
                              mbias : Union[CalibrationImage or List[CalibrationImage],None],
                              mdark : Union[CalibrationImage or List[CalibrationImage],None]):
        """ Generate master bias, dark, flat frames """
        all_filelist = [image.path for image in calib_imagelist]
        all_fileinfo = self.get_imginfo(all_filelist, keywords = ['imagetyp', 'exptime', 'xbinning', 'filter', 'gain', 'date-obs'])
        all_fileinfo['image'] = calib_imagelist
        all_fileinfo_by_group = all_fileinfo.group_by(['xbinning', 'gain']).groups
        master_files = dict()
        for group in all_fileinfo_by_group:
            key = (group['xbinning'][0], group['gain'][0])
            master_files[key] = dict(BIAS = None, DARK = None, FLAT = dict())
            
        if mbias:
            if isinstance(mbias, CalibrationImage):
                mbias = [mbias]
            for bias in mbias:
                header = bias.header
                bias_key = (header['xbinning'], header['gain'])
                master_files[bias_key]['BIAS'] = bias
        if mdark:
            if isinstance(mdark, CalibrationImage):
                mdark = [mdark]
            for dark in mdark:
                header = dark.header
                dark_key = (header['xbinning'], header['gain'])
                master_files[dark_key]['DARK'] = dark
        
        # Run the calibration
        for group in all_fileinfo_by_group:
            # Separate the images by type
            key = (group['xbinning'][0], group['gain'][0])
            if not master_files[key]['BIAS']:
                bias_key = ['BIAS', 'ZERO']
                bias_mask  = np.isin(group['imagetyp'], bias_key)
                bias_fileinfo = group[bias_mask]
                if not bias_fileinfo:
                    raise ValueError("No BIAS or ZERO frames found.")
                date_str = np.mean(Time(bias_fileinfo['date-obs'])).datetime.strftime('%Y%m%d')
                output_name = f'{date_str}-zero.fits'
                output_name = self.combine_img(filelist = bias_fileinfo['file'], 
                                                combine_method = 'median', 
                                                scale = None, 
                                                output_name = output_name,
                                                print_output = True,
                                                clip = 'extrema',
                                                clip_extrema_nlow=1,
                                                clip_extrema_nhigh=1)
                mbias = CalibrationImage(path = output_name, telinfo = bias_fileinfo[0].telinfo)
                master_files[key]['BIAS'] = mbias
                    
            if not master_files[key]['DARK']:
                dark_key = ['DARK']
                dark_mask  = np.isin(group['imagetyp'], dark_key)
                dark_fileinfo = group[dark_mask]
                if not dark_fileinfo:
                    pass
                else:
                    date_str = np.mean(Time(dark_fileinfo['date-obs'])).datetime.strftime('%Y%m%d')
                    output_name = f'{date_str}-dark.fits'
                    b_darkimagelist = []
                    for dark in tqdm(dark_fileinfo['image'], desc = 'BIAS correction on DARK frames...'):
                        b_dark_image = self.correct_bias(tgt_image = dark, bias_image = master_files[key]['BIAS'], remove = False)
                        b_darkimagelist.append(b_dark_image)
                    b_darkfilelist = [image.path for image in b_darkimagelist]
                        
                    output_name = self.combine_img(filelist = b_darkfilelist, 
                                                combine_method = 'median', 
                                                scale = None, 
                                                output_name = output_name,
                                                print_output = True,
                                                clip = 'extrema',
                                                clip_extrema_nlow=1,
                                                clip_extrema_nhigh=1)
                    mdark = CalibrationImage(path = output_name, telinfo = dark_imagelist[0].telinfo)
                    master_files[key]['DARK'] = mdark
            
            if not master_files[key]['FLAT']:    
                if (not master_files[key]['DARK']) or (not master_files[key]['BIAS']):
                    raise ValueError("Master BIAS or DARK frame not found.")

                flat_key = ['FLAT']
                flat_mask = np.isin(group['imagetyp'], flat_key)
                flat_fileinfo = group[flat_mask]
                if not flat_fileinfo:
                    pass
                else:
                    flat_fileinfo_by_filter = flat_fileinfo.group_by('filter').groups
                    for filter_group in flat_fileinfo_by_filter:
                        filter_name = filter_group['filter'][0]
                        date_str = np.mean(Time(filter_group['date-obs'])).datetime.strftime('%Y%m%d')
                        output_name = f'{date_str}-flat-{filter_name}.fits'
                        db_flatimagelist = []
                        for flat in tqdm(filter_group['image'], desc = 'BIAS, DARK correction on FLAT frames...'):
                            db_flat_image = self.correct_db(tgt_image = flat, bias_image = master_files[key]['BIAS'], dark_image = master_files[key]['DARK'], remove = False)
                            db_flatimagelist.append(db_flat_image)
                        db_flatfilelist = [image.path for image in db_flatimagelist]
                    
                        output_name = self.combine_img(filelist = db_flatfilelist, 
                                                       combine_method = 'median', 
                                                       scale = None, 
                                                       output_name = output_name,
                                                       print_output = True,
                                                       clip = 'extrema',
                                                       clip_extrema_nlow=1,
                                                       clip_extrema_nhigh=1)
                        mflat = CalibrationImage(path = output_name, telinfo = flat_imagelist[0].telinfo)
                        master_files[key]['FLAT'][filter_name] = mflat
        return master_files
            
        
        
        
    
# %%
P = Preprocess()
# %%
from tqdm import tqdm
import glob
filelist = glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/test/rawdata/SNU*-g-*.fts')
bias_path = '/data1/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/test/zero/20210914-zero.fits'
dark_path = '/data1/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/test/dark/120-20210804-dark.fits'
flat_path = '/data1/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/test/flat/20210804-ng.fits'
bias = CalibrationImage(bias_path, telinfo = Helper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'))
dark = CalibrationImage(dark_path, telinfo = Helper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'))
flat = CalibrationImage(flat_path, telinfo = Helper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'))
for image in tqdm(filelist):
    sci = ScienceImage(image, telinfo = Helper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'))
    #bsci = P.correct_bias(tgt_image = sci, bias_image = bias, remove = False)
    #dbsci = P.correct_dark(tgt_image = bsci, dark_image = dark, remove = True)
    #fdbsci = P.correct_flat(tgt_image = dbsci, flat_image = flat, remove = True)
    fdbsci = P.correct_bdf(tgt_image = sci, bias_image = bias, dark_image = dark, flat_image = flat, remove = False)
# %%
import glob
tel_name = '7DT02'
biaslist = glob.glob(f'/data1/7DT/S240422ed/obsdata/{tel_name}/*BIAS*.fits')
darklist = glob.glob(f'/data1/7DT/S240422ed/obsdata/{tel_name}/*DARK*.fits')
flatlist = glob.glob(f'/data1/7DT/S240422ed/obsdata/{tel_name}/*FLAT*.fits')
caliblist = biaslist + darklist + flatlist
calib_imagelist = [CalibrationImage(path = file_, telinfo = Helper().get_telinfo(telescope = '7DT', ccd = 'C361K', readoutmode = 'high')) for file_ in caliblist]
# %%
self = Preprocess()
mbias = CalibrationImage(path = '/data1/7DT/S240422ed/obsdata/7DT02/20240423-zero.fits', telinfo = tgt_imagelist[0].telinfo)
# %%
