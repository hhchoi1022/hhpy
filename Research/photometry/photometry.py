#%%
import glob
import os
from datetime import datetime
import json
import multiprocessing
import gc
import numpy as np
from tqdm import tqdm

from Research.helper import PhotometryHelper
from Research.photometry import Image
#%%

class Photometry(PhotometryHelper):
    
    def __init__(self, 
                 imagelist : list, 
                 telescope_info: dict = None,
                 reference_image : str = None):
        super().__init__()
        self.original_imagelist = imagelist
        self.target_imagelist = imagelist
        self.reference_image = reference_image
        if telescope_info == None:
            telescope_info = self.get_telinfo()
            
        self.telinfo = telescope_info
        self.failed_imagelist = dict()
        self.failed_imagelist['Execution Time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.failed_imagelist['calculate'] = list()
        self.failed_imagelist['align'] = list()
        self.failed_imagelist['combine'] = list()
        self.failed_imagelist['photometry'] = list()

    def astrometry(self, 
                   sex_configfile : str,
                   filelist : list = None,
                   ra : float = None,
                   dec : float = None,
                   radius : float = None,
                   scalelow : float = 0.6, 
                   scalehigh : float = 0.8, 
                   overwrite : bool = False, 
                   remove : bool = True
                   ):
        
        if filelist is None:
            filelist = self.target_imagelist
        
        for file_ in filelist:
            im = Image(image=file_, telescope_info=self.telinfo, reference_image=self.reference_image)
            try:
                im.astrometry(sex_configfile = sex_configfile,
                              scalelow = scalelow,
                              scalehigh = scalehigh,
                              overwrite = overwrite,
                              ra = ra,
                              dec = dec,
                              radius = radius)
            except:
                self.failed_imagelist['astrometry'].append(file_)
        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_astrometry_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_astrometry_failed_imagelist.txt")

    def faster_astrometry(self,
                          sex_configfile : str,
                          num_processes : int = 3,
                          filelist : list = None,
                          ra : float = None,
                          dec : float = None,
                          radius : float = None,
                          scalelow : float = 0.6, 
                          scalehigh : float = 0.8, 
                          overwrite : bool = False, 
                          remove : bool = True
                          ):
        if filelist is None:
            filelist = self.target_imagelist
            
        # Split the target image list into groups
        print(f"Splitting image list into {num_processes} groups...")
        num_images = len(self.target_imagelist)
        num_images_per_group = num_images // num_processes
        image_groups = [self.target_imagelist[i * num_images_per_group:(i + 1) * num_images_per_group] for i in range(num_processes - 1)]
        image_groups.append(self.target_imagelist[(num_processes - 1) * num_images_per_group:])
        
        # Create a pool of processes and map the worker function to each group
        print('Run astrometry with {} processes'.format(num_processes))
        with multiprocessing.Pool(num_processes) as pool:
            results = []
            for group in image_groups:
                args = (sex_configfile, group, ra, dec, radius, scalelow, scalehigh, overwrite, False)
                result = pool.apply_async(self._worker_calculate, args=args)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]
        
        if remove:
            image_paths = set([os.path.dirname(file) for file in filelist])
            for path in image_paths:
                os.chdir(path)
                os.system(f'rm tmp* astrometry* *.conv default.nnw *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match check.fits *.param {os.path.basename(sex_configfile)}')

        return output
    
    def _worker_astrometry(self, filelist, sex_configfile, ra, dec, radius, scalelow, scalehigh, overwrite, remove=False):
        self.astrometry(filelist=filelist,
                        sex_configfile=sex_configfile,
                        ra=ra,
                        dec=dec,
                        radius=radius,
                        scalelow=scalelow,
                        scalehigh=scalehigh,
                        overwrite=overwrite,
                        remove=False)

    def calculate(self, 
                  filelist: list = None,
                  sex_configfile: str = None,  # Absolute Path
                  detect_threshold: float = 3.0,
                  aperture_type: str = 'relative',  # relative or absolute
                  aperture_sizes: list = [1.5, 2.5, 3.5],  # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                  ref_catalog_name: str = 'APASS',
                  ref_catalog_conversion : str = 'SDSS',
                  ref_maxmag=16,
                  ref_minmag=12,
                  visualize: bool = True,
                  update_header: bool = True):
        
        if filelist is None:
            filelist = self.target_imagelist
            
        for file_ in filelist:
            im = Image(image=file_, telescope_info=self.telinfo, reference_image=self.reference_image)
            try:
                im.calculate_zeropoint(sex_configfile=sex_configfile, 
                                       detect_threshold=detect_threshold, 
                                       aperture_type=aperture_type, 
                                       aperture_sizes=aperture_sizes, 
                                       ref_catalog_name=ref_catalog_name, 
                                       ref_catalog_conversion=ref_catalog_conversion, 
                                       ref_maxmag=ref_maxmag, 
                                       ref_minmag=ref_minmag, 
                                       visualize=visualize, 
                                       update_header=update_header)
            except:
                self.failed_imagelist['calculate'].append(file_)
        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_calculate_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_calculate_failed_imagelist.txt")

    def faster_calculate(self, 
                         filelist : list = None,
                         sex_configfile: str = None,  # Absolute Path
                         detect_threshold: float = 3.0,
                         aperture_type: str = 'relative',  # relative or absolute
                         aperture_sizes: list = [1.5, 2.5, 3.5],  # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                         ref_catalog_name: str = 'SMSS',
                         ref_catalog_conversion: str = 'SDSS',
                         ref_maxmag=16,
                         ref_minmag=12,
                         visualize: bool = True,
                         update_header: bool = True,
                         num_processes: int = 4
                         ):

        if filelist is None:
            filelist = self.target_imagelist

        # Split the target image list into groups
        print(f"Splitting image list into {num_processes} groups...")
        num_images = len(self.target_imagelist)
        num_images_per_group = num_images // num_processes
        image_groups = [self.target_imagelist[i * num_images_per_group:(i + 1) * num_images_per_group] for i in range(num_processes - 1)]
        image_groups.append(self.target_imagelist[(num_processes - 1) * num_images_per_group:])
        
        # Create a pool of processes and map the worker function to each group
        print('Calculating zeropoints with {} processes'.format(num_processes))
        with multiprocessing.Pool(num_processes) as pool:
            results = []
            for group in image_groups:
                args = (group, sex_configfile, detect_threshold, aperture_type, aperture_sizes, ref_catalog_name, ref_catalog_conversion, ref_maxmag, ref_minmag, visualize, update_header)
                result = pool.apply_async(self._worker_calculate, args=args)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]
            
        return output

    def _worker_calculate(self, filelist, sex_configfile, detect_threshold, aperture_type, aperture_sizes, ref_catalog_name, ref_catalog_conversion, ref_maxmag, ref_minmag, visualize, update_header):
        self.calculate(filelist=filelist, 
                       sex_configfile=sex_configfile, 
                       detect_threshold=detect_threshold, 
                       aperture_type=aperture_type, 
                       aperture_sizes=aperture_sizes, 
                       ref_catalog_name=ref_catalog_name, 
                       ref_catalog_conversion=ref_catalog_conversion, 
                       ref_maxmag=ref_maxmag, 
                       ref_minmag=ref_minmag, 
                       visualize=visualize, 
                       update_header=update_header)
        
    def exclude_outlier(self, 
                        sigma : float = 5,
                        depth_cut : float = None,
                        seeing_cut : float = None,
                        depth_key : str = 'DEPTH5_APER_1_HH',
                        seeing_key : str = 'SEEING_HH',
                        move : bool = False, 
                        visualize : bool = True):
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt
        all_tbl = self.get_imginfo(filelist = self.target_imagelist, keywords = [depth_key, seeing_key])
        masked_index = (all_tbl[depth_key].mask) | (all_tbl[seeing_key].mask)
        masked_tbl = all_tbl[masked_index]
        all_tbl = all_tbl[~masked_index]
        if depth_cut is not None:
            select_index = (np.array(all_tbl[depth_key]).astype(float) ) > depth_cut
        if seeing_cut is not None:
            select_index = (np.array(all_tbl[seeing_key]).astype(float) ) < seeing_cut
        if (depth_cut is not None) & (seeing_cut is not None):
            select_index = (np.array(all_tbl[depth_key].data).astype(float)  > depth_cut) & (np.array(all_tbl[seeing_key].data).astype(float)  < seeing_cut)
        if (depth_cut is None) & (seeing_cut is None):
            select_index = ~(sigma_clip(list(all_tbl[depth_key]),sigma = sigma, maxiters= 5).mask) & ~(sigma_clip(list(all_tbl[seeing_key]),sigma = sigma, maxiters= 5).mask)
        cut_tbl = all_tbl[~select_index]
        selected_tbl = all_tbl[select_index]
        self.failed_imagelist['outlier'] = list(cut_tbl['file'])
        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_outlier_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_outlier_failed_imagelist.txt")
            
        if visualize:
            plt.figure(figsize = (10,5))
            plt.scatter(all_tbl[seeing_key],all_tbl[depth_key], c= 'k', marker = 'o', alpha = 0.6)
            plt.scatter(cut_tbl[seeing_key],cut_tbl[depth_key], c= 'r', label = 'Outlier', marker = 'o', alpha = 0.6)
            plt.xlabel('Seeing[arcsec]')
            plt.ylabel('Depth[AB]')
            plt.grid()
            plt.legend()
            plt.show()
        if (move == True) & (len(cut_tbl) > 0):
            for image in cut_tbl['file']:
                os.system(f'rm {image}')
        if (move == True) & (len(masked_tbl) > 0):
            for image in masked_tbl['file']:
                os.system(f'rm {image}')
        self.target_imagelist = selected_tbl['file']

    def combine(self,
                filelist: list = None,
                group_key: str = 'jd',
                group_tolerance: float = 0.5,
                align : bool = False, 
                align_detection_sigma : float = 5,
                cutout : bool = False,
                cutout_size : float = 0.95,
                zp_key: str = 'ZP_APER_1_HH', 
                **kwargs):
        
        if filelist is None:
            filelist = self.target_imagelist
        group_tbl = self.group_table(tbl=self.get_imginfo(filelist=filelist), key=group_key, tolerance=group_tolerance)
        tbl_groups = group_tbl.group_by('group').groups
        self.target_imagelist = list()
        for tbl_group in tqdm(tbl_groups, desc = 'Processing images...'):
            files_to_combine = list()
            tbl_group.sort('jd')
            for file_ in tbl_group['file']:
                failed = False
                im = Image(image=file_, telescope_info=self.telinfo, reference_image=tbl_group['file'][len(tbl_group['file']) // 2])
                if zp_key not in im.header.keys():
                    self.failed_imagelist['calculate'].append(file_)
                if align:
                    try:
                        im.align(cutout = cutout, cutout_size = cutout_size, detection_sigma = align_detection_sigma)
                    except:
                        self.failed_imagelist['align'].append(file_)
                        failed = True
                if cutout & (not align):
                    try:
                        im.cutout(size = cutout_size)
                    except:
                        self.failed_imagelist['cutout'].append(file_)
                        failed = True
                if not failed:
                    files_to_combine.append(im.target_image)
            files_to_combine = list(files_to_combine)
            if len(files_to_combine) > 0:
                try:
                    print(f'Combining {len(files_to_combine)} images in group {tbl_group["group"][0]}')
                    combined_image = self.combine_img(filelist=files_to_combine, zp_key=zp_key)#, **kwargs)
                    self.target_imagelist.append(combined_image)
                except:
                    self.failed_imagelist['combine'].append(file_)
        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_combine_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_combine_failed_imagelist.txt")

    def faster_combine(self, 
                       filelist : list = None,
                       group_key: str = 'jd',
                       group_tolerance: float = 0.5,
                       align : bool = False,
                       align_detection_sigma : float = 5,
                       cutout : bool = False,
                       cutout_size : float = 0.95,
                       num_processes: int = 4, 
                       zp_key: str = 'ZP_APER_1_HH',
                       **kwargs):
        
        if filelist is None:
            filelist = self.target_imagelist

        # Split the target image list into groups
        print(f"Querying image information...")
        imageinfo = self.get_imginfo(filelist=self.target_imagelist)
        print('Querying image information... Done')
        
        groupped_tbl = self.group_table(tbl=imageinfo, key=group_key, tolerance=group_tolerance).group_by('group')
        groups = groupped_tbl.groups
        
        num_images = len(self.target_imagelist)
        num_images_per_group = num_images // num_processes
        groups_idx = [int(groupped_tbl['group'][i * num_images_per_group]) for i in range(num_processes)]
        
        image_groups = [list(groups[groups_idx[i]:groups_idx[i+1]]['file']) for i in range(len(groups_idx) - 1)]
        image_groups.append(list(groups[groups_idx[-1]:]['file']))
                
        # Create a pool of processes and map the worker function to each group
        print('Combining images... with {} processes'.format(num_processes))
        with multiprocessing.Pool(num_processes) as pool:
            results = []
            for group in image_groups:
                args = (group, group_key, group_tolerance, align, align_detection_sigma, cutout, cutout_size, zp_key)
                result = pool.apply_async(self._combine_worker, args=args, kwds=kwargs)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]

    def _combine_worker(self, filelist, group_key='jd', group_tolerance=0.2, align = False, align_detection_sigma = 5, cutout = False, cutout_size = 0.95, zp_key='ZP_APER_1_HH', **kwargs):
        self.combine(filelist=filelist, 
                     group_key=group_key,
                     group_tolerance=group_tolerance, 
                     align=align,
                     align_detection_sigma=align_detection_sigma,
                     cutout=cutout,
                     cutout_size=cutout_size,
                     zp_key=zp_key,
                     **kwargs)
        
    def photometry(self, 
                   ra, 
                   dec, 
                   filelist : list = None,
                   sex_configfile: str = None,
                   detect_threshold : float  = 3.0,
                   aperture_type : str = 'relative', # relative or absolute
                   aperture_sizes : list = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                   trim_target_image : bool = True,
                   trim_reference_image : bool = False,
                   trim_size : int = 1500,
                   subtract : bool = True,
                   visualize : bool = True):
        
        if filelist is None:
            filelist = self.target_imagelist
            
        for file_ in filelist:
            im = Image(image=file_, telescope_info=self.telinfo, reference_image=self.reference_image)
            try:
                im.photometry(ra = ra,
                                dec = dec,
                                sex_configfile = sex_configfile, 
                                detect_threshold = detect_threshold,
                                aperture_type = aperture_type,
                                aperture_sizes = aperture_sizes,
                                trim_target_image = trim_target_image,
                                trim_reference_image = trim_reference_image,
                                trim_size = trim_size,
                                subtract = subtract,
                                visualize = visualize)
            except:
                self.failed_imagelist['photometry'].append(file_)
        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_photometry_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_photometry_failed_imagelist.txt")

    def faster_photometry(self, 
                          ra, 
                          dec, 
                          filelist : list = None,
                          sex_configfile: str = None,
                          detect_threshold : float  = 3.0,
                          aperture_type : str = 'relative', # relative or absolute
                          aperture_sizes : list = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                          trim_target_image : bool = True,
                          trim_reference_image : bool = False,
                          trim_size : int = 1500,
                          subtract : bool = True,
                          visualize : bool = True,
                          num_processes : int = 4
                          ):

        # Split the target image list into groups
        print(f"Splitting image list into {num_processes} groups...")
        num_images = len(self.target_imagelist)
        num_images_per_group = num_images // num_processes
        image_groups = [self.target_imagelist[i * num_images_per_group:(i + 1) * num_images_per_group] for i in range(num_processes - 1)]
        image_groups.append(self.target_imagelist[(num_processes - 1) * num_images_per_group:])
        
        # Create a pool of processes and map the worker function to each group
        print(f'Running photometry with {num_processes} processes')
        with multiprocessing.Pool(num_processes) as pool:
            results = []
            for group in image_groups:
                args = (ra, dec, group, sex_configfile, detect_threshold, aperture_type, aperture_sizes, trim_target_image, trim_reference_image, trim_size, subtract, visualize)
                result = pool.apply_async(self._worker_photometry, args=args)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]
        return output

    def _worker_photometry(self, ra, dec, filelist, sex_configfile, detect_threshold, aperture_type, aperture_sizes, trim_target_image, trim_reference_image, trim_size, subtract, visualize):
        self.photometry(ra = ra, 
                        dec = dec, 
                        filelist = filelist,
                        sex_configfile = sex_configfile,
                        detect_threshold = detect_threshold,
                        aperture_type = aperture_type,
                        aperture_sizes = aperture_sizes,
                        trim_target_image = trim_target_image,
                        trim_reference_image = trim_reference_image,
                        trim_size = trim_size,
                        subtract = subtract,
                        visualize = visualize)
 
#%% KCT
if __name__ == '__main__':    
    #filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/g/com*120.fits'))
    #A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'), reference_image = '/data1/reference_image/KCT_STX16803/Ref-KCT_STX16803-NGC1566-g-4440.com.fits')
    #filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/r/com*120.fits'))
    #A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'), reference_image = '/data1/reference_image/KCT_STX16803/Ref-KCT_STX16803-NGC1566-r-3360.com.fits')
    filelist = sorted(glob.glob('/data2/supernova_rawdata/SN2021aefx/photometry/KCT_STX16803/i/com*120.fits'))
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'), reference_image = '/data1/reference_image/KCT_STX16803/Ref-KCT_STX16803-NGC1566-i-720.com.fits')
    sex_configfile = None
#%% LSGT
if __name__ == '__main__':    
    # filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/LSGT_SNUCAMII/g/com*180.fits'))
    # A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'LSGT', ccd = 'SNUCAMII'), reference_image = '/data1/reference_image/LSGT_STX16803/Calib-LSGT-NGC1566-20210916-180448-g-540.com.fits')
    # filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/LSGT_SNUCAMII/r/com*180.fits'))
    # A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'LSGT', ccd = 'SNUCAMII'), reference_image = '/data1/reference_image/LSGT_STX16803/Calib-LSGT-NGC1566-20210916-181452-r-540.com.fits')
    filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/LSGT_SNUCAMII/i/com*180.fits'))
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'LSGT', ccd = 'SNUCAMII'), reference_image = '/data1/reference_image/LSGT_STX16803/Calib-LSGT-NGC1566-20220401-100321-i-540.com.fits')
    sex_configfile = None
#%% RASA36
if __name__ == '__main__':    
    filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/RASA36/r/HIGH/*60.fits'))
    filelist_1 = filelist[:1209]
    filelist_2 = filelist[1209:]
#%% RASA36, continued...
if __name__ == '__main__':    
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'RASA36', ccd = 'KL4040', readoutmode = 'HIGH'), reference_image = '/data1/reference_image/RASA36_KL4040/Ref-RASA36-NGC1566-r-3180-HIGH.com.fits')
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_HIGH.config'
#%% RASA36, continued...
if __name__ == '__main__':    
    filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/RASA36/r/MERGE/com*60.fits'))
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'RASA36', ccd = 'KL4040', readoutmode = 'MERGE'), reference_image = '/data1/reference_image/RASA36_KL4040/Ref-RASA36-NGC1566-r-3180-MERGE.com.fits')
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_MERGE.config'
    
#%% Calculate ZP
if __name__ == '__main__':    
    #A.calculate(sex_configfile = sex_configfile)
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 3, ref_catalog_name = 'APASS', ref_catalog_conversion ='PS1')
#%%
if __name__ == '__main__':    
    A.exclude_outlier(sigma = 3)
#%%
if __name__ == '__main__':    
    #A.combine(group_tolerance= 0.3)
    A.faster_combine(group_tolerance= 0.3, num_processes= 3)
#%%
if __name__ == '__main__':    
    #A.photometry(ra=64.9723704, dec=-54.9481347, sex_configfile = sex_configfile, calculate = True, trim_reference_image= True, trim_target_image= True, subtract = True, visualize= True)
    A.faster_photometry(ra = 64.9723704, 
                        dec = -54.9481347, 
                        filelist = filelist,
                        num_processes = 3,
                        sex_configfile = sex_configfile,
                        detect_threshold = 3.0,
                        aperture_type = 'relative', # relative or absolute
                        aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                        trim_target_image = False,
                        trim_reference_image = False,
                        trim_size = 2000,
                        subtract = True,
                        visualize = True)
# %%
if __name__ == '__main__':    
    filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/RASA36/reference_image/Calib*60.fits'))
#%%
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'RASA36', ccd = 'KL4040', readoutmode = 'HIGH'), reference_image = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/RASA36/reference_image/com_align_Calib-RASA36-NGC1097-20211130-040342-r-60.fits')
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_HIGH.config'
# %%
if __name__ == '__main__':
    A.astrometry(sex_configfile)
#%%
if __name__ == '__main__':
    A.exclude_outlier(sigma = 3, depth_cut = 19.3, seeing_cut = 3.5, depth_key = 'UL5_1', seeing_key = 'SEEING', move = True, visualize = True)
    #A.exclude_outlier(sigma = 3, depth_cut = 19.5, seeing_cut = 3.2, depth_key = 'UL5_1', seeing_key = 'SEEING', move = True, visualize = True)

#%%
if __name__ == '__main__':    
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 3, ref_catalog_name = 'APASS', ref_catalog_conversion = None)
#%%
if __name__ == '__main__':    
    A.combine(group_tolerance= 10000,  align = True, cutout = False, cutout_size = 0.95, align_detection_sigma = 3)
    #A.faster_combine(group_tolerance= 0.3, num_processes= 6,  align = True, cutout = True, cutout_size = 0.95, align_detection_sigma = 3)
# %%
if __name__ == '__main__':    
    #A.photometry(ra=64.9723704, dec=-54.9481347, sex_configfile = sex_configfile, calculate = True, trim_reference_image= True, trim_target_image= True, subtract = True, visualize= True)
    A.photometry(ra = 41.575542, 
                dec = -30.239489, 
                #filelist = filelist,
                sex_configfile = sex_configfile,
                detect_threshold = 1.5,
                aperture_type = 'relative', # relative or absolute
                aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                trim_target_image = True,
                trim_reference_image = True,
                trim_size = 2000,
                subtract = True,
                visualize = True)
# %%


# %%
if __name__ == '__main__':    
    filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2023rve/analysis/LSGT_ASI16000/reference_image/*-r-180.fits'))
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'LSGT', ccd = 'ASI1600MM'), reference_image = None)
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_HIGH.config'
# %%
if __name__ == '__main__':
    A.exclude_outlier(sigma = 3, depth_key = 'DEPTH5_APER_1_HH', seeing_key = 'SEEING_HH', move = False, visualize = True)

#%%
if __name__ == '__main__':    
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 6, ref_catalog_name = 'APASS', ref_catalog_conversion = None)
#%%
if __name__ == '__main__':    
    A.faster_combine(group_tolerance= 0.3, num_processes= 3,  align = True, cutout = True, cutout_size = 0.95, align_detection_sigma = 3)
# %%
if __name__ == '__main__':    
    #A.photometry(ra=64.9723704, dec=-54.9481347, sex_configfile = sex_configfile, calculate = True, trim_reference_image= True, trim_target_image= True, subtract = True, visualize= True)
    A.faster_photometry(ra = 41.575542, 
                        dec = -30.239489, 
                        filelist = filelist,
                        num_processes = 6,
                        sex_configfile = sex_configfile,
                        detect_threshold = 3.0,
                        aperture_type = 'relative', # relative or absolute
                        aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                        trim_target_image = True,
                        trim_reference_image = True,
                        trim_size = 2000,
                        subtract = True,
                        visualize = True)
# %%
