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
        self.failed_imagelist['astrometry'] = list()
        self.failed_imagelist['calculate'] = list()
        self.failed_imagelist['cutout'] = list()
        self.failed_imagelist['align'] = list()
        self.failed_imagelist['combine'] = list()
        self.failed_imagelist['photometry'] = list()

    def scamp(self,
              filelist : list = None,
              sex_configfile : str = None,
              scamp_configfile : str = None,
              print_progress : bool = True):
        
        if filelist is None:
            filelist = self.target_imagelist
        
        if sex_configfile == None:
            sex_configfile = f"{os.path.join(self.sexpath,self.telinfo['obs'])}.scampconfig"
            
        result = self.run_scamp(filelist = filelist,
                                sex_configfile = sex_configfile,
                                scamp_configfile = scamp_configfile,
                                update_files = True,
                                print_progress = print_progress
                                )
        
    def faster_scamp(self):
        pass

    ################# START Astrometry #################
    def astrometry(self, 
                   filelist : list = None,
                   sex_configfile : str = None,
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
                              ra = ra,
                              dec = dec,
                              radius = radius,
                              scalelow = scalelow,
                              scalehigh = scalehigh,
                              overwrite = overwrite,
                              remove = remove
                              )
            except:
                self.failed_imagelist['astrometry'].append(file_)
        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_astrometry_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_astrometry_failed_imagelist.txt")
        
        if remove:
            image_paths = set([os.path.dirname(file) for file in filelist])
            for path in image_paths:
                os.chdir(path)
                os.system(f'rm tmp* astrometry* *.conv default.nnw *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match check.fits *.param {os.path.basename(sex_configfile)}')


    def faster_astrometry(self,
                          filelist : list = None,
                          sex_configfile : str = None,
                          num_processes : int = 3,
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
        image_paths = set([os.path.dirname(file) for file in filelist])
            
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
                kwargs = dict(filelist=group,
                              sex_configfile=sex_configfile,
                              ra=ra,
                              dec=dec,
                              radius=radius,
                              scalelow=scalelow,
                              scalehigh=scalehigh,
                              overwrite=overwrite)
                result = pool.apply_async(self._worker_astrometry, kwds=kwargs)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]
        
        if remove:
            for path in image_paths:
                os.chdir(path)
                os.system(f'rm tmp* astrometry* *.conv default.nnw *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match check.fits *.param {os.path.basename(sex_configfile)}')

        return output
    
    def _worker_astrometry(self, filelist, sex_configfile, ra, dec, radius, scalelow, scalehigh, overwrite):
        self.astrometry(filelist=filelist,
                        sex_configfile=sex_configfile,
                        ra=ra,
                        dec=dec,
                        radius=radius,
                        scalelow=scalelow,
                        scalehigh=scalehigh,
                        overwrite=overwrite,
                        remove=False)

    ################# END Astrometry #################
    
    ################# START Calculate ZP #################
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
                kwargs = dict(filelist=group, 
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
                result = pool.apply_async(self._worker_calculate, kwds=kwargs)
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
    ################# END Calculate ZP #################
    
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
                os.makedirs(os.path.join(os.path.dirname(image), "outlier"), exist_ok=True)
                os.system(f'mv {image} {os.path.join(os.path.dirname(image), "outlier")}')
        if (move == True) & (len(masked_tbl) > 0):
            for image in masked_tbl['file']:
                os.makedirs(os.path.join(os.path.dirname(image), "masked"), exist_ok=True)
                os.system(f'mv {image} {os.path.join(os.path.dirname(image), "masked")}')
        self.target_imagelist = selected_tbl['file']

    ################# START Combine #################

    def combine(self,
                filelist: list = None,
                group_key: str = 'jd',
                group_tolerance: float = 0.5,
                align : bool = False, 
                align_cut_outer : bool = False,
                align_outer_size : float = 0.95,
                align_detection_sigma : float = 5,
                zp_key: str = 'ZP_APER_1_HH'):
        
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
                        im.align(cut_outer = align_cut_outer, outer_size = align_outer_size, detection_sigma = align_detection_sigma)
                    except:
                        self.failed_imagelist['align'].append(file_)
                        failed = True
                if not failed:
                    files_to_combine.append(im.target_image)
            files_to_combine = files_to_combine
            if len(files_to_combine) > 0:
                try:
                    print(f'Combining {len(files_to_combine)} images in group {tbl_group["group"][0]}')
                    combined_image = self.combine_img(filelist=files_to_combine, zp_key=zp_key)
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
                       align_cut_outer : bool = False,
                       align_outer_size : float = 0.95,
                       align_detection_sigma : float = 5,
                       num_processes: int = 4, 
                       zp_key: str = 'ZP_APER_1_HH'):
        
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
                kwargs = dict(filelist=group,
                              group_key=group_key, 
                              group_tolerance=group_tolerance, 
                              align=align, 
                              align_cut_outer = align_cut_outer,
                              align_outer_size = align_outer_size,
                              align_detection_sigma = align_detection_sigma, 
                              zp_key=zp_key)
                result = pool.apply_async(self._worker_combine, kwds=kwargs)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]

    def _worker_combine(self, filelist, group_key, group_tolerance, align, align_cut_outer, align_outer_size, align_detection_sigma, zp_key):
        self.combine(filelist=filelist,
                     group_key=group_key, 
                     group_tolerance=group_tolerance, 
                     align=align, 
                     align_cut_outer = align_cut_outer,
                     align_outer_size = align_outer_size,
                     align_detection_sigma = align_detection_sigma, 
                     zp_key=zp_key)

    ################# END Combine #################

    def cutout(self,
               filelist : list = None,
               cutout_size : int = 2000):
        if filelist is None:
            filelist = self.target_imagelist
            
        for file_ in filelist:
            im = Image(image=file_, telescope_info=self.telinfo, reference_image=self.reference_image)
            try:
                im.cutout(cutout_size = cutout_size)
            except:
                self.failed_imagelist['cutout'].append(file_)

        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_cutout_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_cutout_failed_imagelist.txt")

    ################# START Align #################
    
    def align(self, 
              filelist : list = None,
              cut_outer : bool = False,
              outer_size : float = 0.95,
              detection_sigma : float = 5):
        
        if filelist is None:
            filelist = self.target_imagelist
            
        for file_ in filelist:
            im = Image(image=file_, telescope_info=self.telinfo, reference_image=self.reference_image)
            try:
                im.align(cut_outer = cut_outer, outer_size = outer_size, detection_sigma = detection_sigma)
            except:
                self.failed_imagelist['align'].append(file_)

        # Save
        save_path = os.path.join(self.photpath, 'photometry_log',datetime.now().strftime('%Y%m%d'))
        os.makedirs(save_path, exist_ok=True)
        with open(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_align_failed_imagelist.txt", 'w') as file:
            json.dump(self.failed_imagelist, file, indent=4)
            print(f"{save_path}/{datetime.now().strftime('%Y%m%d_%H%M%S')}_align_failed_imagelist.txt")

        
    def faster_align(self,
                     filelist : list = None,
                     cut_outer : bool = False,
                     outer_size : float = 0.95,
                     detection_sigma : float = 5,
                     num_processes : int = 4):
        
        if filelist is None:
            filelist = self.target_imagelist

        # Split the target image list into groups
        print(f"Splitting image list into {num_processes} groups...")
        num_images = len(self.target_imagelist)
        num_images_per_group = num_images // num_processes
        image_groups = [self.target_imagelist[i * num_images_per_group:(i + 1) * num_images_per_group] for i in range(num_processes - 1)]
        image_groups.append(self.target_imagelist[(num_processes - 1) * num_images_per_group:])
        
        # Create a pool of processes and map the worker function to each group
        print('Image alignment with {} processes'.format(num_processes))
        with multiprocessing.Pool(num_processes) as pool:
            results = []
            for group in image_groups:
                kwargs = dict(filelist = group,
                              cut_outer = cut_outer,
                              outer_size = outer_size,
                              detection_sigma = detection_sigma)
                result = pool.apply_async(self._worker_align, kwds=kwargs)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]
        return output

    def _worker_align(self, filelist, cut_outer, outer_size, detection_sigma):
        self.align(filelist = filelist, 
                   cut_outer = cut_outer, 
                   outer_size = outer_size, 
                   detection_sigma = detection_sigma)
        
    ################# END Align #################
    
    ################# START Photometry #################
    def photometry(self, 
                   ra, 
                   dec, 
                   filelist : list = None,
                   sex_configfile: str = None,
                   detect_threshold : float  = 3.0,
                   aperture_type : str = 'relative', # relative or absolute
                   aperture_sizes : list = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                   cutout_target_image : bool = True,
                   cutout_reference_image : bool = False,
                   cutout_size : int = 1500,
                   align : bool = True,
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
                                cutout_target_image = cutout_target_image,
                                cutout_reference_image = cutout_reference_image,
                                cutout_size = cutout_size,
                                align = align,
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
                          cutout_target_image : bool = True,
                          cutout_reference_image : bool = False,
                          cutout_size : int = 1500,
                          align : bool = True,
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
                kwargs = dict(ra = ra, 
                              dec = dec, 
                              filelist = group,
                              sex_configfile = sex_configfile,
                              detect_threshold = detect_threshold,
                              aperture_type = aperture_type,
                              aperture_sizes = aperture_sizes,
                              cutout_target_image = cutout_target_image,
                              cutout_reference_image = cutout_reference_image,
                              cutout_size = cutout_size,
                              align = align,
                              subtract = subtract,
                              visualize = visualize)
                result = pool.apply_async(self._worker_photometry, kwds=kwargs)
                results.append(result)
            
            # Wait for all results to complete
            output = [result.get() for result in results]
        return output

    def _worker_photometry(self, ra, dec, filelist, sex_configfile, detect_threshold, aperture_type, aperture_sizes, cutout_target_image, cutout_reference_image, cutout_size, align, subtract, visualize):
        self.photometry(ra = ra, 
                        dec = dec, 
                        filelist = filelist,
                        sex_configfile = sex_configfile,
                        detect_threshold = detect_threshold,
                        aperture_type = aperture_type,
                        aperture_sizes = aperture_sizes,
                        cutout_target_image = cutout_target_image,
                        cutout_reference_image = cutout_reference_image,
                        cutout_size = cutout_size,
                        align = align,
                        subtract = subtract,
                        visualize = visualize)
    ################# END Photometry #################
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
    filelist = sorted(glob.glob('/data1/supernova_rawdata/SN2021aefx/photometry/RASA36/r/HIGH/Calib*60.fits'))
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
    #A.photometry(ra=64.9723704, dec=-54.9481347, sex_configfile = sex_configfile, calculate = True, cutout_reference_image= True, cutout_target_image= True, subtract = True, visualize= True)
    A.faster_photometry(ra = 64.9723704, 
                        dec = -54.9481347, 
                        filelist = filelist,
                        num_processes = 3,
                        sex_configfile = sex_configfile,
                        detect_threshold = 3.0,
                        aperture_type = 'relative', # relative or absolute
                        aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                        cutout_target_image = False,
                        cutout_reference_image = False,
                        cutout_size = 2000,
                        subtract = True,
                        visualize = True)
# %% #########################################################ls###

if __name__ == '__main__':
    filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/g/Calib*120.fits'))
if __name__ == '__main__':
    #A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'), reference_image = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/reference_image/Ref-KCT_STX16803-NGC1097-r-5400.com.fits')
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'KCT', ccd = 'STX16803'), reference_image = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/KCT_STX16803/reference_image/Ref-KCT_STX16803-NGC1097-g-4440.com.fits')
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/KCT.config'
#%%
if __name__ == '__main__':
    pass
    A.exclude_outlier(sigma = 2, depth_cut=17, seeing_cut=6, depth_key = 'DEPTH5_APER_1_HH', seeing_key = 'SEEING_HH', move = True, visualize = True)

#%%
if __name__ == '__main__':
    pass
    A.astrometry(sex_configfile  = sex_configfile, ra = 41.575542, dec = -30.239489, radius = 2)
    #A.faster_astrometry(sex_configfile  = sex_configfile, ra = 41.575542, dec = -30.239489, radius = 2, num_processes = 4)
#%%
if __name__ == '__main__':
    pass
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 6, ref_catalog_name = 'APASS', ref_catalog_conversion = None)
#%%
if __name__ == '__main__':
    A.combine(group_tolerance= 0.3,  align = True, align_cut_outer = True, align_detection_sigma = 3)
    #A.faster_combine(group_tolerance= 0.3, num_processes= 4,  align = True, align_cut_outer = False, align_outer_size = 0.95, align_detection_sigma = 3)

#%%
if __name__ == '__main__':
    A.align(cut_outer = False, outer_size = 0.95, detection_sigma = 5)
    #A.faster_align(cut_outer = False, outer_size = 0.95, detection_sigma = 5, num_processes =#%%
#%%
if __name__ == '__main__':    
    A.photometry(ra = 41.575542, 
                dec = -30.239489, 
                #filelist = filelist,
                sex_configfile = sex_configfile,
                detect_threshold = 3,
                aperture_type = 'relative', # relative or absolute
                aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                cutout_target_image = True,
                cutout_reference_image = True,
                cutout_size = 2000,
                align = False,
                subtract = True,
                visualize = True)


#%%
# (Reference image with no rotation, but shallower image)Combine -> Cutout with 2500 -> Align with cutouted reference image -> astrometry -> photometry
# (Reference image with rotation) Combine -> Cutout with 2500 -> Align with cutouted reference image -> Cutout with 2000 -> astrometry -> photometry

if __name__ == '__main__':    
    filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/RASA36/r/align_com*fits'))
    #filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/RASA36/reference_image_tmp/acom*fits'))
    filelist_1 = filelist[:1056]
    filelist_2 = filelist[1056:]
#%%
if __name__ == '__main__':
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'RASA36', ccd = 'KL4040', readoutmode = 'HIGH'), reference_image = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/RASA36/reference_image/acom_align_Calib-RASA36-NGC1097-20230811-072108-r-60.fits')
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_HIGH.config'
# %%
if __name__ == '__main__':
    #A.cutout(cutout_size = 2500)
    A.cutout(cutout_size = 2500)

#%%
if __name__ == '__main__':
    A.align(cut_outer = False, outer_size = 0.95, detection_sigma = 5)
    #A.faster_align(cut_outer = False, outer_size = 0.95, detection_sigma = 5, num_processes =#%%
#%%
if __name__ == '__main__':
    A.astrometry(sex_configfile  = sex_configfile, ra = 41.575542, dec = -30.239489, radius = 2)
    #A.faster_astrometry(sex_configfile = sex_configfile, num_processes = 2, ra = 41.575542, dec = -30.239489, radius = 2, overwrite = True)
    #A.scamp(sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/RASA36_HIGH.scampconfig')
    pass
#%%
if __name__ == '__main__':
    pass
    #A.exclude_outlier(sigma = 2, depth_key = 'DEPTH5_APER_1_HH', seeing_key = 'SEEING_HH', move = True, visualize = True)
    #A.exclude_outlier(sigma = 3, depth_key = 'UL5_1', seeing_key = 'SEEING', move = False, visualize = True)
    A.exclude_outlier(sigma = 3, depth_cut = 14, seeing_cut = 6, depth_key = 'DEPTH5_APER_1_HH', seeing_key = 'SEEING_HH', move = True, visualize = True)

#%%
if __name__ == '__main__':
    pass
    #A.calculate(sex_configfile = sex_configfile, ref_catalog_name = 'APASS', ref_catalog_conversion = None)
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 6, ref_catalog_name = 'APASS', ref_catalog_conversion = None)
#%%
if __name__ == '__main__':    
    #A.combine(group_tolerance= 10000,  align = True, align_cut_outer = False, align_detection_sigma = 5)
    A.faster_combine(group_tolerance= 0.3, num_processes= 4,  align = True, align_cut_outer = False, align_outer_size = 0.95, align_detection_sigma = 3)

# %%
if __name__ == '__main__':    
    A.photometry(ra = 41.575542, 
                dec = -30.239489, 
                #filelist = filelist,fin 
                sex_configfile = sex_configfile,
                detect_threshold = 3,
                aperture_type = 'relative', # relative or absolute
                aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                cutout_target_image = True,
                cutout_reference_image = True,
                cutout_size = 2000,
                align = False,
                subtract = True,
                visualize = True)
# %%
if __name__ == '__main__':    
    filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/LSGT_ASI1600MM/reference_image/r/com*-r-180.fits'))
    filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/analysis/LSGT_ASI1600MM/r/com*-r-180.fits'))
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = 'LSGT', ccd = 'ASI1600MM'), reference_image = '/mnt/data1/supernova_rawdata/SN2023rve/analysis/LSGT_ASI1600MM/reference_image/r/com_align_Calib-LSGT_ASI1600MM-NGC1097-20230623-181206-r-180.fits')
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/LSGT_ASI1600MM.config'
# %%
if __name__ == '__main__':
    A.exclude_outlier(sigma = 3, seeing_cut = 3.25, depth_cut = 19.4, depth_key = 'DEPTH5_APER_1_HH', seeing_key = 'SEEING_HH', move = False, visualize = True)

#%%
if __name__ == '__main__':    
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 3, ref_catalog_name = 'APASS', ref_catalog_conversion = None)
#%%
if __name__ == '__main__':    
    #A.combine(group_tolerance= 10000,  align = True, align_cut_outer = False, align_detection_sigma = 5)
    A.faster_combine(group_tolerance= 0.3, num_processes= 4,  align = True, align_cut_outer = False, align_outer_size = 0.95, align_detection_sigma = 3)

# %%
if __name__ == '__main__':    
    A.photometry(ra = 41.575542, 
                dec = -30.239489, 
                #filelist = filelist,fin 
                sex_configfile = sex_configfile,
                detect_threshold = 1.5,
                aperture_type = 'relative', # relative or absolute
                aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                cutout_target_image = False,
                cutout_reference_image = False,
                cutout_size = 2000,
                align = True,
                subtract = True,
                visualize = True)
# %%



#%%

# %%
if __name__ == '__main__':    
    filelist = sorted(glob.glob('/mnt/data1/supernova_rawdata/SN2023rve/7DT/late/r/com*.fits'))
    A = Photometry(filelist, telescope_info = PhotometryHelper().get_telinfo(telescope = '7DT', ccd = 'C361K', readoutmode = 'low'), reference_image = None)
    sex_configfile = '/home/hhchoi1022/hhpy/Research/photometry/sextractor/7DT_gain2750.config'
# %%
if __name__ == '__main__':
    A.exclude_outlier(sigma = 3, depth_key = 'DEPTH5_APER_1_HH', seeing_key = 'SEEING_HH', move = False, visualize = True)

#%%
if __name__ == '__main__':    
    A.faster_calculate(sex_configfile = sex_configfile, num_processes = 2, ref_catalog_name = 'APASS', ref_catalog_conversion = None)

#%%
if __name__ == '__main__':    
    A.combine(group_tolerance= 10000,  align = True, align_cut_outer = False, align_detection_sigma = 5)
   # A.faster_combine(group_tolerance= 0.3, num_processes= 4,  align = True, align_cut_outer = False, align_outer_size = 0.95, align_detection_sigma = 3)
# %%
if __name__ == '__main__':    
    #A.photometry(ra=64.9723704, dec=-54.9481347, sex_configfile = sex_configfile, calculate = True, cutout_reference_image= True, cutout_target_image= True, subtract = True, visualize= True)
    A.faster_photometry(ra = 41.575542, 
                        dec = -30.239489, 
                        filelist = filelist,
                        num_processes = 6,
                        sex_configfile = sex_configfile,
                        detect_threshold = 3.0,
                        aperture_type = 'relative', # relative or absolute
                        aperture_sizes = [1.5, 2.5, 3.5], # relative (1.5*seeing, 2.5*seeing, 3.5*seeing) or absolute (3", 5", 7")
                        cutout_target_image = True,
                        cutout_reference_image = True,
                        cutout_size = 2000,
                        subtract = True,
                        visualize = True)
# %%
