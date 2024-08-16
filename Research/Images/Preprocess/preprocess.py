#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 20:51:01 2022

@author: hhchoi1022
"""

###### Preprocessing ######


# %%
from astropy.io import fits
import ccdproc
import os
import glob
from ccdproc import CCDData
from ccdproc import Combiner
import numpy as np
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astropy.wcs import FITSFixedWarning
import warnings
from ccdproc import ImageFileCollection
from HHsupport_phot import combine_img
from astropy.table import Table
warnings.simplefilter('ignore', category=FITSFixedWarning)
# %%

class Prepcoess:
    
    def __init__(self, 
                 imdir,
                 calibdir = None
                ):
        
        self.imdir = imdir
        self.fitsinfo = self.get_info_images(self.imdir)
        self.calibdir = calibdir
        self.calibinfo = Table()
        if not calibdir == None:
            self.calibinfo = self.get_info_images(self.calibdir)
        
    def get_info_images(self, directory):
        '''
        parameters
        ----------
        {the directory absolute path} containing fits images
        
        returns 
        -------
        fits file information in astropy table format
        
        notes 
        -----
        -----
        '''
        iminfo = ImageFileCollection(directory).summary
        absfiles = []
        for file in iminfo['file']:
            absfiles.append(directory + file)
        iminfo['file'] = absfiles
        return iminfo
    
    def master_bias(self, bias_collection,
                    clip = 'sigma',
                    combine = 'median',
                    write = True,
                    c_filename = 'Mbias',
                    
                    # Clipping
                    clip_sigma_low = 2,
                    clip_sigma_high = 5,
                    clip_minmax_min = 3,
                    clip_minmax_max = 3,
                    clip_extrema_nlow = 1,
                    clip_extrema_nhigh = 1
                    ):
        biasfiles = bias_collection['file']
        
        mbias = combine_img(biasfiles,
                            clip = clip,
                            combine = combine,
                            write = write,
                            c_filename = self.imdir+c_filename+'.fits',
                            clip_sigma_low = clip_sigma_low,
                            clip_sigma_high = clip_sigma_high,
                            clip_minmax_min = clip_minmax_min,
                            clip_minmax_max = clip_minmax_max,
                            clip_extrema_nlow = clip_extrema_nlow,
                            clip_extrema_nhigh = clip_extrema_nhigh
                            )
        return mbias
    
    def master_dark(self, dark_collection, mbias_key = 'Mbias',
                    clip = 'sigma',
                    combine = 'median',
                    write = True,
                    c_filename = 'Mdark',
                    
                    # Clipping
                    clip_sigma_low = 2,
                    clip_sigma_high = 5,
                    clip_minmax_min = 3,
                    clip_minmax_max = 3,
                    clip_extrema_nlow = 1,
                    clip_extrema_nhigh = 1
                    ):
        mbias = self.imdir + mbias_key + '.fits'
        mbias_data = CCDData.read(mbias, unit='adu')
        groups_exp = dark_collection.group_by('exptime')
        list_mdark_data = []
        list_mdark_hdr = []
        for group_exp in groups_exp.groups:
            exp_filelist = group_exp['file']
            exptime = group_exp[0]['exptime']
            exp_datalist = []
            for file in exp_filelist:
                dark = CCDData.read(file, unit='adu')
                dark.header['MBIAS'] = mbias
                b_subtracted = ccdproc.subtract_bias(dark, mbias_data)
                exp_datalist.append(b_subtracted)
            mdark_exp = combine_img(exp_datalist, 
                                    clip = clip,
                                    combine = combine,
                                    write = True,
                                    c_filename = f'{self.imdir}{c_filename}{int(exptime)}.fits',
                                    clip_sigma_low = clip_sigma_low,
                                    clip_sigma_high = clip_sigma_high,
                                    clip_minmax_min = clip_minmax_min,
                                    clip_minmax_max = clip_minmax_max,
                                    clip_extrema_nlow = clip_extrema_nlow,
                                    clip_extrema_nhigh = clip_extrema_nhigh
                                    )
            list_mdark_data.append(mdark_exp)
            list_mdark_hdr.append(mdark_exp.header)
        return list_mdark_data, list_mdark_hdr
    
    def master_flat(self, flat_collection, mbias_key = 'Mbias', mdark_key = 'Mdark',
                    clip = 'sigma',
                    combine = 'median',
                    write = True,
                    c_filename = 'Mflat',
                    scale = True,
                    scale_mode = 'multiply',
                    
                    # Clipping
                    clip_sigma_low = 2,
                    clip_sigma_high = 5,
                    clip_minmax_min = 3,
                    clip_minmax_max = 3,
                    clip_extrema_nlow = 1,
                    clip_extrema_nhigh = 1
                    ):
        mbias = self.imdir + mbias_key + '.fits'
        mdarklist = sorted(glob.glob(self.imdir + mdark_key + '*.fits'))
        #mdark = mdarklist[-1]
        mdark = mdarklist[0]
        mbias_data = CCDData.read(mbias, unit='adu')
        mdark_data = CCDData.read(mdark, unit='adu')
        
        groups_filter = flat_collection.group_by('filter')
        list_flat_data = []
        list_flat_hdr = []
        for group_filter in groups_filter.groups:
            filter_filelist = group_filter['file']
            filter_ = group_filter[0]['filter']
            filter_datalist = []
            for file in filter_filelist:
                flat = CCDData.read(file, unit='adu')
                b_subtracted = ccdproc.subtract_bias(flat, mbias_data)
                db_subtracted = ccdproc.subtract_dark(b_subtracted, mdark_data,
                                                      exposure_time='EXPTIME',
                                                      exposure_unit=u.second,
                                                      scale=True)
                filter_datalist.append(db_subtracted)
            mflat_filter = combine_img(filter_datalist, 
                                        clip = clip,
                                        combine = combine,
                                        write = write,  
                                        scale = scale,
                                        scale_mode = scale_mode,
                                        c_filename = f'{self.imdir}{c_filename}{filter_}.fits',
                                        clip_sigma_low = clip_sigma_low,
                                        clip_sigma_high = clip_sigma_high,
                                        clip_minmax_min = clip_minmax_min,
                                        clip_minmax_max = clip_minmax_max,
                                        clip_extrema_nlow = clip_extrema_nlow,
                                        clip_extrema_nhigh = clip_extrema_nhigh
                                        )
            list_flat_data.append(mflat_filter)
            list_flat_hdr.append(mflat_filter.header)
        return list_flat_data, list_flat_hdr

    def cor_img(self, obj_collection, mbias_key = 'Mbias', mdark_key = 'Mdark', mflat_key = 'Mflat',
                suffix='Calib_',
                ):
        mbias_data = CCDData.read(f'{self.imdir}{mbias_key}.fits', unit='adu')
        groups_filter = obj_collection.group_by('filter')
        for group_filter in groups_filter.groups:
            filter_filelist = group_filter['file']
            filter_ = str(group_filter[0]['filter'])
            mdarklist = sorted(glob.glob(self.imdir + mdark_key + '*.fits'))
            #mdark_name = mdarklist[1]
            mdark_name = mdarklist[0]
            mflat_name = f'{self.imdir}{mflat_key}{filter_}.fits'
            mdark_data = CCDData.read(mdark_name, unit='adu')
            mflat_data = CCDData.read(mflat_name, unit='adu')
            for file in filter_filelist:
                obj_data = CCDData.read(file, unit = 'adu')
                b_subtracted = ccdproc.subtract_bias(obj_data, mbias_data)
                db_subtracted = ccdproc.subtract_dark(b_subtracted, mdark_data,
                                                     exposure_time='EXPTIME',
                                                     exposure_unit=u.second,
                                                     scale=True)
                fdb_subtracted = ccdproc.flat_correct(db_subtracted, mflat_data)
                fdb_subtracted.header = obj_data.header
                fdb_subtracted.header['MBIAS'] = f'{self.imdir}{mbias_key}.fits'
                fdb_subtracted.header['MDARK'] = f'{mdark_name}'
                fdb_subtracted.header['MFLAT'] = f'{self.imdir}{mflat_key}{filter_}.fits'
                fdb_subtracted.write(f'{self.imdir}{suffix}{os.path.basename(file)}', overwrite=True)
    
    def main(self):
        if 'ncombine' in self.fitsinfo.keys():
            self.fitsinfo.remove_rows(self.fitsinfo['ncombine'].mask == False)
        
        Mbias = Table()
        Mdark = Table()
        Mflat = Table()
        if len(self.calibinfo) > 0:
            Mbias = self.calibinfo[(self.calibinfo['imagetyp'] == 'zero') | (self.calibinfo['imagetyp'] == 'bias') | (self.calibinfo['imagetyp'] == 'Bias Frame')]
            Mdark = self.calibinfo[(self.calibinfo['imagetyp'] == 'dark') | (self.calibinfo['imagetyp'] == 'Dark Frame')]
            Mflat = self.calibinfo[(self.calibinfo['imagetyp'] == 'flat') | (self.calibinfo['imagetyp'] == 'Flat Field')]

        # Bias
        if len(Mbias) != 0:
            os.system(f'cp {Mbias["file"][0]} {self.imdir}/Mbias.fits')
        else:        
            bias_collection = self.fitsinfo[(self.fitsinfo['imagetyp'] == 'zero') | (self.fitsinfo['imagetyp'] == 'bias') | (self.fitsinfo['imagetyp'] == 'Bias Frame')]
            self.master_bias(bias_collection = bias_collection, c_filename = 'Mbias')
        
        # Dark
        if len(Mdark) != 0:
            os.system(f'cp {Mdark["file"][0]} {self.imdir}/Mdark.fits')
        else:
            dark_collection = self.fitsinfo[(self.fitsinfo['imagetyp'] == 'dark') | (self.fitsinfo['imagetyp'] == 'Dark Frame')]
            self.master_dark(dark_collection = dark_collection, mbias_key = 'Mbias')
        
        # Flat
        if len(Mflat) != 0:
            for mflatfile in Mflat:
                filter_ = mflatfile['filter']
                os.system(f'cp {mflatfile["file"]} {self.imdir}/Mflat{filter_}.fits')
        else:
            flat_collection = self.fitsinfo[(self.fitsinfo['imagetyp'] == 'flat') | (self.fitsinfo['imagetyp'] == 'Flat Field')]
            exptime = flat_collection['exptime'][0]
            self.master_flat(flat_collection = flat_collection, mbias_key = 'Mbias', mdark_key = 'Mdark')
        
        # Correction
        obj_collection = self.fitsinfo[(self.fitsinfo['imagetyp'] == 'object') | (self.fitsinfo['imagetyp'] == 'Light Frame')]
        self.cor_img(obj_collection, mbias_key = f'Mbias', mdark_key = 'Mdark', mflat_key = 'Mflat', suffix = 'Calib_')
        

# %%
