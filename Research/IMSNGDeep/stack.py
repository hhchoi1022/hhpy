#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from astropy.io import fits
import sep
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from photutils.aperture import CircularAperture
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import ascii
from regions.core import PixCoord
from regions.shapes.circle import CirclePixelRegion
#%%

from astropy.io import fits
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from photutils.background import SExtractorBackground

class Modelingbkg:
    
    """
    Find backgounnd modeld of an image
    The basic structure of this code is ASTRiDE
    
    ***Parameters
    ***
    ***
    ***
    
    
    """
    
    def __init__(self, filename, mode = 'constant', bkg_box_size = 50):
        
        hdulist = fits.open(filename)
        data = hdulist[0].data.astype(np.float64)
        hdr = hdulist[0].header
        
        # Check WCS information 
        try:
            wcsinfo = hdulist[0].header["CTYPE1"]
            if wcsinfo:
                self.wcs = WCS(hdr)
                self.filename = filename
        except:
            self.wcs = False
        
        hdulist.close()
        
        # Raw image.
        self.data = data
        self.hdr = hdr
        # Background structure and background map
        self._bkg = None
        self.background_map = None
        # Background removed image.
        self.image = None
        # Statistics for the image data.
        self._med = None
        self._std = None
        self.bkg_box_size = bkg_box_size
        
        bkg_options = ('constant', 'sex')
        if mode not in bkg_options:
            raise RuntimeError('"remove_bkg" must be the one among: %s' %
                               ', '.join(bkg_options))
        self.mode = mode
        
    def _bkgconst(self):
        # Get background map and subtract.
        sigma_clip = SigmaClip(sigma=3., maxiters=10)
        bkg_estimator = MedianBackground()
        self._bkg = Background2D(self.data,
                           (self.bkg_box_size, self.bkg_box_size),
                           filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        self.background_map = self._bkg.background
        self.image = self.data - self.background_map
    
        self._med = self._bkg.background_median
        self._std = self._bkg.background_rms_median

    def _bkgse(self):
        # Get background map and subtract.
        sigma_clip = SigmaClip(sigma=3., maxiters=10)
        bkg_estimator = SExtractorBackground(sigma_clip)
        self._bkg = Background2D(self.data,
                           (self.bkg_box_size, self.bkg_box_size),
                           filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        self.background_map = self._bkg.background
        self.image = self.data - self.background_map
    
        self._med = self._bkg.background_median
        self._std = self._bkg.background_rms_median
        
    def _plot(self, mode = 'img'):
        if self.wcs:
            if mode == 'img':
                plt.figure(figsize = (5,5), dpi = 100)
                plt.subplot(projection = self.wcs)
                m, s = np.median(self.data), np.std(self.data)
                plt.imshow(self.data, origin='lower', cmap='Greys_r', vmax = m+s, vmin = m-s, interpolation='nearest')
                plt.grid(color = 'grey', ls  ='solid', alpha = 0.6)
                plt.xlabel('RA')
                plt.ylabel('Dec')
            elif mode == 'bkg':
                try:
                    plt.figure(figsize = (5,5), dpi = 100)
                    plt.subplot(projection = self.wcs)
                    m, s = np.median(self.background_map), np.std(self.background_map)
                    plt.imshow(self.background_map, origin='lower', cmap='Greys_r', vmax = m+20*s, vmin = m-20*s, interpolation='nearest')
                    plt.grid(color = 'grey', ls  ='solid', alpha = 0.6)
                    plt.xlabel('RA')
                    plt.ylabel('Dec')
                except:
                    raise RuntimeError('Background map is not calculated. Run Modelingbkg._bkgconst(_bkgse) FIRST ')
            elif mode == 'sub':
                try:
                    plt.figure(figsize = (5,5), dpi = 100)
                    plt.subplot(projection = self.wcs)
                    m, s = np.median(self.image), np.std(self.image)
                    plt.imshow(self.image, origin='lower', cmap='Greys_r', vmax = m+1*s, vmin = m-1*s, interpolation='nearest')
                    plt.grid(color = 'grey', ls  ='solid', alpha = 0.6)
                    plt.xlabel('RA')
                    plt.ylabel('Dec')
                except:
                    raise RuntimeError('Background map is not calculated. Run Modelingbkg._bkgconst(_bkgse) FIRST ')
        else:
            if mode == 'img':
                plt.figure(figsize = (5,5), dpi = 100)
                plt.subplot(projection = self.wcs)
                m, s = np.median(self.data), np.std(self.data)
                plt.imshow(self.data, origin='lower', cmap='Greys_r', vmax = m+s, vmin = m-s, interpolation='nearest')
                plt.grid(color = 'grey', ls  ='solid', alpha = 0.6)
                plt.xlabel('RA')
                plt.ylabel('Dec')
            elif mode == 'bkg':
                try:
                    plt.figure(figsize = (5,5), dpi = 100)
                    plt.title('Background')
                    plt.subplot(projection = self.wcs)
                    m, s = np.median(self.background_map), np.std(self.background_map)
                    plt.imshow(self.background_map, origin='lower', cmap='Greys_r', vmax = m+10*s, vmin = m-10*s, interpolation='nearest')
                    plt.grid(color = 'grey', ls  ='solid', alpha = 0.6)
                    plt.xlabel('RA')
                    plt.ylabel('Dec')
                except:
                    raise RuntimeError('Background map is not calculated. Run Modelingbkg._bkgconst(_bkgse) FIRST ')
            elif mode == 'sub':
                try:
                    plt.figure(figsize = (5,5), dpi = 100)
                    plt.titme('bkg-Subtracted')
                    plt.subplot(projection = self.wcs)
                    m, s = np.median(self.image), np.std(self.image)
                    plt.imshow(self.image, origin='lower', cmap='Greys_r', vmax = m+1*s, vmin = m-1*s, interpolation='nearest')
                    plt.grid(color = 'grey', ls  ='solid', alpha = 0.6)
                    plt.xlabel('RA')
                    plt.ylabel('Dec')
                except:
                    raise RuntimeError('Background map is not calculated. Run Modelingbkg._bkgconst(_bkgse) FIRST ')

#%%

import glob, os, sys
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from astropy.io import fits
import astroalign as aa
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.io import fits
import astropy.units as u
from ccdproc import Combiner


class Stacking:
    
    
    """
    Stacking image to reduce noise
    
    """
    
    def __init__(self, imagelist, combinemode = 'iraf', check = True, remove = False, bkgsubtract = False):
        if type(imagelist) == str:
            imagelist = glob.glob(imagelist)
        
        
        self.imagelist = imagelist
        self._commode = combinemode
        self._checkmode = check
        self._removemode = remove
        self._bkgsubtractmode = bkgsubtract
        self.refim =None
        
    
    def show_outlier(self):

        seeinglist = []
        depthlist = []
        imagelist = []
        skyvallist = []
        skysiglist = []
        
        for image in self.imagelist:
            hdr = fits.getheader(image)
            try:
                depth = hdr['UL5_4']
                seeing = hdr['SEEING']
                skyval = hdr['SKYVAL']
                skysig = hdr['SKYSIG']
                
                seeinglist.append(seeing)
                depthlist.append(depth)
                imagelist.append(image)
                skyvallist.append(skyval)
                skysiglist.append(skysig)
            except:
                pass
        
        all_tbl = Table()
        all_tbl['Image'] = imagelist
        all_tbl['Depth'] = depthlist
        all_tbl['Seeing'] = seeinglist
        all_tbl['Skyval'] = skyvallist
        all_tbl['Skysig'] = skysiglist
        all_tbl = all_tbl[~(sigma_clip(all_tbl['Depth'],sigma = 3).mask) & (~sigma_clip(all_tbl['Seeing'],sigma = 3).mask) & (~sigma_clip(all_tbl['Skyval'],sigma = 3).mask)]
        
        
        if self._checkmode:
            plt.figure(figsize = (10,5), dpi = 150)
            plt.scatter(all_tbl['Seeing'],all_tbl['Depth'], marker = 'o', alpha = 0.6, c = all_tbl['Skyval'], cmap = 'coolwarm')
            plt.xlim(np.median(seeinglist)-1*np.std(seeinglist),np.median(seeinglist)+2*np.std(seeinglist))
            plt.ylim(np.median(depthlist)-2*np.std(depthlist),np.median(depthlist)+2*np.std(depthlist))
            plt.xlabel('Seeing[arcsec]')
            plt.ylabel('Depth[3*Seeing, 5 sigma]')
            plt.colorbar(location = 'top', label = 'Skyval')
            plt.grid()
            plt.legend()
            plt.show()
            
        self.imagelist = all_tbl['Image'].tolist()
        return all_tbl
    
    
    def select_images(self):
        all_tbl = self.show_outlier()
        choose2 = ''
        while choose2 not in ['Y', 'y']:
            cut_seeing = float(input('Give a cut off line for Seeing'))
            cut_depth= float(input('Give a cut off line for Depth'))
            cut_skyval = float(input('GIve a cut off line for Skyvalue'))
            cut_tbl = all_tbl[(all_tbl['Depth']>cut_depth) & (all_tbl['Seeing'] < cut_seeing) & (all_tbl['Skyval'] < cut_skyval)]
            cut_tbl.sort('Depth')
            counts = len(cut_tbl)
            meddepth = np.median(cut_tbl['Depth'])
            expdepth = meddepth+2.5*np.log10(np.sqrt(counts))
            medsigma = np.median(cut_tbl['Skysig'])
            expsigma = medsigma/np.sqrt(counts)
            medseeing = np.median(cut_tbl['Seeing'])
            expsurfacedepth = expdepth + 2.5*np.log10(3*medseeing)
            refimage = cut_tbl[int(counts/2)]['Image']
            plt.figure(figsize = (10,5), dpi = 150)
            plt.xlim(np.median(all_tbl['Seeing'])-1*np.std(all_tbl['Seeing']),np.median(all_tbl['Seeing'])+2*np.std(all_tbl['Seeing']))
            plt.ylim(np.median(all_tbl['Depth'])-2*np.std(all_tbl['Depth']),np.median(all_tbl['Depth'])+2*np.std(all_tbl['Depth']))
            plt.scatter(cut_tbl['Seeing'], cut_tbl['Depth'], marker = 'o', alpha = 0.6, label  =f'Selected[{len(cut_tbl)}]', c = 'r')
            plt.scatter(all_tbl['Seeing'], all_tbl['Depth'], marker = 'o', alpha = 0.1, c = 'k')
            plt.axvline(x= cut_seeing, c= 'r', linestyle = '--', linewidth = 1)
            plt.axhline(y= cut_depth, c= 'r', linestyle = '--', linewidth = 1)
            plt.xlabel('Seeing[arcsec]')
            plt.ylabel('Depth[3*Seeing, 5 sigma]')
            plt.text(cut_seeing, cut_depth-0.1, 'Expected depth = %.2f'%expdepth)
            plt.grid()
            plt.legend()
            plt.show()
            choose2 = input('Make a stacked image with %d images?'%counts)
        self.imagelist = cut_tbl['Image'].tolist()
        self.refim = refimage
        
    def scale_bkg(self, image):
        
        tgtim_hdu = fits.open(self.refim)
        image_hdu = fits.open(image)
        
        if ('SKYVAL' in image_hdu[0].header.keys()) & ('SKYVAL' in tgtim_hdu[0].header.keys()):
            bkg_offset = image_hdu[0].header['SKYVAL'] - tgtim_hdu[0].header['SKYVAL']
            
        else:
            # Referenced by Gregory Paek
            tgt_bkg = np.median(tgtim_hdu[0].data[~np.isnan(tgtim_hdu[0].data)].flatten())
            img_bkg = np.median(image_hdu[0].data[~np.isnan(image_hdu[0].data)].flatten())
            bkg_offset = tgt_bkg-img_bkg
            
        scaled_data = image_hdu[0].data-bkg_offset
        image_hdu[0].data = scaled_data
        
        return image_hdu
        
    def align_images(self, image, scale = True):
        
        ref_hdu = fits.open(self.refim)
        ref_data = ref_hdu[0].data
        ref_hdr = ref_hdu[0].header
        
        image_hdu = fits.open(image)
        image_data = image_hdu[0].data
        image_data = image_data.byteswap().newbyteorder()
        if scale == True:
            image_data = self.scale_bkg(image)[0].data
        image_hdr = image_hdu[0].header
        for key in ['JD','DATE-OBS','JD-OBS','MJD','SEEING','PEEING','APER_3','ZP_4','UL5_4']:
            if (key in ref_hdr.keys()) & (key in image_hdr.keys()):
                ref_hdr[key] = image_hdr[key]
        
        aligned_data, footprint = aa.register(image_data, ref_data, fill_value = 0)
        
        aligned_hdu = fits.PrimaryHDU(data = aligned_data, header = ref_hdr)
        aligned_hdu = fits.HDUList([aligned_hdu])
        imdir = os.path.dirname(image)
        os.makedirs(f'{imdir}/aligned', exist_ok=True)
        os.chdir(f'{imdir}/aligned')
        fits.writeto(f'{imdir}/aligned/aligned_{os.path.basename(image)}', aligned_hdu[0].data, aligned_hdu[0].header, overwrite = True)
        
        return imdir
    
    def combine_ccdproc(self, alignedlist):
        
        # Constructing data to combine
        comdataarray = []
        if type(alignedlist[0]) == str:
            for i, image in enumerate(alignedlist):
                hdu = fits.open(image)
                data = CCDData(hdu[0].data, unit = u.adu)
                comdataarray.append(data)
                self.printProgress(i, len(alignedlist), f'Stacking Progress[{i}]:', f'Complete[{len(alignedlist)}]', 1, 50)
        else:
            for i, hdu in enumerate(alignedlist):
                data = CCDData(hdu[0].data, unit = u.adu)
                comdataarray.append(data)
                self.printProgress(i, len(alignedlist), f'Stackikng Progress[{i}]:', f'Complete[{len(alignedlist)}]', 1, 50)
        hdr = hdu[0].header
        print('\nCombining...')
        combiner = Combiner(comdataarray)
        print('\nAlmost done...')
        if len(comdataarray) > 3:
            combiner.clip_extrema(nhigh=2, nlow = 1)
        hdu = combiner.median_combine()
        
        # Constructing header

        hdr['NCOMBINE'] = len(alignedlist)
        for n, inim in enumerate(alignedlist):
            hdr[f'COMB{n}'] = os.path.basename(inim)

        combined_hdu = fits.PrimaryHDU(data = hdu.data, header = hdr)
        combined_hdu = fits.HDUList([combined_hdu])
        
        # Writing combined output
        fits.writeto(f'{hdr["OBJECT"]}.com.fits', combined_hdu[0].data, combined_hdu[0].header, overwrite=True)
        
    def combine_iraf(self, alignedlist):
        import os
        os.chdir('/data2/iraf')
        from pyraf import iraf
        
        hdu = fits.open(alignedlist[0])
        hdr = hdu[0].header
        imdir = os.path.dirname(alignedlist[0])
        n_im = len(alignedlist)
        
        os.chdir(imdir)
        iraf.chdir(imdir)
        with open('image.list','w',encoding='UTF-8') as f:
            for name in alignedlist:
                f.write(name+'\n')
                
        if not n_im == 0:
            if n_im > 30:
                iraf.imcombine('@image.list', f'{hdr["OBJECT"]}.com.fits', combine = 'median', reject = 'minmax', zero = 'mode',nlow = 2, nhigh= 4)
            elif n_im > 15:
                iraf.imcombine('@image.list', f'{hdr["OBJECT"]}.com.fits', combine = 'median', reject = 'minmax', zero = 'mode',nlow = 1, nhigh= 3)
            else:
                iraf.imcombine('@image.list', f'{hdr["OBJECT"]}.com.fits', combine = 'median', reject = 'minmax', zero = 'mode', nhigh= 1)

            os.system('rm image.list') 
    
    def printProgress (self, iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100): 
        formatStr = "{0:." + str(decimals) + "f}" 
        percent = formatStr.format(100 * (iteration / float(total))) 
        filledLength = int(round(barLength * iteration / float(total))) 
        bar = '#' * filledLength + '-' * (barLength - filledLength) 
        sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)), 
        if iteration == total: 
            sys.stdout.write('\n') 
        sys.stdout.flush() 
        
    def run(self):
        self.select_images()
        imdir = os.path.dirname(self.imagelist[0])
        for i, image in enumerate(self.imagelist):
            imdir = self.align_images(image, scale = True)
            self.printProgress(i, len(self.imagelist), f'Alignment Progress[{i}]:', f'Complete[{len(self.imagelist)}]', 1, 50)
        alignedlist = glob.glob(f'{imdir}/aligned/*.fits')
        if self._commode == 'iraf':
            self.combine_iraf(alignedlist)
        elif self._commode == 'ccdproc':
            self.combine_ccdproc(alignedlist)
        
    #%%
    
imkey = '/data2/IMSNGDEEP/NGC3147/MAO_SNUCAM/*.fits'
stack  =Stacking(imkey)
all_tbl = stack.run()
